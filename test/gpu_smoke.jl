# GPU smoke test: runs the full physics step! on the CPU and on the GPU with
# identical synthetic forcing and compares all model state, field by field.
#
# This script is run manually (it is NOT part of runtests.jl):
#
#   1. Install CUDA once into your default Julia environment (it is not a
#      dependency of this package; environment stacking makes it loadable):
#
#        julia -e "using Pkg; Pkg.add(\"CUDA\")"
#
#   2. From the package directory:
#
#        julia --project=. test/gpu_smoke.jl
#
# Without a functional CUDA setup the script falls back to the
# KernelAbstractions CPU backend for the "device" run, so the script logic
# itself can be exercised on any machine (the comparison is then trivially
# bit-identical).
#
# Expected result on a GPU: every field agrees within the Float32 round-off
# tolerances already used for the Fortran comparisons (rtol 1e-4, atol 1e-5);
# Nsnow agrees exactly except possibly in rare borderline cells where a
# one-ulp difference flips a layering threshold (any such cell also shows up
# in the layer fields).
#
# The printed timings are informational only - at this grid size, and on
# bandwidth-limited cards, they say nothing about full-domain performance.

using FlexibleSnowModelOSHD
using Dates
using Printf
import KernelAbstractions

# Load CUDA from the stacked default environment if available
cuda_functional = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end

if cuda_functional
    CUDA.allowscalar(false)
    device_arch = GPU(CUDABackend())
    device_name = "GPU: " * CUDA.name(CUDA.device())
else
    device_arch = GPU(KernelAbstractions.CPU())
    device_name = "KernelAbstractions CPU backend (no functional CUDA found - script-logic test only)"
end

println("Device architecture: ", device_name)
println("Host threads:        ", Threads.nthreads())

# ---------------------------------------------------------------------------
# Synthetic case: terrain, three tile configurations, diurnal forcing
# ---------------------------------------------------------------------------

const Tf = Float32
const Ti = Int32
const Nx, Ny = 256, 256
const nsteps = 48
const t0 = DateTime(2026, 1, 15, 0)

function make_landuse()
    dem = [Tf(1000 + 1500 * (i - 1) / (Nx - 1) + 200 * sin(2π * j / Ny)) for i in 1:Nx, j in 1:Ny]
    d(x) = Dict("data" => x)
    return Dict(
        "elevation" => d(dem),
        "skyvf" => d(fill(0.95f0, Nx, Ny)),
        "prec_multi" => d(ones(Float64, Nx, Ny)),
        "slopemu" => d(fill(0.2f0, Nx, Ny)),
        "xi" => d(fill(150.0f0, Nx, Ny)),
        "Ld" => d(fill(250.0f0, Nx, Ny)),
        "forest" => d(fill(0.6f0, Nx, Ny)),
        "glacier" => d(fill(0.5f0, Nx, Ny)),
        "fveg" => d(fill(0.5f0, Nx, Ny)),
        "hcan" => d(fill(12.0f0, Nx, Ny)),
        "lai" => d(fill(2.5f0, Nx, Ny)),
        "vfhp" => d(fill(0.5f0, Nx, Ny)),
        "fves" => d(fill(0.5f0, Nx, Ny)),
    )
end

# Mirrors the configurations of the regression tests: exercises ebalsrf and
# ebalfor, both snow cover fraction code paths, and the glacier branches
const configs = [
    ("open", Dict("tile" => "open", "config" => Dict("SNFRAC" => 0))),
    (
        "forest", Dict(
            "tile" => "forest",
            "config" => Dict("CANMOD" => 1, "EXCHNG" => 2, "SNFRAC" => 4, "ZOFFST" => 1),
            "params" => Dict("hfsn" => 0.3, "z0_snow" => 0.01),
        ),
    ),
    ("glacier", Dict("tile" => "glacier", "config" => Dict("SNFRAC" => 0))),
]

# Initial snowpack, set on the CPU structure before moving it to the device
function init_snowpack!(fsm)
    for j in 1:Ny, i in 1:Nx
        fsm.state.Nsnow[i, j] = 2
        fsm.state.fsnow[i, j] = 1.0f0
        for (k, ds) in enumerate((0.1f0, 0.2f0))
            fsm.state.Ds[k, i, j] = ds
            fsm.state.Sice[k, i, j] = (150.0f0 + 10.0f0 * (i % 5)) * ds
            fsm.state.Sliq[k, i, j] = 0.0f0
            fsm.state.Tsnow[k, i, j] = 263.0f0 + k
            fsm.state.histowet[k, i, j] = 0.0f0
        end
    end
    return nothing
end

# Hourly forcing precomputed as host arrays (applied to the device with
# copyto!, which works across architectures)
function make_forcing(landuse)
    dem = landuse["elevation"]["data"]
    forcing = []
    for h in 0:(nsteps - 1)
        hod = mod(h, 24)
        sun = max(0.0f0, sin(π * (hod - 6) / 12))
        Ta = @. Tf(272.0 - 0.0065 * (dem - 1500) + 4 * sun + 0.5 * sin(2π * dem / 300))
        push!(
            forcing, Dict(
                :Ta => Ta,
                :RH => fill(Tf(75), Nx, Ny),
                :Ua => fill(Tf(3 + 2 * sun), Nx, Ny),
                :Ps => @.(Tf(101325 * exp(-dem / 8000))),
                :Sdir => fill(Tf(400 * sun), Nx, Ny),
                :Sdif => fill(Tf(80 * sun), Nx, Ny),
                :Sdird => fill(Tf(350 * sun), Nx, Ny),
                :LW => @.(Tf(250 + 2 * (Ta - 260))),
                :Sf => fill(Tf(hod < 6 ? 2.0e-4 : 0.0), Nx, Ny),   # snowfall at night
                :Rf => fill(Tf(hod == 14 ? 5.0e-5 : 0.0), Nx, Ny), # a little rain once a day
                :Sf24h => fill(Tf(4.3), Nx, Ny),
                :Tv => fill(Tf(0.5), Nx, Ny),
            )
        )
    end
    return forcing
end

function apply_forcing!(met, f)
    for (name, value) in f
        copyto!(getfield(met, name), value)
    end
    return nothing
end

function run_case(arch, landuse, settings, forcing)
    fsm = setup(FlexibleSnowModelOSHD.CPU(), Tf, Ti, landuse, Nx, Ny, settings)
    init_snowpack!(fsm)
    fsm = on_architecture(arch, fsm)
    met = on_architecture(arch, MET{Tf, Ti}(Nx = Nx, Ny = Ny))
    # Time only the second half of the steps: the first steps pay the
    # (config-dependent) kernel compilation on both CPU and GPU, which would
    # otherwise dominate the average
    ntimed = nsteps - nsteps ÷ 2
    elapsed = 0.0
    for h in 1:nsteps
        apply_forcing!(met, forcing[h])
        if h > nsteps ÷ 2
            elapsed += @elapsed step!(fsm, met, t0 + Hour(h - 1))
        else
            step!(fsm, met, t0 + Hour(h - 1))
        end
    end
    return on_architecture(FlexibleSnowModelOSHD.CPU(), fsm), elapsed / ntimed
end

# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

const compare_fields = [
    # (sub-struct, field) — state
    (:state, :Tsrf), (:state, :fsnow), (:state, :albs), (:state, :Ds), (:state, :Sice),
    (:state, :Sliq), (:state, :Tsnow), (:state, :Tsoil), (:state, :theta), (:state, :histowet),
    (:state, :Sveg), (:state, :Tveg), (:state, :Tcan), (:state, :Qcan),
    (:state, :swemin), (:state, :swemax), (:state, :swehist),
    (:state, :snowdepthmin), (:state, :snowdepthmax), (:state, :snowdepthhist),
    # fluxes and diagnostics
    (:diag, :H), (:diag, :LE), (:diag, :G), (:diag, :Rnet), (:diag, :Esrf), (:diag, :Eveg),
    (:diag, :Melt), (:diag, :Roff), (:diag, :meltflux_out), (:diag, :Sbsrf), (:diag, :Gsoil),
]

function compare_case(name, fsm_cpu, fsm_dev; rtol = 1.0f-4, atol = 1.0f-5)
    println("\n--- $name ---------------------------------------------------")
    ok = true

    ndiff = count(fsm_cpu.state.Nsnow .!= fsm_dev.state.Nsnow)
    @printf("  %-15s %s (%d differing cells)\n", "Nsnow", ndiff == 0 ? "PASS" : "FAIL", ndiff)
    ok &= ndiff == 0

    for (sub, field) in compare_fields
        a = getfield(getfield(fsm_cpu, sub), field)
        b = getfield(getfield(fsm_dev, sub), field)
        maxabs = maximum(abs.(a .- b))
        pass = isapprox(a, b; rtol = rtol, atol = atol)
        @printf("  %-15s %s (max abs diff %.3e)\n", field, pass ? "PASS" : "FAIL", maxabs)
        ok &= pass
    end
    return ok
end

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

landuse = make_landuse()
forcing = make_forcing(landuse)

all_ok = true
for (name, settings) in configs
    fsm_cpu, dt_cpu = run_case(FlexibleSnowModelOSHD.CPU(), landuse, settings, forcing)
    fsm_dev, dt_dev = run_case(device_arch, landuse, settings, forcing)
    global all_ok &= compare_case(name, fsm_cpu, fsm_dev)
    @printf(
        "  timing (informational): CPU %.1f ms/step (%d threads), device %.1f ms/step\n",
        1000 * dt_cpu, Threads.nthreads(), 1000 * dt_dev
    )
end

println()
if all_ok
    println("GPU SMOKE TEST PASSED ($device_name)")
else
    println("GPU SMOKE TEST FAILED ($device_name) - see FAIL lines above")
    exit(1)
end
