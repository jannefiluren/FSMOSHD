# GPU test: run the full physics step! on the CPU and on the GPU with identical synthetic
# forcing and compare every model-state and diagnostic field. CUDA is not a dependency of this
# package, so this is skipped unless CUDA is functional in the default (shared) environment.
# To enable it locally, add CUDA to your default environment once and rerun the suite:
#
#     julia -e 'using Pkg; Pkg.add("CUDA")'

using Dates

# Everything is kept local to this function so nothing (e.g. Nx/Ny) leaks into the test module and
# clashes with other test files. `device_arch` is the GPU architecture to compare the CPU against.
function run_gpu_smoke(device_arch)
    host = FlexibleSnowModelOSHD.CPU()
    Tf = Float32
    Nx, Ny = 256, 256
    nsteps = 48
    t0 = DateTime(2026, 1, 15, 0)

    d(x) = Dict("data" => x)
    dem = [Tf(1000 + 1500 * (i - 1) / (Nx - 1) + 200 * sin(2π * j / Ny)) for i in 1:Nx, j in 1:Ny]
    landuse = Dict(
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

    # Mirrors the regression configs: ebalsrf/ebalfor, both snow-fraction paths, glacier branches.
    configs = [
        ("open", Dict("tile" => "open", "physics" => Dict("snow_fraction" => SeasonalSnowFraction))),
        (
            "forest", Dict(
                "tile" => "forest",
                "physics" => Dict("canopy" => OneLayerCanopy, "snow_fraction" => TanhSnowFraction{Tf}(; hfsn = 0.3), "reference_height" => AboveCanopy),
                "params" => Dict("z0_snow" => 0.01),
            ),
        ),
        ("glacier", Dict("tile" => "glacier", "physics" => Dict("snow_fraction" => SeasonalSnowFraction))),
    ]

    # Diurnal forcing precomputed as host arrays (applied to the device with copyto!, which works
    # across architectures).
    forcing = map(0:(nsteps - 1)) do h
        hod = mod(h, 24)
        sun = max(0.0f0, sin(π * (hod - 6) / 12))
        Ta = @. Tf(272.0 - 0.0065 * (dem - 1500) + 4 * sun + 0.5 * sin(2π * dem / 300))
        Dict(
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
    end

    # Set up on the CPU, seed a snowpack, move to `arch`, run all steps, bring the result back.
    function run_case(arch, settings)
        fsm = setup(host, Grid(Tf; Nx = Nx, Ny = Ny), landuse, settings)
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
        fsm = on_architecture(arch, fsm)
        met = on_architecture(arch, MET{Tf}(Nx = Nx, Ny = Ny))
        for h in 1:nsteps
            for (name, value) in forcing[h]
                copyto!(getfield(met, name), value)
            end
            step!(fsm, met, t0 + Hour(h - 1))
        end
        return on_architecture(host, fsm)
    end

    compare_fields = [
        (:state, :Tsrf), (:state, :fsnow), (:state, :albs), (:state, :Ds), (:state, :Sice),
        (:state, :Sliq), (:state, :Tsnow), (:state, :Tsoil), (:state, :theta), (:state, :histowet),
        (:state, :Sveg), (:state, :Tveg), (:state, :Tcan), (:state, :Qcan),
        (:state, :swemin), (:state, :swemax), (:state, :swehist),
        (:state, :snowdepthmin), (:state, :snowdepthmax), (:state, :snowdepthhist),
        (:diag, :H), (:diag, :LE), (:diag, :G), (:diag, :Rnet), (:diag, :Esrf), (:diag, :Eveg),
        (:diag, :Melt), (:diag, :Roff), (:diag, :meltflux_out), (:diag, :Sbsrf), (:diag, :Gsoil),
    ]

    for (name, settings) in configs
        fsm_cpu = run_case(host, settings)
        fsm_dev = run_case(device_arch, settings)
        @testset "$name" begin
            # Layer count must match exactly; the float fields within Float32 round-off (a one-ulp
            # difference that flips a layering threshold would also show up in the layer fields).
            @test fsm_cpu.state.Nsnow == fsm_dev.state.Nsnow
            for (sub, field) in compare_fields
                a = getfield(getfield(fsm_cpu, sub), field)
                b = getfield(getfield(fsm_dev, sub), field)
                @test isapprox(a, b; rtol = 1.0f-4, atol = 1.0f-5)
            end
        end
    end
    return nothing
end

# CUDA is not a dependency of this package. `Pkg.test` runs in an isolated environment that does
# not stack the default (shared) `@v#.#` env, so add it to LOAD_PATH to make an already-installed
# CUDA loadable (this restores the load path a standalone `julia --project=.` would have). On a
# machine without CUDA the `using` fails and the test is skipped.
"@v#.#" in LOAD_PATH || push!(LOAD_PATH, "@v#.#")

# @eval so the `using` can sit inside the try (a bare `using` is not allowed in a control-flow block).
_cuda_functional = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end

if _cuda_functional
    CUDA.allowscalar(false)
    @info "GPU smoke test: running on " * CUDA.name(CUDA.device())
    @testset "GPU smoke" begin
        run_gpu_smoke(GPU(CUDABackend()))
    end
else
    @info "GPU smoke test skipped: no functional CUDA in the default environment"
end
