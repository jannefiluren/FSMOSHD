# Synthetic case definitions for test/gpu_smoke.jl: terrain, three tile
# configurations, initial snowpack and diurnal forcing.

using Dates

const Tf = Float32
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
    ("open", Dict("tile" => "open", "physics" => Dict("snow_fraction" => SeasonalSnowFraction))),
    (
        "forest", Dict(
            "tile" => "forest",
            "physics" => Dict("canopy" => OneLayerCanopy, "snow_fraction" => TanhSnowFraction, "reference_height" => AboveCanopy),
            "params" => Dict("hfsn" => 0.3, "z0_snow" => 0.01),
        ),
    ),
    ("glacier", Dict("tile" => "glacier", "physics" => Dict("snow_fraction" => SeasonalSnowFraction))),
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
