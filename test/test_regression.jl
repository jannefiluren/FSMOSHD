# Regression tests comparing current simulation results against reference data

using FlexibleSnowModelOSHD
using Test
using NCDatasets
using Dates
using Serialization

function load_domain_data()

    variables = [
        "easting", "northing", "elevation", "Ld", "dhdxdy", "sd", "prec_multi", "skyvf",
        "fveg", "hcan", "lai", "vfhp", "fves", "forest", "glacier",
    ]

    # Read domain data
    landuse = Dict()
    NCDataset(joinpath(pkgdir(FlexibleSnowModelOSHD), "data/domain_data.nc")) do ds
        for variable in variables
            landuse[variable] = Dict(
                "data" => ds[variable][:, :]
            )
        end
    end

    # Compute derived fields using the data arrays
    landuse["slopemu"] = Dict("data" => sqrt.((landuse["dhdxdy"]["data"] ./ 2)))
    landuse["xi"] = Dict("data" => (sqrt(2) * landuse["sd"]["data"]) ./ landuse["slopemu"]["data"])

    return landuse

end


function interpolate_meteo(Tf, landuse)

    var_mapping = [
        ("Sdir", 0.0, 1000.0),
        ("Sdif", 0.0, 1000.0),
        ("Sdird", 0.0, 1000.0),
        ("LW", 150.0, 500.0),
        ("Sf", 0.0, 100.0),
        ("Rf", 0.0, 100.0),
        ("Ta", 273.15 - 40, 273.15 + 40),
        ("RH", 20.0, 100.0),
        ("Ua", 0.05, 30.0),
        ("Ps", 0.0, 150000.0),
    ]

    times, data_meteo = NCDataset(joinpath(pkgdir(FlexibleSnowModelOSHD), "data/meteo_data.nc")) do ds
        z_stat = ds["elevation"][:]
        dz = landuse["elevation"]["data"] .- z_stat[1]
        times_local = ds["time"][:]

        Nx, Ny = size(landuse["elevation"]["data"])
        Nt = length(times_local)

        data_meteo_local = Dict()

        for (var_name, lower_bound, upper_bound) in var_mapping

            data_grid = Array{Tf}(undef, Nx, Ny, Nt)
            data_stat = ds[var_name][:, :]

            for i in eachindex(times_local)
                lapse_rate = (data_stat[2, i] - data_stat[1, i]) / (z_stat[2] - z_stat[1])
                data_grid[:, :, i] = data_stat[1, i] .+ lapse_rate .* dz
            end

            data_grid[data_grid .< lower_bound] .= lower_bound
            data_grid[data_grid .> upper_bound] .= upper_bound

            data_meteo_local[var_name] = data_grid

        end

        times_local, data_meteo_local
    end

    return times, data_meteo

end

function run_simulations(settings, Tf = Float32, Ti = Int32)

    # Read landuse data
    landuse = load_domain_data()

    # Generate meteo data
    times, data_meteo = interpolate_meteo(Tf, landuse)

    # Pre-scale wind speed (replaces former wind_scaling=0.7 parameter)
    data_meteo["Ua"] .= Tf(0.7) .* data_meteo["Ua"]

    # Setup
    Nx = size(landuse["elevation"]["data"], 1)
    Ny = size(landuse["elevation"]["data"], 2)
    Nt = length(times)

    fsm = setup(Tf, Ti, landuse, Nx, Ny, settings)
    met = MET{Tf, Ti}(Nx = Nx, Ny = Ny)

    # Preallocate arrays to store simulation results
    simulation_results = Dict{String, Any}()
    simulation_results["timestamps"] = Vector{String}(undef, Nt)

    state_vars = ["Tsrf", "Tsnow", "Sice", "Sliq", "fsnow", "albs", "Sveg", "Tveg", "Tcan"]
    flux_vars = ["SWsrf", "H", "LE", "G", "Melt", "Esrf", "Eveg", "Rnet"]
    diag_vars = ["Ds", "snowdepth"]

    for var in [state_vars; flux_vars; diag_vars]
        simulation_results[var] = Array{Tf, 3}(undef, Nx, Ny, Nt)
    end

    # Run simulation
    for (timestep, t) in enumerate(times)

        # Retrieve forcing data
        met.Sdir[:, :] = data_meteo["Sdir"][:, :, timestep]
        met.Sdif[:, :] = data_meteo["Sdif"][:, :, timestep]
        met.Sdird[:, :] = data_meteo["Sdird"][:, :, timestep]
        met.LW[:, :] = data_meteo["LW"][:, :, timestep]
        met.Sf[:, :] = data_meteo["Sf"][:, :, timestep] ./ fsm.dt  # Convert accumulation to rate (kg/m^2/s)
        met.Rf[:, :] = data_meteo["Rf"][:, :, timestep] ./ fsm.dt  # Convert accumulation to rate (kg/m^2/s)
        met.Ta[:, :] = data_meteo["Ta"][:, :, timestep]
        met.RH[:, :] = data_meteo["RH"][:, :, timestep]
        met.Ua[:, :] = data_meteo["Ua"][:, :, timestep]
        met.Ps[:, :] = data_meteo["Ps"][:, :, timestep]

        # Use constant value for time-varying transmissivity
        met.Tv[:, :] .= 0.5

        # Update snowfall tracking
        curr_hour = Dates.value(Hour(t)) + 1
        met.Sf24h_f64 .+= met.Sf .* fsm.dt  # Convert rate back to accumulation for 24h sum
        met.Sf24h_f64 .-= met.Sf_history_f64[:, :, curr_hour]
        met.Sf_history_f64[:, :, curr_hour] .= met.Sf .* fsm.dt
        met.Sf24h[:, :] .= met.Sf24h_f64

        # Run model
        step!(fsm, met, t)

        # Store results
        simulation_results["timestamps"][timestep] = string(t)

        # Store state variables
        simulation_results["Tsrf"][:, :, timestep] .= fsm.Tsrf
        simulation_results["Tsnow"][:, :, timestep] .= fsm.Tsnow[1, :, :]  # First snow layer
        simulation_results["Sice"][:, :, timestep] .= dropdims(sum(fsm.Sice, dims = 1), dims = 1)  # Total snow ice
        simulation_results["Sliq"][:, :, timestep] .= dropdims(sum(fsm.Sliq, dims = 1), dims = 1)  # Total snow liquid
        simulation_results["fsnow"][:, :, timestep] .= fsm.fsnow
        simulation_results["albs"][:, :, timestep] .= fsm.albs
        simulation_results["Sveg"][:, :, timestep] .= fsm.Sveg
        simulation_results["Tveg"][:, :, timestep] .= fsm.Tveg
        simulation_results["Tcan"][:, :, timestep] .= fsm.Tcan

        # Store flux variables
        simulation_results["SWsrf"][:, :, timestep] .= fsm.SWsrf
        simulation_results["H"][:, :, timestep] .= fsm.H
        simulation_results["LE"][:, :, timestep] .= fsm.LE
        simulation_results["G"][:, :, timestep] .= fsm.G
        simulation_results["Melt"][:, :, timestep] .= fsm.Melt
        simulation_results["Esrf"][:, :, timestep] .= fsm.Esrf
        simulation_results["Eveg"][:, :, timestep] .= fsm.Eveg
        simulation_results["Rnet"][:, :, timestep] .= fsm.Rnet

        # Store diagnostic variables
        simulation_results["Ds"][:, :, timestep] .= dropdims(sum(fsm.Ds, dims = 1), dims = 1)  # Total snow depth
        simulation_results["snowdepth"][:, :, timestep] .= simulation_results["Ds"][:, :, timestep] .* simulation_results["fsnow"][:, :, timestep]

    end

    return simulation_results

end

# Test data paths
projdir = dirname(dirname(@__FILE__))
ref_file = joinpath(projdir, "test", "simulation_results.jls")

# Configuration matrix
settings = [
    Dict(
        "tile" => "open",
        "config" => Dict("SNFRAC" => 0),
    ),
    Dict(
        "tile" => "forest",
        "config" => Dict("EXCHNG" => 2, "SNFRAC" => 4, "ZOFFST" => 1),
        "params" => Dict("hfsn" => 0.3, "z0_snow" => 0.01)
    ),
    Dict(
        "tile" => "glacier",
        "config" => Dict("SNFRAC" => 0),
    ),
]

# Load or create reference data
if !isfile(ref_file)
    simulation_refs = []
    for setting in settings
        push!(simulation_refs, run_simulations(setting))
    end
    serialize(ref_file, simulation_refs)
end

simulation_refs = deserialize(ref_file)

# Run regression tests for each configuration
tolerance = 1.0e-10

for (setting, simulation_ref) in zip(settings, simulation_refs)

    tile = setting["tile"]

    @testset "Tile: $tile" begin
        simulation_tst = run_simulations(setting)

        # Get list of variables to compare (exclude metadata)
        variables_to_compare = filter(k -> !(k in ["timestamps", "filenames"]), keys(simulation_tst))

        for variable in variables_to_compare
            @testset "$variable" begin
                @test size(simulation_tst[variable]) == size(simulation_ref[variable])
                max_diff = maximum(abs.(simulation_tst[variable] - simulation_ref[variable]))
                @test max_diff <= tolerance
            end
        end
    end

end
