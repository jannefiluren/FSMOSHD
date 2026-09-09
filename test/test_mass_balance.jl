using Dates
using CSV
using DataFrames
using FlexibleSnowModelOSHD
using Test

const path = dirname(@__FILE__)

function setup_open_example(SNFRAC)

    # set landuse properties
    lus = Dict()
    lus["skyvf"] = Dict("data" => [1.0;;])
    lus["elevation"] = Dict("data" => [2540.0;;])
    lus["slopemu"] = Dict("data" => [1.0;;])
    lus["xi"] = Dict("data" => [1.0;;])
    lus["Ld"] = Dict("data" => [1.0;;])
    lus["prec_multi"] = Dict("data" => [1.0;;])

    # define custom settings
    settings = Dict("tile" => "open", "config" => Dict("SNFRAC" => SNFRAC))

    # create fsm struct
    fsm = setup(Float32, lus, 1, 1, settings)

    # define meteo data struct
    met = MET{Float32}()

    # read meteo file
    df_meteo = CSV.read(joinpath(path, "../data/input_SLF_5WJ.txt"), DataFrame)

    return fsm, met, df_meteo

end

function setup_forest_example(SNFRAC)

    # set landuse properties
    lus = Dict()
    lus["skyvf"] = Dict("data" => [1.0;;])
    lus["elevation"] = Dict("data" => [2540.0;;])
    lus["slopemu"] = Dict("data" => [1.0;;])
    lus["xi"] = Dict("data" => [1.0;;])
    lus["Ld"] = Dict("data" => [1.0;;])
    lus["prec_multi"] = Dict("data" => [1.0;;])

    # add forest properties
    lus["forest"] = Dict("data" => [1;;]) # Forest cover fraction
    lus["fveg"] = Dict("data" => [0.6;;]) # Canopy cover fraction
    lus["fves"] = Dict("data" => [0.6;;]) # Stand-scale canopy cover fraction
    lus["hcan"] = Dict("data" => [20;;])  # Canopy height (m)
    lus["lai"] = Dict("data" => [2.5;;])  # Leaf area index
    lus["vfhp"] = Dict("data" => [0.5;;]) # Hemispherical sky-view fraction including canopy

    # define custom settings
    settings = Dict(
        "tile" => "forest",
        "config" => Dict(
            "EXCHNG" => 2,
            "ZOFFST" => 1,
            "SNFRAC" => SNFRAC,
            # No preferential deposition in canopy gaps
            "CANMOD" => OneLayerCanopy{Float32}(psr = 0, psf = 1),
        ),
    )

    # create fsm struct
    fsm = setup(Float32, lus, 1, 1, settings)

    # define meteo data struct
    met = MET{Float32}()

    # read meteo file
    df_meteo = CSV.read(joinpath(path, "../data/input_SLF_5WJ.txt"), DataFrame)

    return fsm, met, df_meteo

end

function run_fsm(fsm, met, df_meteo)

    # allocate output variable-wise
    prec = zeros(nrow(df_meteo))
    Roff = zeros(nrow(df_meteo))
    Sbsrf = zeros(nrow(df_meteo))
    Sbveg = zeros(nrow(df_meteo))

    # change in storage
    dSWE = -sum(fsm.state.Sice[:, 1, 1] .+ fsm.state.Sliq[:, 1, 1])
    dSveg = -fsm.state.Sveg[1, 1]

    # time loop
    for (i, row) in zip(1:nrow(df_meteo), eachrow(df_meteo))

        # record precipitation before model run as the forest tile modifies these fluxes
        prec[i] = row["Sf"] + row["Rf"]

        # assign input
        met.Sdir .= row["Sdir"]
        met.Sdif .= row["Sdif"]
        met.Sdird .= row["Sdir"]
        met.LW .= row["LW"]
        met.Sf .= row["Sf"] / fsm.params.dt  # Convert accumulation (kg/m^2) to rate (kg/m^2/s)
        met.Rf .= row["Rf"] / fsm.params.dt  # Convert accumulation (kg/m^2) to rate (kg/m^2/s)
        met.Ta .= row["Ta"]
        met.RH .= row["RH"]
        met.Ua .= row["Ua"]
        met.Ps .= row["Ps"]
        met.Sf24h .= row["Sf24h"]
        met.Tv .= 1

        # set time
        t = DateTime(row["year"], row["month"], row["day"], row["hour"])

        # run model and update states
        step!(fsm, met, t)

        # record mass fluxes
        Roff[i] = fsm.diag.Roff[1, 1]
        Sbsrf[i] = fsm.diag.Sbsrf[1, 1]
        Sbveg[i] = fsm.diag.Sbveg[1, 1]

    end

    # change in storage
    dSWE += sum(fsm.state.Sice[:, 1, 1] .+ fsm.state.Sliq[:, 1, 1])
    dSveg += fsm.state.Sveg[1, 1]

    return (
        prec = sum(prec),
        Roff = sum(Roff),
        Sbsrf = sum(Sbsrf),
        Sbveg = sum(Sbveg),
        dSWE = dSWE,
        dSveg = dSveg,
    )

end

function test_results(results, verbose = false)

    mass_actual_change = results.prec - results.Roff - results.Sbsrf - results.Sbveg
    mass_expected_change = results.dSWE + results.dSveg

    if verbose
        println("prec = " * string(results.prec))
        println("Roff = " * string(results.Roff))
        println("Sbsrf = " * string(results.Sbsrf))
        println("Sbveg = " * string(results.Sbveg))
        println("dSWE = " * string(results.dSWE))
        println("dSveg = " * string(results.dSveg))
        println("residual = ", mass_actual_change - mass_expected_change)
    end

    # Tolerance scaled to total precipitation
    @test isapprox(mass_actual_change, mass_expected_change, atol = 1.0e-3 * results.prec)
    return nothing

end

@testset "Mass balance tests" begin

    # Test open tile for all snow cover fraction schemes
    for SNFRAC in 0:4
        fsm, met, df_meteo = setup_open_example(SNFRAC)
        results = run_fsm(fsm, met, df_meteo)
        test_results(results)
    end

    # Test forest tile for all snow cover fraction schemes
    for SNFRAC in 0:4
        fsm, met, df_meteo = setup_forest_example(SNFRAC)
        results = run_fsm(fsm, met, df_meteo)
        test_results(results)
    end

end
