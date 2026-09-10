using Dates
using CSV
using DataFrames
using FlexibleSnowModelOSHD

const path = dirname(@__FILE__)

function setup_example()

    # set landuse properties
    lus = Dict()
    lus["skyvf"] = Dict("data" => [1.0;;])
    lus["elevation"] = Dict("data" => [2540.0;;])
    lus["slopemu"] = Dict("data" => [1.0;;])
    lus["xi"] = Dict("data" => [1.0;;])
    lus["Ld"] = Dict("data" => [1.0;;])
    lus["prec_multi"] = Dict("data" => [1.0;;])

    # define custom settings
    settings = Dict("tile" => "open")

    # create fsm struct
    fsm = setup(Grid(Float32; Nx = 1, Ny = 1), lus, settings)

    # define meteo data struct
    met = MET{Float32}()

    # read meteo file
    df_meteo = CSV.read(joinpath(path, "../data/input_SLF_5WJ.txt"), DataFrame)

    return fsm, met, df_meteo

end


function run_fsm(fsm, met, df_meteo)

    # allocate output variable-wise
    hs = zeros(nrow(df_meteo))

    # time loop
    for (i, row) in zip(1:nrow(df_meteo), eachrow(df_meteo))

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

        # set time
        t = DateTime(row["year"], row["month"], row["day"], row["hour"])

        # run model and update states
        step!(fsm, met, t)

        # write output
        hs[i] = dropdims(sum(fsm.state.Ds, dims = 1), dims = 1)[1, 1]

    end

    # write results to dataframe
    time = DateTime.(df_meteo[!, "year"], df_meteo[!, "month"], df_meteo[!, "day"], df_meteo[!, "hour"])
    df_results = DataFrame(time = time, hs = hs)

    return df_results

end

fsm, met, df_meteo = setup_example()

df_results = run_fsm(fsm, met, df_meteo)

describe(df_results)
