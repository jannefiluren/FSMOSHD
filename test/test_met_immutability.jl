# Tests that step! never modifies the meteorological input structure,
# so the same MET can safely be shared between tiles and threads
# (https://github.com/jannefiluren/FlexibleSnowModelOSHD.jl/issues/15)

using Dates
using FlexibleSnowModelOSHD
using Test

function setup_immutability_example(tile)

    lus = Dict()
    lus["skyvf"] = Dict("data" => [0.9;;])
    lus["elevation"] = Dict("data" => [2540.0;;])
    lus["slopemu"] = Dict("data" => [1.0;;])
    lus["xi"] = Dict("data" => [1.0;;])
    lus["Ld"] = Dict("data" => [1.0;;])
    lus["prec_multi"] = Dict("data" => [1.0;;])

    if tile == "forest"
        lus["forest"] = Dict("data" => [1;;])
        lus["fveg"] = Dict("data" => [0.6;;])
        lus["fves"] = Dict("data" => [0.6;;])
        lus["hcan"] = Dict("data" => [20;;])
        lus["lai"] = Dict("data" => [2.5;;])
        lus["vfhp"] = Dict("data" => [0.5;;])
        settings = Dict("tile" => "forest", "config" => Dict("CANMOD" => 1, "EXCHNG" => 2, "ZOFFST" => 1))
    else
        if tile == "glacier"
            lus["glacier"] = Dict("data" => [1;;])
        end
        settings = Dict("tile" => tile)
    end

    fsm = setup(Float32, lus, 1, 1, settings)
    met = MET{Float32}()

    # Wind speed below the 0.1 m/s minimum to exercise the clamping in drive!,
    # snowfall with Sf24h above Sfmin to exercise the albedo refresh in radiation!
    met.Sdir .= 400
    met.Sdif .= 100
    met.Sdird .= 350
    met.LW .= 280
    met.Sf .= 5.0f-4
    met.Rf .= 1.0f-4
    met.Ta .= 271
    met.RH .= 85
    met.Ua .= 0.05
    met.Ps .= 75000
    met.Sf24h .= 12
    met.Tv .= 0.5

    return fsm, met

end

@testset "MET not modified by step! ($tile tile)" for tile in ["open", "forest", "glacier"]

    fsm, met = setup_immutability_example(tile)

    met_ref = deepcopy(met)

    for t in DateTime(2023, 12, 1, 0):Hour(1):DateTime(2023, 12, 1, 5)
        step!(fsm, met, t)
    end

    for field in fieldnames(typeof(met))
        # isequal treats NaN entries (unset fields) as equal
        @test isequal(getfield(met, field), getfield(met_ref, field))
    end

end
