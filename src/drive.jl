"""
    drive!(fsm, meteo)

Meteorological data preprocessing and unit conversions.

Computes derived meteorological variables (saturation vapour pressure,
specific humidity, wind speed with lower bound, effective snowfall) and
stores them in the model state structure.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions (read-only)
"""
function drive!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    @unpack_constants(Tf)

    @unpack es, Qa, Ua_eff, Sf_eff = fsm

    @unpack Ua, Sf, Ta, RH, Ps = meteo

    Ua_eff .= max.(Ua, Tf(0.1))

    es .= e0 .* exp.(Tf(17.5043) .* (Ta .- Tm) ./ (Tf(241.3) .+ (Ta .- Tm)))
    Qa .= (RH ./ 100) .* eps_fsm .* es ./ Ps

    # Snowfall reaching the surface; further adjusted by canopy! for forest tiles
    Sf_eff .= Sf

    return nothing

end
