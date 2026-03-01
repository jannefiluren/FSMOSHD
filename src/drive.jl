"""
    drive!(fsm, meteo)

Meteorological data preprocessing and unit conversions.

# Arguments
- `fsm::FSM`: Model state structure
- `meteo::MET`: Current meteorological conditions (modified in-place)
"""
function drive!(fsm::FSM, meteo::MET{Tf,Ti}) where {Tf<:Real,Ti<:Integer}

  @unpack_constants(Tf)

  @unpack dt, wind_scaling = fsm

  @unpack es, Qa, Ua, Sf, Rf, Ta, RH, Ps = meteo

  Ua .= wind_scaling .* Ua
  Ua .= max.(Ua, Tf(0.1))

  es .= e0 .* exp.(Tf(17.5043) .* (Ta .- Tm) ./ (Tf(241.3) .+ (Ta .- Tm)))
  Qa .= (RH ./ 100) .* eps_fsm .* es ./ Ps

end
