"""
    qsat(P, T)

Saturation specific humidity calculation using Tetens formula.

Computes saturation specific humidity for given pressure and temperature using empirical 
relations for saturation vapor pressure over water (T > 0°C) and ice (T ≤ 0°C).

# Arguments
- `P`: Atmospheric pressure (Pa)
- `T`: Temperature (K)

# Returns
- Saturation specific humidity (kg/kg)
"""
function qsat(P::Tf, T::Tf) where {Tf <: Real}

    @unpack_constants(Tf)

    Tc = T - Tm
    if (Tc > Tf(0))
        es = e0 * exp(Tf(17.5043) * Tc / (Tf(241.3) + Tc))
    else
        es = e0 * exp(Tf(22.4422) * Tc / (Tf(272.186) + Tc))
    end
    Qs = eps_fsm * es / P

    return Qs

end
