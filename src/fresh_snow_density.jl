"""
    fresh_snow_density(FSNRHO, rho0, rhob, rhoc, rhof, rhos_min, Ta, Ua, dem)

Density of fresh snow (kg/m^3) for one grid cell.

Pure scalar function (safe inside KernelAbstractions kernels); the caller
passes the `FSNRHO` configuration and the density parameters from `FSM`.

# Arguments
- `FSNRHO`: Fresh snow density configuration (0: fixed, 1: climate-dependent,
  2: climate-dependent with elevation-dependent decompaction)
- `rho0`, `rhob`, `rhoc`, `rhof`, `rhos_min`: Density parameters (see `FSM`)
- `Ta`: Air temperature (K)
- `Ua`: Wind speed (m/s)
- `dem`: Grid elevation (m)
"""
@inline function fresh_snow_density(
        FSNRHO::Integer, rho0::Tf, rhob::Tf, rhoc::Tf, rhof::Tf, rhos_min::Tf,
        Ta::Tf, Ua::Tf, dem::Tf
    ) where {Tf <: Real}

    @unpack_constants(Tf)

    if (FSNRHO == 0)
        # Fixed fresh snow density
        rhonew = rho0
    elseif (FSNRHO == 1)
        # Climate-dependent fresh snow density
        rhonew = max(rhof + rhob * (Ta - Tm) + rhoc * Ua^Tf(0.5), rhos_min)
    else # FSNRHO == 2
        # Climate-dependent with elevation-dependent decompaction
        rhonew = rhof + rhob * (Ta - Tm) + rhoc * Ua^Tf(0.5)
        if (dem <= Tf(1000))
            t_decompaction = Tf(24.0)
        elseif (dem > Tf(4000))
            t_decompaction = Tf(0.0)
        else
            t_decompaction = Tf(24) + (dem - Tf(1000)) / (Tf(4000) - Tf(1000)) * (Tf(0) - Tf(24))
        end
        rhonew = Tf(300) + (rhonew - Tf(300)) * exp(t_decompaction / Tf(100))
        rhonew = max(rhonew, rhos_min)
    end

    return rhonew

end

"""
    fresh_snow_density!(fsm, Ta, Ua, dem)

Convenience method taking the configuration and parameters from `fsm`.
"""
function fresh_snow_density!(fsm::FSM{Tf, Ti}, Ta, Ua, dem) where {Tf <: Real, Ti <: Integer}
    (; FSNRHO, rho0, rhob, rhoc, rhof, rhos_min) = fsm.params
    return fresh_snow_density(FSNRHO, rho0, rhob, rhoc, rhof, rhos_min, Tf(Ta), Tf(Ua), Tf(dem))
end
