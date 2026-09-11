# Fresh snow density parameterizations.

struct FixedFreshSnowDensity{Tf} <: AbstractFreshSnowDensity{Tf} end
struct ClimateFreshSnowDensity{Tf} <: AbstractFreshSnowDensity{Tf} end
struct ElevationFreshSnowDensity{Tf} <: AbstractFreshSnowDensity{Tf} end

FixedFreshSnowDensity{Tf}(grid::Grid; kwargs...) where {Tf} = FixedFreshSnowDensity{Tf}()
ClimateFreshSnowDensity{Tf}(grid::Grid; kwargs...) where {Tf} = ClimateFreshSnowDensity{Tf}()
ElevationFreshSnowDensity{Tf}(grid::Grid; kwargs...) where {Tf} = ElevationFreshSnowDensity{Tf}()

"""
    fresh_snow_density(scheme, rho0, rhob, rhoc, rhof, rhos_min, Ta, Ua, dem)

Density of fresh snow (kg/m^3) for one cell (air temperature `Ta`, wind `Ua`,
elevation `dem`).
"""
function fresh_snow_density end

# Fixed fresh snow density
@inline fresh_snow_density(::FixedFreshSnowDensity, rho0, rhob, rhoc, rhof, rhos_min, Ta, Ua, dem) = rho0

# Climate-dependent fresh snow density
@inline function fresh_snow_density(::ClimateFreshSnowDensity{Tf}, rho0, rhob, rhoc, rhof, rhos_min, Ta, Ua, dem) where {Tf}
    @unpack_constants(Tf)
    return max(rhof + rhob * (Ta - Tm) + rhoc * Ua^Tf(0.5), rhos_min)
end

# Climate-dependent with elevation-dependent decompaction
@inline function fresh_snow_density(::ElevationFreshSnowDensity{Tf}, rho0, rhob, rhoc, rhof, rhos_min, Ta, Ua, dem) where {Tf}
    @unpack_constants(Tf)
    rhonew = rhof + rhob * (Ta - Tm) + rhoc * Ua^Tf(0.5)
    if (dem <= Tf(1000))
        t_decompaction = Tf(24.0)
    elseif (dem > Tf(4000))
        t_decompaction = Tf(0.0)
    else
        t_decompaction = Tf(24) + (dem - Tf(1000)) / (Tf(4000) - Tf(1000)) * (Tf(0) - Tf(24))
    end
    rhonew = Tf(300) + (rhonew - Tf(300)) * exp(t_decompaction / Tf(100))
    return max(rhonew, rhos_min)
end
