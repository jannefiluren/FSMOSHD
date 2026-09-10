# Stability correction for open/glacier surfaces.

struct NoStabilityCorrection{Tf} <: AbstractStabilityCorrection{Tf} end

@kwdef struct LouisStabilityCorrection{Tf} <: AbstractStabilityCorrection{Tf}
    bstb::Tf = 5                         # Atmospheric stability parameter (-)
end

NoStabilityCorrection{Tf}(grid::Grid; kwargs...) where {Tf} = NoStabilityCorrection{Tf}()
LouisStabilityCorrection{Tf}(grid::Grid; kwargs...) where {Tf} = LouisStabilityCorrection{Tf}(; kwargs...)

"""
    stability_factor(scheme, CD, z0, Ta, Tsrf, Ua, zU1, zT1)

Atmospheric stability correction applied to the open-terrain eddy diffusivity, following
Louis et al. (1982).
"""
function stability_factor end

@inline stability_factor(::NoStabilityCorrection{Tf}, CD, z0, Ta, Tsrf, Ua, zU1, zT1) where {Tf} = Tf(1)

@inline function stability_factor(sc::LouisStabilityCorrection{Tf}, CD, z0, Ta, Tsrf, Ua, zU1, zT1) where {Tf}
    @unpack_constants(Tf)
    RiB = grav * (Ta - Tsrf) * zU1^Tf(2) / (zT1 * Ta * Ua^Tf(2))
    if (RiB > Tf(0.2))
        RiB = Tf(0.2)
    end
    if (RiB > Tf(0))
        fh = Tf(1) / (Tf(1) + Tf(3) * sc.bstb * RiB * sqrt(Tf(1) + sc.bstb * RiB))
    else
        fh = Tf(1) - Tf(3) * sc.bstb * RiB / (Tf(1) + Tf(3) * sc.bstb^Tf(2) * CD * sqrt(-RiB * zU1 / z0))
    end
    return fh
end
