# Surface-layer structure.

struct OpenSurfaceLayer{Tf, S <: AbstractStabilityCorrection{Tf}} <: AbstractSurfaceLayer{Tf}
    # Atmospheric stability correction for the open-terrain eddy diffusivity. A type parameter
    # (so the scheme stays isbits); its default is on the constructor, since S follows the value.
    stability::S
end

@kwdef struct ForestSurfaceLayer{Tf} <: AbstractSurfaceLayer{Tf}
    rchd::Tf = 0.67                      # Ratio of displacement height to canopy height (-)
    rchz::Tf = 0.2                       # Ratio of roughness length to canopy height (-)
    zgf::Tf = 1                          # Roughness length adjustment factor vs vegetation fraction (-)
    zgr::Tf = 0                          # Roughness length adjustment range vs vegetation fraction (-)
    wcan::Tf = 2.5                       # Parameter of exponential wind profile (-)
    khcf::Tf = 3                         # Diffusivity adjustment for canopy effects (-)
    cveg::Tf = 20                        # Vegetation turbulent transfer coefficient ((s/m)^0.5)
    gsnf::Tf = 0                         # Snow-free vegetation moisture conductance (m/s)
    zsub::Tf = 2                         # Sub-canopy reference height (m)
end

# Default to the Louis (1982) correction, the realistic choice used across the model; pass
# stability = NoStabilityCorrection{Tf}() for a neutral surface layer.
OpenSurfaceLayer{Tf}(; stability = LouisStabilityCorrection{Tf}()) where {Tf} =
    OpenSurfaceLayer{Tf, typeof(stability)}(stability)
OpenSurfaceLayer{Tf}(grid::Grid; kwargs...) where {Tf} = OpenSurfaceLayer{Tf}(; kwargs...)
ForestSurfaceLayer{Tf}(grid::Grid; kwargs...) where {Tf} = ForestSurfaceLayer{Tf}(; kwargs...)

"""
    exchange_coefficients!(surface_layer, i, j, state, diag, surface, params, meteo, zU1, zT1, z0g)

Eddy diffusivities for turbulent heat and moisture transfer at cell `(i, j)`, implemented
for every `AbstractSurfaceLayer`: `diag.KH` and `diag.KWg` over open terrain, and
`diag.KHa`, `diag.KHg`, `diag.KHv`, `diag.KWg`, `diag.KWv` under a canopy. The wind and
temperature reference heights `zU1`, `zT1` and the ground roughness length `z0g` are
resolved by the caller.
"""
function exchange_coefficients! end

# Open/glacier terrain
@inline function exchange_coefficients!(sl::OpenSurfaceLayer{Tf}, i, j, state, diag, surface, params, meteo, zU1, zT1, z0g) where {Tf}
    @unpack_constants(Tf)
    (; Sice, Tsrf) = state
    (; KH, KWg, gs1, Qa, Uaeff) = diag
    (; Ta, Ps) = meteo

    # Roughness lengths and friction velocity
    z0 = z0g
    z0h = Tf(0.1) * z0
    CD = (vkman / log(zU1 / z0))^Tf(2)
    ustar = sqrt(CD) * Uaeff[i, j]

    fh = stability_factor(sl.stability, CD, z0, Ta[i, j], Tsrf[i, j], Uaeff[i, j], zU1, zT1)

    # Eddy diffusivities
    KH[i, j] = fh * vkman * ustar / log(zT1 / z0h)
    Qs = qsat(Ps[i, j], Tsrf[i, j])
    if (Sice[1, i, j] > eps(Tf) || Qa[i, j] > Qs)
        KWg[i, j] = KH[i, j]
    else
        KWg[i, j] = gs1[i, j] * KH[i, j] / (gs1[i, j] + KH[i, j])
    end
    return nothing
end

# Forest terrain
@inline function exchange_coefficients!(sl::ForestSurfaceLayer{Tf}, i, j, state, diag, surface, params, meteo, zU1, zT1, z0g) where {Tf}
    @unpack_constants(Tf)
    (; zU) = params
    (; zsub, gsnf) = sl
    (; fveg, fves, VAI, hcan) = surface
    (; Sveg, Tsrf, Tveg, Qcan) = state
    (; KHa, KHg, KHv, KWg, KWv, Usc, gs1, Uaeff) = diag
    (; Ps) = meteo

    # Roughness lengths, friction velocity and canopy wind profile
    z0g = (sl.zgf + sl.zgr * fveg[i, j]) * z0g
    z0h = Tf(0.1) * z0g
    dh = sl.rchd * hcan[i, j]
    z0v = sl.rchz * hcan[i, j]
    ustar = vkman * Uaeff[i, j] / log((zU1 - dh) / z0v)
    Uh = (ustar / vkman) * log((hcan[i, j] - dh) / z0v)
    KHh = vkman * ustar * (hcan[i, j] - dh)
    Usf = exp(sl.wcan * (zsub / hcan[i, j] - Tf(1))) * Uh

    Uso = Uaeff[i, j] * log(zsub / z0g) / log(zU / z0g)

    # Eddy diffusivities
    rad = (log((zT1 - dh) / (hcan[i, j] - dh)) / (vkman * ustar) + hcan[i, j] * (exp(sl.wcan * (Tf(1) - (z0v + dh) / hcan[i, j])) - Tf(1)) / (sl.wcan * KHh)) / sl.khcf
    KHa[i, j] = sqrt(fves[i, j]) / rad
    Usub = sqrt(fves[i, j]) * Usf + (Tf(1) - sqrt(fves[i, j])) * Uso
    Usub = max(Usub, Tf(0.1))
    rgd = Tf(1) / (vkman^Tf(2) * Usub) * log(zsub / z0h) * log(zsub / z0g)
    KHg[i, j] = Tf(1) / rgd
    Uc = exp(sl.wcan * ((z0v + dh) / hcan[i, j] - Tf(1))) * Uh
    KHv[i, j] = VAI[i, j] * sqrt(Uc) / sl.cveg
    # Usc leaves the model only through the OSHDinternal output catalog (uaca)
    Usc[i, j] = Usub

    Qs = qsat(Ps[i, j], Tsrf[i, j])
    if (Qcan[i, j] > Qs)
        KWg[i, j] = KHg[i, j]
    else
        KWg[i, j] = gs1[i, j] * KHg[i, j] / (gs1[i, j] + KHg[i, j])
    end
    Qs = qsat(Ps[i, j], Tveg[i, j])
    if (Sveg[i, j] > eps(Tf) || Qcan[i, j] > Qs)
        KWv[i, j] = KHv[i, j]
    else
        KWv[i, j] = gsnf * KHv[i, j] / (gsnf + KHv[i, j])
    end
    return nothing
end
