
# Reference height
struct AboveGround{Tf} <: AbstractReferenceHeight{Tf} end
struct AboveCanopy{Tf} <: AbstractReferenceHeight{Tf} end
AboveGround{Tf}(Nx, Ny; kwargs...) where {Tf} = AboveGround{Tf}()
AboveCanopy{Tf}(Nx, Ny; kwargs...) where {Tf} = AboveCanopy{Tf}()

@inline reference_heights(::AboveGround, zU, zT, hcan) = (zU, zT)
@inline reference_heights(::AboveCanopy, zU, zT, hcan) = (zU + hcan, zT + hcan)

# Stability correction for open/glacier surfaces
struct NoStabilityCorrection{Tf} <: AbstractStabilityCorrection{Tf} end
NoStabilityCorrection{Tf}(Nx, Ny; kwargs...) where {Tf} = NoStabilityCorrection{Tf}()

@kwdef struct LouisStabilityCorrection{Tf} <: AbstractStabilityCorrection{Tf}
    bstb::Tf = 5                         # Atmospheric stability parameter (-)
end
LouisStabilityCorrection{Tf}(Nx, Ny; kwargs...) where {Tf} = LouisStabilityCorrection{Tf}(; kwargs...)

@inline stability_factor(::NoStabilityCorrection{Tf}, CD, z0, Ta, Tsrf, Ua, zU1, zT1) where {Tf} = Tf(1)

# Stability correction factor following Louis et al. (1982)
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

# Surface-layer structure
struct OpenSurfaceLayer{Tf, S <: AbstractStabilityCorrection{Tf}} <: AbstractSurfaceLayer{Tf}
    stability::S
end
OpenSurfaceLayer{Tf}(; stability = NoStabilityCorrection{Tf}()) where {Tf} =
    OpenSurfaceLayer{Tf, typeof(stability)}(stability)
OpenSurfaceLayer{Tf}(Nx, Ny; stability = NoStabilityCorrection{Tf}()) where {Tf} =
    OpenSurfaceLayer{Tf}(; stability = stability)

@kwdef struct ForestSurfaceLayer{Tf} <: AbstractSurfaceLayer{Tf}
    rchd::Tf = 0.67                      # Ratio of displacement height to canopy height (-)
    rchz::Tf = 0.2                       # Ratio of roughness length to canopy height (-)
    zgf::Tf = 1                          # Roughness length adjustment factor vs vegetation fraction (-)
    zgr::Tf = 0                          # Roughness length adjustment range vs vegetation fraction (-)
    wcan::Tf = 2.5                       # Parameter of exponential wind profile (-)
    khcf::Tf = 3                         # Diffusivity adjustment for canopy effects (-)
    cveg::Tf = 20                        # Vegetation turbulent transfer coefficient ((s/m)^0.5)
end
ForestSurfaceLayer{Tf}(Nx, Ny; kwargs...) where {Tf} = ForestSurfaceLayer{Tf}(; kwargs...)

# Surface exchange coefficients for open/glacier terrain
@inline function exchange_coefficients!(sl::OpenSurfaceLayer{Tf}, state, diag, landuse, params, meteo, i, j, zU1, zT1, z0g) where {Tf}
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

# Surface exchange coefficients for forest terrain
@inline function exchange_coefficients!(sl::ForestSurfaceLayer{Tf}, state, diag, landuse, params, meteo, i, j, zU1, zT1, z0g) where {Tf}
    @unpack_constants(Tf)
    (; zU, zsub, gsnf) = params
    (; fveg, fves, VAI, hcan) = landuse
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

function sfexch!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    (; reference_height, surface_layer, SNFRAC) = fsm.physics

    backend = get_backend(fsm.state.Tsrf)
    kernel! = sfexch_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.landuse, fsm.params, meteo,
        reference_height, surface_layer, SNFRAC;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function sfexch_kernel!(
        state, diag, landuse, params::Parameters{Tf}, meteo,
        reference_height::AbstractReferenceHeight{Tf},
        surface_layer::AbstractSurfaceLayer{Tf},
        SNFRAC::AbstractSnowFraction{Tf},
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; tthresh, zT, zU) = params
    (; tilefrac, z0_snow, z0sf, hcan) = landuse
    (; Ds, fsnow) = state

    if (tilefrac[i, j] >= tthresh)

        zU1, zT1 = reference_heights(reference_height, zU, zT, hcan[i, j])

        # Ground roughness length
        z0g = z0_snow[i, j]
        if SNFRAC isa PointSnowFraction
            sumtmp = column_sum(Ds, i, j)
            if (sumtmp <= Tf(0.05))
                z0g = z0sf[i, j]
            end
        else
            if (fsnow[i, j] <= eps(Tf))
                z0g = z0sf[i, j]
            end
        end

        exchange_coefficients!(surface_layer, state, diag, landuse, params, meteo, i, j, zU1, zT1, z0g)

    end
end
