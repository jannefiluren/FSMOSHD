# Surface-exchange parameterizations. The former EXCHNG flag is two axes -
# surface-layer structure (Open/Forest) and stability correction (None/Louis) -
# and ZOFFST is the reference-height axis. Each is an isbits scheme dispatched by
# small @inline point functions the kernel calls per cell.

# --- ZOFFST: where reference heights are measured from --------------------
struct AboveGround{Tf} <: AbstractReferenceHeight{Tf} end
struct AboveCanopy{Tf} <: AbstractReferenceHeight{Tf} end
AboveGround{Tf}(Nx, Ny; kwargs...) where {Tf} = AboveGround{Tf}()
AboveCanopy{Tf}(Nx, Ny; kwargs...) where {Tf} = AboveCanopy{Tf}()

@inline reference_heights(::AboveGround, zU, zT, hcan) = (zU, zT)
@inline reference_heights(::AboveCanopy, zU, zT, hcan) = (zU + hcan, zT + hcan)

# --- Surface-layer structure (roughness lengths and friction velocity) ----
# Each returns its own bundle of intermediates; nothing downstream sees both.
@kwdef struct OpenSurfaceLayer{Tf} <: AbstractSurfaceLayer{Tf}
    rchd::Tf = 0.67                      # Ratio of displacement height to canopy height (-)
    rchz::Tf = 0.2                       # Ratio of roughness length to canopy height (-)
    cden::Tf = 0.004                     # Dense canopy turbulent transfer coefficient (-)
    cveg::Tf = 20                        # Vegetation turbulent transfer coefficient ((s/m)^0.5)
end

@kwdef struct ForestSurfaceLayer{Tf} <: AbstractSurfaceLayer{Tf}
    rchd::Tf = 0.67                      # Ratio of displacement height to canopy height (-)
    rchz::Tf = 0.2                       # Ratio of roughness length to canopy height (-)
    zgf::Tf = 1                          # Roughness length adjustment factor vs vegetation fraction (-)
    zgr::Tf = 0                          # Roughness length adjustment range vs vegetation fraction (-)
    wcan::Tf = 2.5                       # Parameter of exponential wind profile (-)
    khcf::Tf = 3                         # Diffusivity adjustment for canopy effects (-)
    cveg::Tf = 20                        # Vegetation turbulent transfer coefficient ((s/m)^0.5)
end

OpenSurfaceLayer{Tf}(Nx, Ny; kwargs...) where {Tf} = OpenSurfaceLayer{Tf}(; kwargs...)
ForestSurfaceLayer{Tf}(Nx, Ny; kwargs...) where {Tf} = ForestSurfaceLayer{Tf}(; kwargs...)

@inline function surface_layer_state(sl::OpenSurfaceLayer{Tf}, s, z0g, zU1, zT1) where {Tf}
    z0v = sl.rchz * s.hcan
    z0 = (z0v^s.fveg) * (z0g^(Tf(1) - s.fveg))
    z0h = Tf(0.1) * z0
    dh = s.fveg * sl.rchd * s.hcan
    CD = (s.vkman / log((zU1 - dh) / z0))^Tf(2)
    ustar = sqrt(CD) * s.Ua
    return (z0g = z0g, z0 = z0, z0h = z0h, dh = dh, CD = CD, ustar = ustar)
end

@inline function surface_layer_state(sl::ForestSurfaceLayer{Tf}, s, z0g, zU1, zT1) where {Tf}
    forested = s.fveg > eps(Tf)
    # canopy height is only meaningful where there is canopy; substitute a safe
    # value where fveg == 0 so the (unselected) forest intermediates stay finite
    hcan = ifelse(forested, s.hcan, Tf(1))

    z0gf = (sl.zgf + sl.zgr * s.fveg) * z0g
    dh = sl.rchd * hcan
    z0v = sl.rchz * hcan
    z0g = ifelse(forested, z0gf, z0g)
    z0h = Tf(0.1) * z0g
    ustar = ifelse(forested, s.vkman * s.Ua / log((zU1 - dh) / z0v),
                             s.vkman * s.Ua / log(s.zU / z0g))
    Uh = (ustar / s.vkman) * log((hcan - dh) / z0v)
    KHh = s.vkman * ustar * (hcan - dh)
    Usf = exp(sl.wcan * (s.zsub / hcan - Tf(1))) * Uh
    return (z0g = z0g, z0h = z0h, z0v = z0v, dh = dh, ustar = ustar, Uh = Uh, KHh = KHh, Usf = Usf)
end

# --- Stability correction. "Needs nothing" (neutral) is the default -------
struct NoStabilityCorrection{Tf} <: AbstractStabilityCorrection{Tf} end
NoStabilityCorrection{Tf}(Nx, Ny; kwargs...) where {Tf} = NoStabilityCorrection{Tf}()

@kwdef struct LouisStabilityCorrection{Tf} <: AbstractStabilityCorrection{Tf}
    bstb::Tf = 5                         # Atmospheric stability parameter (-)
end
LouisStabilityCorrection{Tf}(Nx, Ny; kwargs...) where {Tf} = LouisStabilityCorrection{Tf}(; kwargs...)

@inline stability_factor(::NoStabilityCorrection{Tf}, sl, S, s, zU1, zT1) where {Tf} = Tf(1)
@inline canopy_richardson(::NoStabilityCorrection{Tf}, sl, S, s) where {Tf} = Tf(0)

# Louis is defined only against OpenSurfaceLayer: it needs S.CD/S.z0, which the
# forest bundle does not carry, so Forest + Louis is a MethodError, not garbage.
@inline function stability_factor(sc::LouisStabilityCorrection{Tf}, ::OpenSurfaceLayer, S, s, zU1, zT1) where {Tf}
    Tint = s.fveg * s.Tveg + (Tf(1) - s.fveg) * s.Tsrf
    RiB = s.grav * (s.Ta - Tint) * (zU1 - S.dh)^Tf(2) / ((zT1 - S.dh) * s.Ta * s.Ua^Tf(2))
    if (RiB > Tf(0.2))
        RiB = Tf(0.2)
    end
    if (RiB > Tf(0))
        fh = Tf(1) / (Tf(1) + Tf(3) * sc.bstb * RiB * sqrt(Tf(1) + sc.bstb * RiB))
    else
        fh = Tf(1) - Tf(3) * sc.bstb * RiB / (Tf(1) + Tf(3) * sc.bstb^Tf(2) * S.CD * sqrt(-RiB * zU1 / S.z0))
    end
    return fh
end

@inline function canopy_richardson(::LouisStabilityCorrection{Tf}, ::OpenSurfaceLayer, S, s) where {Tf}
    Ric = s.grav * (s.Tcan - s.Tsrf) * s.hcan / (s.Tcan * S.ustar^Tf(2))
    return max(min(Ric, Tf(10)), Tf(0))
end

# --- Eddy diffusivities, dispatched on the surface layer ------------------
@inline function eddy_diffusivities(sl::OpenSurfaceLayer{Tf}, S, s, fh, Ric, Uso, zT1) where {Tf}
    KHa = fh * s.vkman * S.ustar / log((zT1 - S.dh) / S.z0)
    KHg = s.vkman * S.ustar * ((Tf(1) - s.fveg) * fh / log(S.z0 / S.z0h) + s.fveg * sl.cden / (Tf(1) + Tf(0.5) * Ric))
    KHv = sqrt(S.ustar) * s.VAI / sl.cveg
    return (KHa = KHa, KHg = KHg, KHv = KHv, Usc = Uso)
end

@inline function eddy_diffusivities(sl::ForestSurfaceLayer{Tf}, S, s, fh, Ric, Uso, zT1) where {Tf}
    rad = (log((zT1 - S.dh) / (s.hcan - S.dh)) / (s.vkman * S.ustar) + s.hcan * (exp(sl.wcan * (Tf(1) - (S.z0v + S.dh) / s.hcan)) - Tf(1)) / (sl.wcan * S.KHh)) / sl.khcf
    KHa = sqrt(s.fves) / rad
    Usub = sqrt(s.fves) * S.Usf + (Tf(1) - sqrt(s.fves)) * Uso
    Usub = max(Usub, Tf(0.1))
    rgd = Tf(1) / (s.vkman^Tf(2) * Usub) * log(s.zsub / S.z0h) * log(s.zsub / S.z0g)
    KHg = Tf(1) / rgd
    Uc = exp(sl.wcan * ((S.z0v + S.dh) / s.hcan - Tf(1))) * S.Uh
    KHv = s.VAI * sqrt(Uc) / sl.cveg
    return (KHa = KHa, KHg = KHg, KHv = KHv, Usc = Usub)
end

"""
    sfexch!(fsm, meteo)

Surface exchange coefficients and turbulent transfer calculations.

The per-cell physics lives in `sfexch_kernel!`, a KernelAbstractions kernel
launched over the whole grid (see `ebalsrf!` for the pattern).
"""
function sfexch!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    (; SNFRAC, tthresh, zT, zU, gsnf, zsub) = fsm.params
    (; Nx, Ny) = fsm.grid
    (; z0_snow, z0sf, VAI, fveg, fves, hcan, tilefrac) = fsm.landuse
    (; Qcan, fsnow, Sice, Sveg, Tcan, Tsrf, Tveg, Ds) = fsm.state
    (; KH, KHa, KHg, KHv, KWg, KWv, Usc, gs1, Qa, Uaeff) = fsm.diag
    (; reference_height, surface_layer, stability) = fsm.physics

    @unpack Ta, Ps = meteo

    backend = get_backend(Tsrf)
    kernel! = sfexch_kernel!(backend)
    kernel!(
        KH, KHa, KHg, KHv, KWg, KWv, Usc,
        z0_snow, z0sf, VAI, Qcan, fsnow, Sice, Sveg, Tcan, Tsrf, Tveg, Ds,
        fveg, fves, hcan, tilefrac, gs1, Qa, Uaeff, Ta, Ps,
        tthresh, zT, zU, gsnf, zsub,
        reference_height, surface_layer, stability, SNFRAC;
        ndrange = (Int(Nx), Int(Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function sfexch_kernel!(
        KH, KHa, KHg, KHv, KWg, KWv, Usc,
        z0_snow, z0sf, VAI, Qcan, fsnow,
        Sice, Sveg, Tcan, Tsrf, Tveg, Ds,
        fveg, fves, hcan, tilefrac, gs1,
        Qa, Uaeff, Ta, Ps,
        tthresh::Tf, zT::Tf, zU::Tf, gsnf::Tf, zsub::Tf,
        reference_height::AbstractReferenceHeight{Tf},
        surface_layer::AbstractSurfaceLayer{Tf},
        stability::AbstractStabilityCorrection{Tf},
        SNFRAC::Ti,
    ) where {Tf, Ti}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        zU1, zT1 = reference_heights(reference_height, zU, zT, hcan[i, j])

        # Ground roughness length (SNFRAC stays a runtime axis until its own stage)
        z0g = z0_snow[i, j]
        if (SNFRAC == 3)
            sumtmp = column_sum(Ds, i, j)
            if (sumtmp <= Tf(0.05))
                z0g = z0sf[i, j]
            end
        else
            if (fsnow[i, j] <= eps(Tf))
                z0g = z0sf[i, j]
            end
        end

        s = (fveg = fveg[i, j], fves = fves[i, j], VAI = VAI[i, j], hcan = hcan[i, j],
             Ua = Uaeff[i, j], Tsrf = Tsrf[i, j], Tveg = Tveg[i, j], Tcan = Tcan[i, j],
             Ta = Ta[i, j], zU = zU, zsub = zsub, vkman = vkman, grav = grav)

        S = surface_layer_state(surface_layer, s, z0g, zU1, zT1)
        fh = stability_factor(stability, surface_layer, S, s, zU1, zT1)
        Ric = canopy_richardson(stability, surface_layer, S, s)

        Uso = Uaeff[i, j] * log(zsub / S.z0g) / log(zU / S.z0g)

        if (fveg[i, j] == 0)
            KH[i, j] = fh * vkman * S.ustar / log(zT1 / S.z0h)
            Qs = qsat(Ps[i, j], Tsrf[i, j])
            if (Sice[1, i, j] > eps(Tf) || Qa[i, j] > Qs)
                KWg[i, j] = KH[i, j]
            else
                KWg[i, j] = gs1[i, j] * KH[i, j] / (gs1[i, j] + KH[i, j])
            end
            Usc[i, j] = Uso
        else
            K = eddy_diffusivities(surface_layer, S, s, fh, Ric, Uso, zT1)
            KHa[i, j] = K.KHa
            KHg[i, j] = K.KHg
            KHv[i, j] = K.KHv
            Usc[i, j] = K.Usc

            Qs = qsat(Ps[i, j], Tsrf[i, j])
            if (Qcan[i, j] > Qs)
                KWg[i, j] = K.KHg
            else
                KWg[i, j] = gs1[i, j] * K.KHg / (gs1[i, j] + K.KHg)
            end
            Qs = qsat(Ps[i, j], Tveg[i, j])
            if (Sveg[i, j] > eps(Tf) || Qcan[i, j] > Qs)
                KWv[i, j] = K.KHv
            else
                KWv[i, j] = gsnf * K.KHv / (gsnf + K.KHv)
            end
        end

    end
end
