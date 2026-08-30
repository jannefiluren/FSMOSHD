# Snow albedo parameterizations. All schemes hold only scalars (isbits); the
# per-cell afs/adc maps live on Landuse and are passed in as arrays.

@kwdef struct DiagnosticAlbedo{Tf} <: AbstractAlbedo{Tf}
    amin::Tf = 0.6                       # Minimum albedo for melting snow (-)
    amax::Tf = 0.86                      # Maximum albedo for fresh snow (-)
    Talb::Tf = -2                        # Albedo decay temperature threshold (C)
end

@kwdef struct DecayAlbedo{Tf} <: AbstractAlbedo{Tf}
    amin::Tf = 0.6                       # Minimum albedo for melting snow (-)
    tcld::Tf = 3600 * 1000               # Cold snow albedo decay time scale (s)
    tmlt::Tf = 3600 * 100                # Melting snow albedo decay time scale (s)
    adfs::Tf = 3                         # Albedo adjustment, shortwave (-)
    adfl::Tf = 2                         # Albedo adjustment, longwave (-)
    Sfmin::Tf = 10                       # Min 24h snowfall to refresh albedo (kg/m^2)
end

@kwdef struct PrognosticAlbedo{Tf} <: AbstractAlbedo{Tf}
    ALRADT::Bool = true                  # Aspect-dependent decay tuning
    adm::Tf = 100                        # Melting snow albedo decay time (h)
    amin::Tf = 0.6                       # Minimum albedo for melting snow (-)
    Sfmin::Tf = 10                       # Min 24h snowfall to refresh albedo (kg/m^2)
end

DiagnosticAlbedo{Tf}(Nx, Ny; kwargs...) where {Tf} = DiagnosticAlbedo{Tf}(; kwargs...)
DecayAlbedo{Tf}(Nx, Ny; kwargs...) where {Tf} = DecayAlbedo{Tf}(; kwargs...)
PrognosticAlbedo{Tf}(Nx, Ny; kwargs...) where {Tf} = PrognosticAlbedo{Tf}(; kwargs...)

"""
    snow_albedo(scheme, albs, s)

Return the updated snow albedo for one cell, given the current albedo `albs` and
a NamedTuple `s` of that cell's scalar drivers. A pure function; the kernel does
the `[i, j]` extraction and writes the result back. Every `AbstractAlbedo`
implements it.
"""
function snow_albedo end

@inline function snow_albedo(c::DiagnosticAlbedo, albs, s)
    afs_loc = c.amax
    a = c.amin + (afs_loc - c.amin) * (s.Tsrf - s.Tm) / c.Talb
    a = max(a, min(afs_loc, c.amin))
    a = min(a, max(afs_loc, c.amin))
    return a
end

@inline function snow_albedo(c::DecayAlbedo{Tf}, albs, s) where {Tf}
    (; Tsrf, fveg, trcn, fsky, afs, Tm, dt, summer_decay, Sdir, Sdif, Sf, Tv) = s

    tau = c.tcld
    if (Tsrf >= Tm)
        tau = c.tmlt
    end
    # Forest adjustments -> not yet properly tested for OSHD but option currently unused
    if summer_decay
        tau = Tf(70.0) * Tf(3600.0)
    end

    if fveg > Tf(0) && Sdir > eps(Tf)
        tau = tau / ((Tf(1) - trcn * fsky) * (Tf(1) + c.adfl * Tv) + c.adfs * Tv)
    elseif fveg > Tf(0) && Sdif > eps(Tf)
        tau = tau / ((Tf(1) - trcn * fsky) + c.adfs * trcn * fsky)
    elseif (fveg > Tf(0) && (Sdir + Sdif <= eps(Tf)))
        tau = tau / (Tf(2.0) - trcn * fsky)
    end

    rt = Tf(1) / tau + Sf / c.Sfmin
    alim = (c.amin / tau + Sf * afs / c.Sfmin) / rt
    a = alim + (albs - alim) * exp(-rt * dt)
    if (a < min(afs, c.amin))
        a = min(afs, c.amin)
    end
    if (a > max(afs, c.amin))
        a = max(afs, c.amin)
    end
    return a
end

@inline function snow_albedo(c::PrognosticAlbedo{Tf}, albs, s) where {Tf}
    (; Tsrf, afs, adc, SWEtmp, Tm, dt, Sdir, Sdird, Sf, Sf24h) = s

    adc_loc = adc
    adm_loc = c.adm
    afs_loc = afs

    # BC 08.23: aspect-dependent albedo tuning. Activated for oper season 2024 or optionally.
    # BC Oct 23: Jan's suggestion: modify only when the decay rate should be increased
    # (ad* DECREASE), not decreased
    if c.ALRADT
        if ((Sdir > eps(Tf)) && (Sdird < Sdir))
            adm_loc = adm_loc * (Sdird) / (Sdir)
            adc_loc = adc_loc * (Sdird) / (Sdir)
            if (adm_loc < eps(Tf))
                adm_loc = eps(Tf)
            end
            if (adc_loc < eps(Tf))
                adc_loc = eps(Tf)
            end
        end
    end

    a = albs
    if (Tsrf >= Tm)
        a = (a - c.amin) * exp(-(dt / Tf(3600)) / adm_loc) + c.amin
    else
        a = a - (dt / Tf(3600)) / adc_loc
    end
    if (SWEtmp < Tf(75.0)) # more stuff showing on and up through snow
        afs_loc *= Tf(0.8)
    end
    # Reset to fresh snow albedo (wasn't originally available; only else term)
    if ((Sf * dt) > Tf(0.0) && Sf24h > c.Sfmin)
        a = afs_loc
    else
        a = a + (afs_loc - a) * Sf * dt / c.Sfmin
    end
    ## End Adjustments
    if (a > afs_loc)
        a = afs_loc
    end
    if (a < c.amin)
        a = c.amin
    end
    return a
end

"""
    radiation!(fsm, meteo, t)

Snow albedo calculations, surface and canopy net shortwave radiation.

The per-cell physics lives in `radiation_kernel!`, a KernelAbstractions
kernel launched over the whole grid (see `ebalsrf!` for the pattern). The
two former grid loops (albedo, net radiation) are fused into one kernel;
they only communicate through values of the same grid cell, so the fusion is
exact. Calendar tests on `t` are resolved on the host, since `Dates`
operations cannot run inside kernels.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions (read-only)
- `t`: Current simulation time
"""
function radiation!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}, t) where {Tf <: Real, Ti <: Integer}

    (; CANOPY, ALBEDO) = fsm.physics

    # Dates cannot cross into kernels: resolve the calendar test here
    summer_decay = Dates.value(Month(t)) > 4 && Dates.value(Month(t)) < 10

    backend = get_backend(fsm.state.albs)
    kernel! = radiation_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.landuse, fsm.params, meteo,
        CANOPY, ALBEDO, summer_decay;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function radiation_kernel!(
        state, diag, landuse, params::Parameters{Tf}, meteo,
        CANOPY::AbstractCanopy{Tf},
        ALBEDO::AbstractAlbedo{Tf}, summer_decay::Bool,
    ) where {Tf}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    (; dt, tthresh) = params
    (; fsky_terr, fveg, tilefrac, alb0, fsky, scap, trcn, afs, adc) = landuse
    (; albs, Sice, Sliq, fsnow, Sveg, Tsrf) = state
    (; alb, asrf_out, SWveg, SWsrf, SWsci, LWt, LWeff) = diag
    (; LW, Sdif, Sdir, Sdird, Sf, Sf24h, Ta, Tv) = meteo

    if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        # Snow albedo

        SWEtmp = zero(Tf)
        for si in 1:size(Sice, 1)
            SWEtmp += Sice[si, i, j] + Sliq[si, i, j]
        end
        alb_s = (; Tsrf = Tsrf[i, j], fveg = fveg[i, j], trcn = trcn[i, j], fsky = fsky[i, j],
            afs = afs[i, j], adc = adc[i, j], Sdir = Sdir[i, j], Sdif = Sdif[i, j], Sdird = Sdird[i, j],
            Sf = Sf[i, j], Sf24h = Sf24h[i, j], Tv = Tv[i, j], SWEtmp, Tm, dt, summer_decay)
        albs[i, j] = snow_albedo(ALBEDO, albs[i, j], alb_s)

        # Surface and canopy net shortwave radiation

        # Surface albedo
        asrf = albs[i, j] * (Tf(1) - fveg[i, j] * canopy_fsar(CANOPY))
        if (fsnow[i, j] <= eps(Tf))
            asrf = alb0[i, j]
            albs[i, j] = alb0[i, j]
        end

        # Partial snowcover on canopy
        fcans = Tf(0.0)
        if (scap[i, j] > eps(Tf))
            fcans = Sveg[i, j] / scap[i, j]
        end
        aveg = (Tf(1) - fcans) * canopy_avg0(CANOPY) + fcans * canopy_avgs(CANOPY)
        acan = fveg[i, j] * aveg
        # Canopy surface albedo for computing terrain radiation over canopy
        alb[i, j] = fveg[i, j] * aveg + (Tf(1) - fveg[i, j]) * asrf

        # Surface albedo is stored in asurf_out to write in results
        asrf_out[i, j] = alb[i, j]

        # Solar radiation trasmission
        if CANOPY isa NoCanopy
            # No canopy: open-sky transmission
            SWveg[i, j] = Tf(0)
            SWsrf[i, j] = (Tf(1) - alb[i, j]) * (Sdir[i, j] + Sdif[i, j])
            SWsci[i, j] = Sdif[i, j] + Sdir[i, j]
        else
            Sdif_aux = fsky[i, j] / fsky_terr[i, j] * Sdif[i, j]
            tdif = trcn[i, j]
            tdir = Tv[i, j]

            # Effective albedo and net radiation
            alb[i, j] = acan + (Tf(1) - acan) * asrf * tdif^Tf(2)
            if (Sdif_aux + Sdir[i, j] > eps(Tf))
                alb[i, j] = (acan * (Sdif_aux + tdir * Sdir[i, j]) + asrf * tdif * (tdif * Sdif_aux + tdir * Sdir[i, j])) / (Sdif_aux + Sdir[i, j])
            end
            SWsrf[i, j] = (Tf(1) - asrf) * (tdif * Sdif_aux + tdir * Sdir[i, j])
            SWveg[i, j] = ((Tf(1) - tdif) * (Tf(1) - aveg) + tdif * asrf * (Tf(1) - tdif)) * Sdif_aux + (tdir * fveg[i, j] * (Tf(1) - aveg) + tdir * asrf * (Tf(1) - tdif)) * Sdir[i, j]   # local SWR absorption by vegetation correlates with local tdir
            SWsci[i, j] = tdif * Sdif_aux + tdir * Sdir[i, j]
        end

        # Thermal emissions from surroundings
        # Terrain LWR if not calculated later;
        LWt[i, j] = fsky_terr[i, j] * LW[i, j] + (Tf(1) - fsky_terr[i, j]) * sb * Ta[i, j]^Tf(4)

        # LWeff equals LWt except when EBALFOR is used, where terrain impacts are accounted for already
        if CANOPY isa NoCanopy
            LWeff[i, j] = LWt[i, j]
        else
            LWeff[i, j] = LW[i, j]
        end

    end
end
