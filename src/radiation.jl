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

    @unpack Nx, Ny, dt = fsm

    @unpack tthresh, fsky_terr, fveg, tilefrac = fsm

    @unpack asmn, avg0, avgs, Talb, tcld, tmlt, adfs, adfl, fsar, Sfmin = fsm

    @unpack alb0, fsky, scap, trcn = fsm

    @unpack albs, Sice, Sliq, fsnow, Sveg, Tsrf = fsm

    @unpack ALBEDO, CANMOD, ALRADT = fsm

    @unpack alb, asrf_out, SWveg, SWsrf, SWsci, LWt, LWeff = fsm

    @unpack adm, adc, afs = fsm

    @unpack LW, Sdif, Sdir, Sdird, Sf, Sf24h, Ta, Tv = meteo

    # Dates cannot cross into kernels: resolve the calendar test here
    # (forest adjustment of the prognostic albedo decay time, ALBEDO == 1)
    summer_decay = Dates.value(Month(t)) > 4 && Dates.value(Month(t)) < 10

    backend = get_backend(albs)
    kernel! = radiation_kernel!(backend)
    kernel!(
        albs, alb, asrf_out, SWveg, SWsrf, SWsci, LWt, LWeff,
        fsky_terr, fveg, tilefrac, alb0, fsky, scap, trcn,
        Sice, Sliq, fsnow, Sveg, Tsrf, adc, afs,
        LW, Sdif, Sdir, Sdird, Sf, Sf24h, Ta, Tv,
        dt, tthresh, asmn, avg0, avgs, Talb, tcld, tmlt, adfs, adfl, fsar, Sfmin, adm,
        ALBEDO, CANMOD, ALRADT, summer_decay;
        ndrange = (Int(Nx), Int(Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function radiation_kernel!(
        albs, alb, asrf_out, SWveg, SWsrf, SWsci, LWt, LWeff,
        fsky_terr, fveg, tilefrac, alb0,
        fsky, scap, trcn,
        Sice, Sliq, fsnow, Sveg, Tsrf,
        adc, afs,
        LW, Sdif, Sdir, Sdird, Sf,
        Sf24h, Ta, Tv,
        dt::Tf, tthresh::Tf, asmn::Tf, avg0::Tf, avgs::Tf, Talb::Tf, tcld::Tf,
        tmlt::Tf, adfs::Tf, adfl::Tf, fsar::Tf, Sfmin::Tf, adm::Tf,
        ALBEDO::Ti, CANMOD::Ti, ALRADT::Ti, summer_decay::Bool,
    ) where {Tf, Ti}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        # Snow albedo

        adc_loc = adc[i, j]
        adm_loc = adm
        afs_loc = afs[i, j]

        if (ALBEDO == 0)
            # Diagnostic
            albs[i, j] = asmn + (afs_loc - asmn) * (Tsrf[i, j] - Tm) / Talb
            if (albs[i, j] < min(afs_loc, asmn))
                albs[i, j] = min(afs_loc, asmn)
            end
            if (albs[i, j] > max(afs_loc, asmn))
                albs[i, j] = max(afs_loc, asmn)
            end

        elseif (ALBEDO == 1)
            # Prognostic
            tau = tcld
            if (Tsrf[i, j] >= Tm)
                tau = tmlt
            end
            # Forest adjustments -> not yet properly tested for OSHD but option currently unused
            if summer_decay
                tau = Tf(70.0) * Tf(3600.0)
            end

            if fveg[i, j] > Tf(0) && Sdir[i, j] > eps(Tf)
                tau = tau / ((Tf(1) - trcn[i, j] * fsky[i, j]) * (Tf(1) + adfl * Tv[i, j]) + adfs * Tv[i, j])
            elseif fveg[i, j] > Tf(0) && Sdif[i, j] > eps(Tf)
                tau = tau / ((Tf(1) - trcn[i, j] * fsky[i, j]) + adfs * trcn[i, j] * fsky[i, j])
            elseif (fveg[i, j] > Tf(0) && (Sdir[i, j] + Sdif[i, j] <= eps(Tf)))
                tau = tau / (Tf(2.0) - trcn[i, j] * fsky[i, j])
            end

            rt = Tf(1) / tau + Sf[i, j] / Sfmin
            alim = (asmn / tau + Sf[i, j] * afs_loc / Sfmin) / rt
            albs[i, j] = alim + (albs[i, j] - alim) * exp(-rt * dt)
            if (albs[i, j] < min(afs_loc, asmn))
                albs[i, j] = min(afs_loc, asmn)
            end
            if (albs[i, j] > max(afs_loc, asmn))
                albs[i, j] = max(afs_loc, asmn)
            end

        else # ALBEDO == 2
            # Prognostic, tuned, copied from JIM
            SWEtmp = Tf(0.0)
            for si in 1:size(Sice, 1)
                SWEtmp += Sice[si, i, j] + Sliq[si, i, j]
            end

            # BC 08.23: aspect-dependent albedo tuning. Activated for oper season 2024
            # or optionally.
            if (ALRADT == 1)
                # BC Oct 23: Jan's suggestion: modify only when the decay rate should be increased (ad* DECREASE), not decreased
                if ((Sdir[i, j] > eps(Tf)) && (Sdird[i, j] < Sdir[i, j]))
                    adm_loc = adm_loc * (Sdird[i, j]) / (Sdir[i, j])
                    adc_loc = adc_loc * (Sdird[i, j]) / (Sdir[i, j])
                    if (adm_loc < eps(Tf))
                        adm_loc = eps(Tf)
                    end
                    if (adc_loc < eps(Tf))
                        adc_loc = eps(Tf)
                    end
                end
            end

            if (Tsrf[i, j] >= Tm)
                albs[i, j] = (albs[i, j] - asmn) * exp(-(dt / Tf(3600)) / adm_loc) + asmn
            else
                albs[i, j] = albs[i, j] - (dt / Tf(3600)) / adc_loc
            end
            if (SWEtmp < Tf(75.0)) # more stuff showing on and up through snow
                afs_loc *= Tf(0.8)
            end
            # Reset to fresh snow albedo (wasn't originally available; only else term)
            if ((Sf[i, j] * dt) > Tf(0.0) && Sf24h[i, j] > Sfmin)
                albs[i, j] = afs_loc
            else
                albs[i, j] = albs[i, j] + (afs_loc - albs[i, j]) * Sf[i, j] * dt / Sfmin
            end
            ## End Adjustments
            if (albs[i, j] > afs_loc)
                albs[i, j] = afs_loc
            end
            if (albs[i, j] < asmn)
                albs[i, j] = asmn
            end

        end

        # Surface and canopy net shortwave radiation

        # Surface albedo
        asrf = albs[i, j] * (Tf(1) - fveg[i, j] * fsar)
        if (fsnow[i, j] <= eps(Tf))
            asrf = alb0[i, j]
            albs[i, j] = alb0[i, j]
        end

        # Partial snowcover on canopy
        fcans = Tf(0.0)
        if (scap[i, j] > eps(Tf))
            fcans = Sveg[i, j] / scap[i, j]
        end
        aveg = (Tf(1) - fcans) * avg0 + fcans * avgs
        acan = fveg[i, j] * aveg
        # Canopy surface albedo for computing terrain radiation over canopy
        alb[i, j] = fveg[i, j] * aveg + (Tf(1) - fveg[i, j]) * asrf

        # Surface albedo is stored in asurf_out to write in results
        asrf_out[i, j] = alb[i, j]

        # Solar radiation trasmission
        if (CANMOD == 0)
            SWveg[i, j] = Tf(0)
            SWsrf[i, j] = (Tf(1) - alb[i, j]) * (Sdir[i, j] + Sdif[i, j])
            SWsci[i, j] = Sdif[i, j] + Sdir[i, j]
        end

        if (CANMOD == 1)
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
        if (CANMOD == 0 || fveg[i, j] == 0)
            LWeff[i, j] = LWt[i, j]
        else
            LWeff[i, j] = LW[i, j]
        end

    end
end
