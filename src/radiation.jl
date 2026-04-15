"""
    radiation!(fsm, meteo, t)

Snow albedo calculations, surface and canopy net shortwave radiation.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions  
- `t`: Current simulation time
"""
function radiation!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}, t) where {Tf<:Real, Ti<:Integer}

  @unpack_constants(Tf)

  @unpack Nx, Ny, dt = fsm
  
  @unpack tthresh, fsky_terr, fveg, dem, tilefrac = fsm

  @unpack asmx, asmn, avg0, avgs, Talb, tcld, tmlt, adfs, adfl, fsar, Sfmin = fsm

  @unpack alb0,fsky,scap, trcn = fsm

  @unpack albs, Sice, Sliq, fsnow,Sveg, Tsnow, Tsrf = fsm

  @unpack ALBEDO, CANMOD, ALRADT, ALPERT = fsm

  @unpack alb, asrf_out, SWveg, SWsrf, SWsci, LWt = fsm

  @unpack adm, adc, afs = fsm

  @unpack LW, Sdif, Sdir, Sdird, Sf, Sf24h, Ta, Tv = meteo
  
  # Snow albedo
  for j = 1:Ny
    for i = 1:Nx

      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        adc_loc = adc[i,j]
        adm_loc = adm
        afs_loc = afs[i,j]

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
          if (Dates.value(Month(t)) > 4 && Dates.value(Month(t)) < 10)
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
            if ((Sdir[i,j] > eps(Tf)) && (Sdird[i,j] < Sdir[i,j]))
              adm_loc = adm_loc * (Sdird[i,j])/(Sdir[i,j])
              adc_loc = adc_loc * (Sdird[i,j])/(Sdir[i,j])
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
            afs_loc *= Tf(0.80)
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
      end
    end
  end

  # Surface and canopy net shortwave radiation
  for j = 1:Ny
    for i = 1:Nx
      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        # Surface albedo
        asrf = albs[i, j] * (Tf(1) - fveg[i, j] * fsar)
        if (fsnow[i, j] <= eps(Tf))
          asrf = alb0[i, j]
          albs[i,j] = alb0[i,j]
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

        # LW overwritten by LWt only if EBALFOR is used, where terrain impacts are accounted for already 
        if (CANMOD == 0 || fveg[i, j] == 0)
          LW[i, j] = LWt[i, j]
        end

      end

    end

  end

end