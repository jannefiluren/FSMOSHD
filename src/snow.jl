# Gsoil = zeros(Nx, Ny)
# Roff = zeros(Nx, Ny)
# meltflux_out = zeros(Nx, Ny)
# Sbsrf = zeros(Nx, Ny)
# Roff_bare = zeros(Nx, Ny)
# Roff_snow = zeros(Nx, Ny)
# fsnow_thres = zeros(Nx, Ny)
# unload = zeros(Nx, Ny)

function snow(fsm::FSM, Rf, Sf, Ta, Ua)

  @unpack HYDROL, DENSTY, OSHDTN, HN_ON, SNFRAC = fsm

  @unpack tthresh = fsm

  @unpack dt = fsm

  @unpack Dzsnow, Dzsoil, Nsmax, Nsoil, Nx, Ny = fsm

  @unpack a_eta, b_eta, c_eta, eta0, eta1, rgr0, rho0, rhob, rhoc, rhof, rcld, rmlt, snda, trho, Wirr = fsm

  @unpack Ds, Nsnow, fsnow, rgrn, Sice, Sliq, Tsnow, Tsoil, Tsrf = fsm
  
  @unpack dem, tilefrac = fsm

  @unpack ksnow, ksoil = fsm

  @unpack Esrf, G, Melt = fsm

  @unpack Gsoil, Roff, meltflux_out, Sbsrf, Roff_bare, Roff_snow, fsnow_thres, unload = fsm

  Gsoil[:, :] .= G[:, :]
  Roff[:, :] .= 0
  meltflux_out[:, :] .= 0
  Roff_bare[:, :] .= Rf[:, :] .* dt .* (1 .- fsnow[:, :])
  Roff_snow[:, :] .= Rf[:, :] .* dt .* fsnow[:, :]

  a = zeros(Nsmax)
  b = zeros(Nsmax)
  c = zeros(Nsmax)
  csnow = zeros(Nsmax)
  dTs = zeros(Nsmax)
  D = zeros(Nsmax)
  E = zeros(Nsmax)
  Gs = zeros(Nsmax)
  rhs = zeros(Nsmax)
  R = zeros(Nsmax)
  S = zeros(Nsmax)
  U = zeros(Nsmax)
  W = zeros(Nsmax)


  # Points with existing snowpack
  for j = 1:Ny
    for i = 1:Nx

      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        Sbsrf[i, j] = 0

        if (fsnow[i, j] > eps(Float64)) # This condition should be equivalent to Nsnow[i,j] > 0

          # Except for point case, apply a minimum threshold of 0.1 to fsnow
          # to avoid 'long tails' in SWE due to slowing down depletion rates
          if (SNFRAC == 3)
            fsnow_thres[i, j] = fsnow[i, j]
          else
            fsnow_thres[i, j] = max(fsnow[i, j], 0.1)
          end

          # Heat conduction
          for k = 1:Nsnow[i, j]
            csnow[k] = (Sice[k, i, j] * hcap_ice + Sliq[k, i, j] * hcap_wat) / fsnow[i, j]
          end
          if (Nsnow[i, j] == 1)
            Gs[1] = 2 / (Ds[1, i, j] / ksnow[1, i, j] + Dzsoil[1] / ksoil[1, i, j])
            dTs[1] = (G[i, j] + Gs[1] * (Tsoil[1, i, j] - Tsnow[1, i, j])) * dt / (csnow[1] + Gs[1] * dt)
          else
            for k = 1:Nsnow[i, j]-1
              Gs[k] = 2 / (Ds[k, i, j] / ksnow[k, i, j] + Ds[k+1, i, j] / ksnow[k+1, i, j])
            end
            a[1] = 0
            b[1] = csnow[1] + Gs[1] * dt
            c[1] = -Gs[1] * dt
            rhs[1] = (G[i, j] - Gs[1] * (Tsnow[1, i, j] - Tsnow[2, i, j])) * dt
            for k = 2:Nsnow[i, j]-1
              a[k] = c[k-1]
              b[k] = csnow[k] + (Gs[k-1] + Gs[k]) * dt
              c[k] = -Gs[k] * dt
              rhs[k] = Gs[k-1] * (Tsnow[k-1, i, j] - Tsnow[k, i, j]) * dt + Gs[k] * (Tsnow[k+1, i, j] - Tsnow[k, i, j]) * dt
            end
            k = Nsnow[i, j]
            Gs[k] = 2 / (Ds[k, i, j] / ksnow[k, i, j] + Dzsoil[1] / ksoil[1, i, j])
            a[k] = c[k-1]
            b[k] = csnow[k] + (Gs[k-1] + Gs[k]) * dt
            c[k] = 0
            rhs[k] = Gs[k-1] * (Tsnow[k-1, i, j] - Tsnow[k, i, j]) * dt + Gs[k] * (Tsoil[1, i, j] - Tsnow[k, i, j]) * dt
            ### call TRIDIAG(Nsnow(i,j),Nsmax,a,b,c,rhs,dTs)
            tridiag!(Nsnow[i, j], Nsmax, a, b, c, rhs, dTs)
          end
          for k = 1:Nsnow[i, j]
            Tsnow[k, i, j] = Tsnow[k, i, j] + dTs[k]
            if (HN_ON)
              Tsnow[k, i, j] = max(Tsnow[k, i, j], (Tm - 40))
            end
          end
          k = Nsnow[i, j]
          Gsoil[i, j] = Gs[k] * (Tsnow[k, i, j] - Tsoil[1, i, j])

          # Convert melting ice to liquid water
          dSice = Melt[i, j] * fsnow_thres[i, j] * dt

          meltflux_out[i, j] = dSice

          for k = 1:Nsnow[i, j]
            coldcont = csnow[k] * (Tm - Tsnow[k, i, j])
            if (coldcont < 0)
              dSice = dSice - fsnow[i, j] * coldcont / Lf
              Tsnow[k, i, j] = Tm
            end
            if (dSice > eps(Float64))
              if (dSice > Sice[k, i, j])  # Layer melts completely
                dSice = dSice - Sice[k, i, j]
                Ds[k, i, j] = 0
                Sliq[k, i, j] = Sliq[k, i, j] + Sice[k, i, j]
                Sice[k, i, j] = 0
              else                       # Layer melts partially
                Ds[k, i, j] = (1 - dSice / Sice[k, i, j]) * Ds[k, i, j]
                Sice[k, i, j] = Sice[k, i, j] - dSice
                Sliq[k, i, j] = Sliq[k, i, j] + dSice
                dSice = 0                # Melt exhausted
              end
            end
          end

          # Remove snow by sublimation 
          dSice = max(Esrf[i, j] * fsnow_thres[i, j], 0.0) * dt
          if (dSice > eps(Float64))
            for k = 1:Nsnow[i, j]
              if (dSice > Sice[k, i, j])  # Layer sublimates completely
                dSice = dSice - Sice[k, i, j]
                Ds[k, i, j] = 0
                Sbsrf[i, j] = Sbsrf[i, j] + Sice[k, i, j]
                Sice[k, i, j] = 0
              else                       # Layer sublimates partially
                Ds[k, i, j] = (1 - dSice / Sice[k, i, j]) * Ds[k, i, j]
                Sice[k, i, j] = Sice[k, i, j] - dSice
                Sbsrf[i, j] = Sbsrf[i, j] + dSice
                dSice = 0                # Sublimation exhausted
              end
            end
          end

          # Snow hydraulics
          # First, unloading snow is added to liquid water if Ta above freezing point
          if (Ta[i, j] >= Tm)
            Roff_bare[i, j] = Roff_bare[i, j] + unload[i, j] * (1 - fsnow[i, j]) # Bare soil fraction
            Roff_snow[i, j] = Roff_snow[i, j] + unload[i, j] * fsnow[i, j] # Snow covered ground fraction
          end
          if (HYDROL == 0)
            # Free-draining snow 
            meltflux_out[i, j] = 0
            for k = 1:Nsnow[i, j]
              Roff_snow[i, j] = Roff_snow[i, j] + Sliq[k, i, j]
              meltflux_out[i, j] = meltflux_out[i, j] + Sliq[k, i, j]
              Sliq[k, i, j] = 0
            end
          elseif (HYDROL == 1)
            # Bucket storage 
            for k = 1:Nsnow[i, j]
              phi = 0
              if (Ds[k, i, j] > eps(Float64))
                phi = 1 - Sice[k, i, j] / (rho_ice * Ds[k, i, j] * fsnow[i, j])
              end
              SliqMax = fsnow[i, j] * rho_wat * Ds[k, i, j] * phi * Wirr
              Sliq[k, i, j] = Sliq[k, i, j] + Roff_snow[i, j]
              Roff_snow[i, j] = 0
              if (Sliq[k, i, j] > SliqMax)       # Liquid capacity exceeded
                Roff_snow[i, j] = Sliq[k, i, j] - SliqMax   # so drainage to next layer
                Sliq[k, i, j] = SliqMax
              end
              coldcont = csnow[k] * (Tm - Tsnow[k, i, j])
              if (coldcont > 0)       # Liquid can freeze
                dSice = min(Sliq[k, i, j], fsnow[i, j] * coldcont / Lf)
                Sliq[k, i, j] = Sliq[k, i, j] - dSice
                Sice[k, i, j] = Sice[k, i, j] + dSice
                meltflux_out[i, j] = meltflux_out[i, j] - dSice
                Tsnow[k, i, j] = Tsnow[k, i, j] + Lf * dSice / csnow[k] / fsnow[i, j]
              end
            end

            if (meltflux_out[i, j] < 0)
              meltflux_out[i, j] = 0
            end

          else  # HYDROL == 2
            # Density-dependent bucket storage
            for k = 1:Nsnow[i, j]

              SliqCap = 0   ###hackhack

              if (Ds[k, i, j] > eps(Float64))
                rhos = Sice[k, i, j] / Ds[k, i, j] / fsnow[i, j]
                SliqCap = 0.03 + 0.07 * (1 - rhos / 200)
                SliqCap = max(SliqCap, 0.03)
              end
              SliqMax = SliqCap * Sice[k, i, j]
              Sliq[k, i, j] = Sliq[k, i, j] + Roff_snow[i, j]
              Roff_snow[i, j] = 0
              if (Sliq[k, i, j] > SliqMax)       # Liquid capacity exceeded
                Roff_snow[i, j] = Sliq[k, i, j] - SliqMax   # so drainage to next layer
                Sliq[k, i, j] = SliqMax
              end
              coldcont = csnow[k] * (Tm - Tsnow[k, i, j])
              if (coldcont > eps(Float64))       # Liquid can freeze
                dSice = min(Sliq[k, i, j], fsnow[i, j] * coldcont / Lf)
                Sliq[k, i, j] = Sliq[k, i, j] - dSice
                Sice[k, i, j] = Sice[k, i, j] + dSice
                # to account for refreezing of melt
                meltflux_out[i, j] = meltflux_out[i, j] - dSice
                Tsnow[k, i, j] = Tsnow[k, i, j] + Lf * dSice / csnow[k] / fsnow[i, j]
              end
            end

            if (meltflux_out[i, j] < 0)
              meltflux_out[i, j] = 0
            end

          end

          # Snow compaction
          if (DENSTY == 0)
            # Fixed snow density
            for k = 1:Nsnow[i, j]
              Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rho0 / fsnow[i, j]
            end
          elseif (DENSTY == 1)
            # Snow compaction with age
            for k = 1:Nsnow[i, j]
              if (Ds[k, i, j] > eps(Float64))
                rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
                if (Tsnow[k, i, j] >= Tm)
                  if (rhos < rmlt)
                    rhos = rmlt + (rhos - rmlt) * exp(-dt / trho)
                  end
                else
                  if (rhos < rcld)
                    rhos = rcld + (rhos - rcld) * exp(-dt / trho)
                  end
                end
                Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rhos / fsnow[i, j]
              end
            end
          elseif (DENSTY == 2)
            # Snow compaction by overburden
            mass = 0
            for k = 1:Nsnow[i, j]
              mass = mass + 0.5 * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
              if (Ds[k, i, j] > eps(Float64))
                rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
                rhos = rhos + (rhos * grav * mass * dt / (eta0 * exp(-(Tsnow[k, i, j] - Tm) / 12.4 + rhos / 55.6)) + dt * rhos * snda * exp((Tsnow[k, i, j] - Tm) / 23.8 - max(rhos - 150, 0.0) / 21.7))
                Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rhos / fsnow[i, j]
              end
              mass = mass + 0.5 * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
            end
          else # DENSTY == 3
            # Snow compaction by overburden, dependent on liquid water content (Crocus B92)
            mass = 0
            for k = 1:Nsnow[i, j]
              mass = mass + 0.5 * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
              if (Ds[k, i, j] > eps(Float64))
                rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
                f1 = 1 / (1 + 600 * Sliq[k, i, j] / (rho_wat * Ds[k, i, j] * fsnow[i, j]))
                f2 = 1.0
                eta = f1 * f2 * eta1 * (rhos / c_eta) * exp(a_eta * (Tm - Tsnow[k, i, j]) + b_eta * rhos)
                rhos = rhos + rhos * grav * mass * dt / eta
                Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rhos / fsnow[i, j]
              end
              mass = mass + 0.5 * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
            end
          end

          # Snow grain growth --> *GM for now, this code feature is not functional because the state variable rgrn is not tracked (no bin output)
          for k = 1:Nsnow[i, j]
            ggr = 2e-13
            if (Tsnow[k, i, j] < Tm)
              if (rgrn[k, i, j] < 1.50e-4)
                ggr = 2e-14
              else
                ggr = 7.3e-8 * exp(-4600 / Tsnow[k, i, j])
              end
            end
            rgrn[k, i, j] = rgrn[k, i, j] + dt * ggr / rgrn[k, i, j]
          end
        end  # Existing snowpack

        # Important check: meltflux_out doesn't track liquid water retention during percolation, while Roff_snow does, so here we use Roff_snow to 
        # constrain meltflux_out; relevant for HYDROL = [1,2]. Minor inconsistencies remain in case of rainfall and percolation through the entire 
        # snowpack in the same timestep. Verified by LQ and GM, Nov 2021
        if (meltflux_out[i, j] > Roff_snow[i, j])
          meltflux_out[i, j] = Roff_snow[i, j]
        end

        # Add bare soil runoff to snowmelt runoff for total runoff
        Roff[i, j] = Roff_snow[i, j] + Roff_bare[i, j]

      end

    end
  end


  # Add snow, recalculate layers
  for j = 1:Ny
    for i = 1:Nx

      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        # Falling snow temperature
        Tsnow0 = min(Ta[i, j], Tm)
        if (HN_ON)
          Tsnow0 = max(Tsnow0, (Tm - 40))
        end


        # Add snowfall and frost to layer 1 with fresh snow density and grain size
        Esnow = 0
        if (Esrf[i, j] < 0 && Tsrf[i, j] < Tm)
          Esnow = fsnow[i, j] * Esrf[i, j]
          Sbsrf[i, j] = Esnow * dt
        end
        dSice = (Sf[i, j] - Esnow) * dt  # Think about how to scale for fsnow...


        # Catch to round infinitesimally small new snow amounts.
        # The small amounts were due to EnKF-assimilated daily precip being downscaled to hourly
        # based on the hourly (COSMO) to daily (COSMO) mass ratio.  Since min daily EnKF was 1mm,
        # when 1mm was downscaled to hourly this mass could become
        # quite small.  When this occurred on bare ground, Tsnow calculation would blow up,
        # place NaN in Tsnow and snow would never again melt through the season
        # Do we still need this Catch in FSM??
        if (Nsnow[i, j] <= 1 && dSice < 0.001 && Sice[1, i, j] < 0.001)
          dSice = round(dSice * 1000 + 0.5) / 1000  ### hack see original code
        end

        if (DENSTY == 0)
          rhonew = rho0
        else
          if (OSHDTN == 0)
            # Initial formulation
            rhonew = max(rhof + rhob * (Ta[i, j] - Tm) + rhoc * Ua[i, j]^0.5, 50.0)
          else # OSHDTN == 1
            # New formulation with decompaction
            rhonew = rhof + rhob * (Ta[i, j] - Tm) + rhoc * Ua[i, j]^0.5
            if (dem[i, j] <= 1000)
              t_decompaction = 24
            elseif (dem[i, j] > 2000)
              t_decompaction = 0
            else
              t_decompaction = 24 + (dem[i, j] - 1000) / (2000 - 1000) * (0 - 24)
            end
            rhonew = 300 + (rhonew - 300) * exp(t_decompaction / 100)
            rhonew = max(rhonew, 50.0)
          end
        end
        if (Sice[1, i, j] + dSice > eps(Float64))
          rgrn[1, i, j] = (Sice[1, i, j] * rgrn[1, i, j] + dSice * rgr0) / (Sice[1, i, j] + dSice)
        end
        Sice[1, i, j] = Sice[1, i, j] + dSice
        snowdepth = sum(Ds[:, i, j]) * fsnow[i, j] + dSice / rhonew

        # Add canopy unloading to layer 1 with bulk snow density and grain size
        rhos = rhof
        if (Ta[i, j] < Tm)  # only if it's cold enough
          mass = sum(Sice[:, i, j]) + sum(Sliq[:, i, j])
          if (snowdepth > eps(Float64))
            rhos = mass / snowdepth
          end
          Sice[1, i, j] = Sice[1, i, j] + unload[i, j]
          snowdepth = snowdepth + unload[i, j] / rhos
        end

        # Store previous snow cover fraction
        fold = fsnow[i, j]
        # Updated Fractional Snow-Covered Area 
        SWEtmp = sum(Sice[:, i, j] + Sliq[:, i, j])
        ### call SNOWCOVERFRACTION(snowdepth,SWEtmp,i,j)



        ##### hachachkacshdfdasfsdafdf

        if SWEtmp > eps(Float64)
          fsnow[i, j] = 1
        else
          fsnow[i, j] = 0
        end




        # Rescale Ds with new snow cover fraction
        if (fsnow[i, j] > eps(Float64))
          # Update surface layer thickness based on new fsnow
          Ds[1, i, j] = Ds[1, i, j] * fold / fsnow[i, j] + dSice / rhonew / fsnow[i, j]
          if (Ta[i, j] < Tm)
            Ds[1, i, j] = Ds[1, i, j] + unload[i, j] / rhos / fsnow[i, j]
          end
        end
        if (Nsnow[i, j] > 1)
          for k = 2:Nsnow[i, j]
            Ds[k, i, j] = Ds[k, i, j] * fold / fsnow[i, j]
          end
        end

        # New snowpack
        if (Nsnow[i, j] == 0 && Sice[1, i, j] > eps(Float64))
          Nsnow[i, j] = 1
          Tsnow[1, i, j] = Tsnow0
        end

        # Store state of old layers
        D[:] .= Ds[:, i, j]
        R[:] .= rgrn[:, i, j]
        S[:] .= Sice[:, i, j]
        W[:] .= Sliq[:, i, j]
        if (fsnow[i, j] > eps(Float64))
          csnow[1] = (Sice[1, i, j] * hcap_ice + Sliq[1, i, j] * hcap_wat) / fsnow[i, j]
          E[1] = csnow[1] * (Tsnow[1, i, j] - Tm) + ((dSice + unload[i, j]) * hcap_ice / fsnow[i, j]) * (Tsnow0 - Tsnow[1, i, j]) # Adjustment given that csnow[1] already includes the new snow
        else
          E[:] .= 0
        end
        if (Nsnow[i, j] > 1)
          for k = 2:Nsnow[i, j]
            csnow[k] = (Sice[k, i, j] * hcap_ice + Sliq[k, i, j] * hcap_wat) / fsnow[i, j]
            E[k] = csnow[k] * (Tsnow[k, i, j] - Tm)
          end
        end
        Nold = Nsnow[i, j]

        # Initialise new layers
        Ds[:, i, j] .= 0
        rgrn[:, i, j] .= 0
        Sice[:, i, j] .= 0
        Sliq[:, i, j] .= 0
        Tsnow[:, i, j] .= Tm
        U[:] .= 0
        Nsnow[i, j] = 0

        if (fsnow[i, j] > eps(Float64))  # Existing or new snowpack

          # Re-assign and count snow layers
          dnew = snowdepth / fsnow[i, j]
          Ds[1, i, j] = dnew
          #k = 1   #hackhackh
          if (Ds[1, i, j] > Dzsnow[1])
            for k = 1:Nsmax
              Ds[k, i, j] = Dzsnow[k]
              dnew = dnew - Dzsnow[k]
              if (dnew <= Dzsnow[k] || k == Nsmax)
                Ds[k, i, j] = Ds[k, i, j] + dnew
                break
              end
            end
          end
          #Nsnow[i, j] = k
          Nsnow[i, j] = sum(Ds[:, i, j] .> 0)   ### hackhack

          # Fill new layers from the top downwards
          knew = 1
          dnew = Ds[1, i, j]
          for kold = 1:Nold
            while true
              if (D[kold] < dnew)
                # All snow from old layer partially fills new layer
                rgrn[knew, i, j] = rgrn[knew, i, j] + S[kold] * R[kold]
                Sice[knew, i, j] = Sice[knew, i, j] + S[kold]
                Sliq[knew, i, j] = Sliq[knew, i, j] + W[kold]
                U[knew] = U[knew] + E[kold]
                dnew = dnew - D[kold]
                break
              else
                # Some snow from old layer fills new layer
                wt = dnew / D[kold]
                rgrn[knew, i, j] = rgrn[knew, i, j] + wt * S[kold] * R[kold]
                Sice[knew, i, j] = Sice[knew, i, j] + wt * S[kold]
                Sliq[knew, i, j] = Sliq[knew, i, j] + wt * W[kold]
                U[knew] = U[knew] + wt * E[kold]
                D[kold] = (1 - wt) * D[kold]
                E[kold] = (1 - wt) * E[kold]
                S[kold] = (1 - wt) * S[kold]
                W[kold] = (1 - wt) * W[kold]
                knew = knew + 1
                if (knew > Nsnow[i, j])
                  break
                end
                dnew = Ds[knew, i, j]
              end
            end
          end

          # Diagnose snow layer temperatures
          for k = 1:Nsnow[i, j]
            csnow[k] = (Sice[k, i, j] * hcap_ice + Sliq[k, i, j] * hcap_wat) / fsnow[i, j]
            Tsnow[k, i, j] = Tm + U[k] / csnow[k]
            if (HN_ON)
              Tsnow[k, i, j] = max(Tsnow[k, i, j], (Tm - 40))
            end
            rgrn[k, i, j] = rgrn[k, i, j] / Sice[k, i, j]
          end
        end # Existing or new snowpack

      end

    end
  end

end
