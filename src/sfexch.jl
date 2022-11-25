# KH = zeros(Nx, Ny)
# KHa = zeros(Nx, Ny)
# KHg = zeros(Nx, Ny)
# KHv = zeros(Nx, Ny)
# KWg = zeros(Nx, Ny)
# KWv = zeros(Nx, Ny)
# Usc = zeros(Nx, Ny)

function sfexch(fsm::FSM)

  @unpack CANMOD, ZOFFST, EXCHNG, OSHDTN, SNFRAC = fsm

  @unpack TILE, tthresh = fsm

  @unpack zT, zU = fsm

  @unpack Nx, Ny = fsm

  @unpack bstb, cden, cveg, gsnf, rchd, rchz, z0sn, wcan, zsub, zgf, zgr, khcf = fsm

  @unpack VAI, z0sf = fsm

  @unpack Qcan, fsnow, Sice, Sveg, Tcan, Tsrf, Tveg, Ds = fsm

  @unpack dem, fveg, fves, hcan, tilefrac = fsm

  @unpack KH, KHa, KHg, KHv, KWg, KWv, Usc = fsm

  @unpack gs1 = fsm
  
  for j = 1:Ny
    for i = 1:Nx
      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        if (ZOFFST == 0)
          # Heights specified above ground
          zU1 = zU
          zT1 = zT
        else # ZOFFST == 1
          # Heights specified above canopy top
          zU1 = zU + hcan[i, j]
          zT1 = zT + hcan[i, j]
        end

        # Ground roughness length
        if (TILE == "glacier")
          z0loc = 0.0009
        else
          if (OSHDTN == 0 || TILE == "forest")
            z0loc = z0sn
          else # OSHDTN == 1
            if (dem[i, j] >= 2300)
              z0loc = 0.003
            elseif (dem[i, j] >= 1500)
              z0loc = 0.03 + (dem[i, j] - 1500) / (2300 - 1500) * (0.003 - 0.03)
            elseif (dem[i, j] >= 1200)  # simple linear b/w two above values
              z0loc = 0.2 + (dem[i, j] - 1200) / (1500 - 1200) * (0.03 - 0.2)
            else
              z0loc = 0.2
            end
          end
        end
        z0g = z0loc

        # BC, stabilize the tuning point runs by using a Ds threshold instead of fsnow.
        # TODO: test the impact for the grid points.
        if (SNFRAC == 3)
          if (sum(Ds[:, i, j]) <= 0.05)
            z0g = z0sf[i, j]
          end
        else
          if (fsnow[i, j] <= eps(Float64))
            z0g = z0sf[i, j]
          end
        end

        # Additional roughness lengths and friction velocity
        if (EXCHNG == 2) # Forest - specific adjustment *GM
          # Open
          z0h = 0.1 * z0g
          ustar = vkman * Ua[i, j] / log(zU / z0g)
          rgo = log(zT / z0h) / (vkman * ustar)

          # Forest
          if (fveg[i, j] > eps(Float64))
            z0g = (zgf + zgr * fveg[i, j]) * z0g
            z0h = 0.1 * z0g
            dh = rchd * hcan[i, j]
            z0v = rchz * hcan[i, j]
            ustar = vkman * Ua[i, j] / log((zU1 - dh) / z0v)
            Uh = (ustar / vkman) * log((hcan[i, j] - dh) / z0v)
            KHh = vkman * ustar * (hcan[i, j] - dh)
            Usf = exp(wcan * (zsub / hcan[i, j] - 1)) * Uh
          end
        else
          z0v = rchz * hcan[i, j]
          z0 = (z0v^fveg[i, j]) * (z0g^(1 - fveg[i, j]))
          z0h = 0.1 * z0
          dh = fveg[i, j] * rchd * hcan[i, j]
          CD = (vkman / log((zU1 - dh) / z0))^2
          ustar = sqrt(CD) * Ua[i, j]
        end
        Uso = Ua[i, j] * log(zsub / z0g) / log(zU / z0g)

        if (EXCHNG == 0)
          # No stability adjustment
          fh = 1
          Ric = 0
        end
        if (EXCHNG == 1)
          # Stability adjustment (Louis et al. 1982, quoted by Beljaars 1992)
          Tint = fveg[i, j] * Tveg[i, j] + (1 - fveg[i, j]) * Tsrf[i, j]
          RiB = grav * (Ta[i, j] - Tint) * (zU1 - dh)^2 / ((zT1 - dh) * Ta[i, j] * Ua[i, j]^2)
          if (RiB > 0.2)
            RiB = 0.2 # New maximum threshold for RiB
          end
          if (RiB > 0)
            fh = 1 / (1 + 3 * bstb * RiB * sqrt(1 + bstb * RiB))
          else
            fh = 1 - 3 * bstb * RiB / (1 + 3 * bstb^2 * CD * sqrt(-RiB * zU1 / z0))
          end
          Ric = grav * (Tcan[i, j] - Tsrf[i, j]) * hcan[i, j] / (Tcan[i, j] * ustar^2)
          Ric = max(min(Ric, 10.0), 0.0)
        end
        # Note that currently, fh and Ric are not used in EXCHNG == 2, i.e. no stability correction

        # Eddy diffusivities
        if (fveg[i, j] == 0)
          KH[i, j] = fh * vkman * ustar / log(zT1 / z0h)
          Qs = qsat(Ps[i, j], Tsrf[i, j])  #call QSAT(Ps[i,j],Tsrf[i,j],Qs)
          if (Sice[1, i, j] > eps(Float64) || Qa[i, j] > Qs)
            KWg[i, j] = KH[i, j]
          else
            KWg[i, j] = gs1[i, j] * KH[i, j] / (gs1[i, j] + KH[i, j])
          end
          Usc[i, j] = Uso
        else
          if (EXCHNG == 2)
            rad = (log((zT1 - dh) / (hcan[i, j] - dh)) / (vkman * ustar) + hcan[i, j] * (exp(wcan * (1 - (z0v + dh) / hcan[i, j])) - 1) / (wcan * KHh)) / khcf
            KHa[i, j] = sqrt(fves[i, j]) / rad
            Usub = sqrt(fves[i, j]) * Usf + (1 - sqrt(fves[i, j])) * Uso
            Usub = max(Usub, 0.1)
            rgd = 1 / (vkman^2 * Usub) * log(zsub / z0h) * log(zsub / z0g)
            KHg[i, j] = 1 / rgd
            Uc = exp(wcan * ((z0v + dh) / hcan[i, j] - 1)) * Uh
            KHv[i, j] = VAI[i, j] * sqrt(Uc) / cveg
            Usc[i, j] = Usub
          else
            KHa[i, j] = fh * vkman * ustar / log((zT1 - dh) / z0)
            KHg[i, j] = vkman * ustar * ((1 - fveg[i, j]) * fh / log(z0 / z0h) + fveg[i, j] * cden / (1 + 0.5 * Ric))
            KHv[i, j] = sqrt(ustar) * VAI[i, j] / cveg
          end
          Qs = qsat(Ps[i, j], Tsrf[i, j])  #call QSAT(Ps[i,j],Tsrf[i,j],Qs)
          if (Qcan[i, j] > Qs)
            KWg[i, j] = KHg[i, j]
          else
            KWg[i, j] = gs1[i, j] * KHg[i, j] / (gs1[i, j] + KHg[i, j])
          end
          Qs = qsat(Ps[i, j], Tveg[i, j])  #call QSAT(Ps[i,j],Tveg[i,j],Qs)
          if (Sveg[i, j] > eps(Float64) || Qcan[i, j] > Qs)
            KWv[i, j] = KHv[i, j]
          else
            KWv[i, j] = gsnf * KHv[i, j] / (gsnf + KHv[i, j])
          end

          if (CANMOD == 0)
            # Combined resistances for 0-layer canopy model
            KH[i, j] = KHg[i, j] * (KHa[i, j] + KHv[i, j]) / (KHa[i, j] + KHg[i, j] + KHv[i, j])
            KWg[i, j] = KWg[i, j] * (KHa[i, j] + KWv[i, j]) / (KHa[i, j] + KWg[i, j] + KWv[i, j])
          end
        end

      end
    end
  end

end