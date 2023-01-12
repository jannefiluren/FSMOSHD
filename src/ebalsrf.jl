# Esrf = zeros(Nx,Ny)
# Eveg = zeros(Nx,Ny)
# G = zeros(Nx,Ny)
# H = zeros(Nx,Ny)
# Hsrf = zeros(Nx,Ny)
# LE = zeros(Nx,Ny)
# LEsrf = zeros(Nx,Ny)
# LWsci = zeros(Nx,Ny)
# LWveg = zeros(Nx,Ny)
# Melt = zeros(Nx,Ny)
# Rnet = zeros(Nx,Ny)
# Rsrf = zeros(Nx,Ny)

function ebalsrf(fsm::FSM, LW, Ps, Qa, Ta)

  @unpack CANMOD = fsm

  @unpack TILE, tthresh = fsm

  @unpack dt = fsm

  @unpack Nx, Ny = fsm

  @unpack trcn = fsm

  @unpack Sice, Tcan, Tsrf, Tveg = fsm

  @unpack fveg, tilefrac = fsm

  @unpack SWsrf = fsm

  @unpack Ds1, Ts1, ks1 = fsm

  @unpack dTs, Esrf, Eveg, G, H, Hsrf, LE, LEsrf, LWsci, LWveg, Melt, Rnet, Rsrf = fsm

  @unpack KH, KWg = fsm

  for j = 1:Ny
    for i = 1:Nx

      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        if ((CANMOD == 1 && fveg[i, j] == 0) || CANMOD == 0)

          local D #### hack
          #global dTs #### hack

          Tveg[i, j] = Ta[i, j]
          Tcan[i, j] = Ta[i, j]

          # Saturation humidity and density of air
          Qs = qsat(Ps[i, j], Tsrf[i, j])  #call QSAT(Ps[i,j],Tsrf[i,j],Qs)
          Lh = Lv
          if (Tsrf[i, j] < Tm || Sice[1, i, j] > eps(Float64))
            Lh = Ls
          end
          D = Lh * Qs / (Rwat * Tsrf[i, j]^2)
          rho = Ps[i, j] / (Rair * Ta[i, j])

          # Explicit fluxes
          Esrf[i, j] = rho * KWg[i, j] * (Qs - Qa[i, j])
          G[i, j] = 2 * ks1[i, j] * (Tsrf[i, j] - Ts1[i, j]) / Ds1[i, j]
          H[i, j] = cp * rho * KH[i, j] * (Tsrf[i, j] - Ta[i, j])
          LE[i, j] = Lh * Esrf[i, j]
          Melt[i, j] = 0
          Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LW[i, j] - sb * Tsrf[i, j]^4 + (1 - trcn[i, j]) * sb * Tveg[i, j]^4

          # Surface energy balance increments without melt
          dTs[i, j] = (Rnet[i, j] - G[i, j] - H[i, j] - LE[i, j]) / (4 * sb * Tsrf[i, j]^3 + 2 * ks1[i, j] / Ds1[i, j] + rho * (cp * KH[i, j] + Lh * D * KWg[i, j]))
          dE = rho * KWg[i, j] * D * dTs[i, j]
          dG = 2 * ks1[i, j] * dTs[i, j] / Ds1[i, j]
          dH = cp * rho * KH[i, j] * dTs[i, j]
          dR = -4 * sb * Tsrf[i, j]^3 * dTs[i, j]

          # Surface melting
          if (Tsrf[i, j] + dTs[i, j] > Tm && Sice[1, i, j] > eps(Float64))
            Melt[i, j] = sum(Sice[:, i, j]) / dt
            dTs[i, j] = (Rnet[i, j] - G[i, j] - H[i, j] - LE[i, j] - Lf * Melt[i, j]) / (4 * sb * Tsrf[i, j]^3 + 2 * ks1[i, j] / Ds1[i, j] + rho * (cp * KH[i, j] + Ls * D * KWg[i, j]))
            dE = rho * KWg[i, j] * D * dTs[i, j]
            dG = 2 * ks1[i, j] * dTs[i, j] / Ds1[i, j]
            dH = cp * rho * KH[i, j] * dTs[i, j]
            dR = -4 * sb * Tsrf[i, j]^3 * dTs[i, j]
            if (Tsrf[i, j] + dTs[i, j] < Tm)
              Qs = qsat(Ps[i, j], Tm)  #call QSAT(Ps[i,j],Tm,Qs)
              Esrf[i, j] = rho * KWg[i, j] * (Qs - Qa[i, j])
              G[i, j] = 2 * ks1[i, j] * (Tm - Ts1[i, j]) / Ds1[i, j]
              H[i, j] = cp * rho * KH[i, j] * (Tm - Ta[i, j])
              LE[i, j] = Ls * Esrf[i, j]
              Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LW[i, j] - sb * Tm^4 + (1 - trcn[i, j]) * sb * Tveg[i, j]^4
              Melt[i, j] = (Rnet[i, j] - H[i, j] - LE[i, j] - G[i, j]) / Lf
              Melt[i, j] = max(Melt[i, j], 0.0)
              dE = 0.0
              dG = 0.0
              dH = 0.0
              dR = 0.0
              dTs[i, j] = Tm - Tsrf[i, j]
            end
          end

          # In case of glacier without snow, cap Tsrf to 0°C
          # This adjustment:
          #     - assumes the glacier is an infinite heat reservoir.
          #     - does not conserve energy.
          # The excess energy would correspond to glacier melting, which we don't track.
          if (TILE == "glacier")
            if (Tsrf[i, j] + dTs[i, j] > Tm && Sice[1, i, j] <= eps(Float64))
              Qs = qsat(Ps[i, j], Tm)  #call QSAT(Ps[i,j],Tm,Qs)
              Esrf[i, j] = rho * KWg[i, j] * (Qs - Qa[i, j])
              G[i, j] = 2 * ks1[i, j] * (Tm - Ts1[i, j]) / Ds1[i, j]
              H[i, j] = cp * rho * KH[i, j] * (Tm - Ta[i, j])
              LE[i, j] = Ls * Esrf[i, j]
              Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LW[i, j] - sb * Tm^4 + (1 - trcn[i, j]) * sb * Tveg[i, j]^4
              dE = 0.0
              dG = 0.0
              dH = 0.0
              dR = 0.0
              dTs[i, j] = Tm - Tsrf[i, j]
            end
          end

          # Update surface temperature and fluxes
          Esrf[i, j] = Esrf[i, j] + dE
          G[i, j] = G[i, j] + dG
          H[i, j] = H[i, j] + dH
          LE[i, j] = Lh * Esrf[i, j]
          Rnet[i, j] = Rnet[i, j] + dR
          Tsrf[i, j] = Tsrf[i, j] + dTs[i, j]

          # Sublimation limited by amount of snow after melt
          Ssub = sum(Sice[:, i, j]) - Melt[i, j] * dt
          if (Ssub > eps(Float64) && Esrf[i, j] * dt > Ssub)
            Esrf[i, j] = Ssub / dt
            LE[i, j] = Ls * Esrf[i, j]
            H[i, j] = Rnet[i, j] - G[i, j] - LE[i, j] - Lf * Melt[i, j]
          end
          Hsrf[i, j] = H[i, j]
          LEsrf[i, j] = LE[i, j]
          Rsrf[i, j] = Rnet[i, j]

          # Ensure LWsci and LWveg exist as variable even in open runs
          LWsci[i, j] = LW[i, j]
          LWveg[i, j] = 0

          if (CANMOD == 0)
            # Add fluxes from canopy in zero-layer model
            Eveg[i, j] = 0
            if (fveg[i, j] > eps(Float64))
              Eveg[i, j] = -KWv[i, j] * Esrf[i, j] / (KHa[i, j] + KWv[i, j])
              H[i, j] = KHa[i, j] * H[i, j] / (KHa[i, j] + KHv[i, j])
              Lh = Ls
              if (Tveg[i, j] > Tm)
                Lh = Lv
              end
              LE[i, j] = LE[i, j] + Lh * Eveg[i, j]
              Rnet[i, j] = Rnet[i, j] + SWveg[i, j] + (1 - trcn[i, j]) * (LW[i, j] + sb * Tsrf[i, j]^4 - 2 * sb * Tveg[i, j]^4)
            end
          end
        end

      end

    end
  end

end