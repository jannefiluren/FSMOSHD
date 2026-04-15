"""
    ebalsrf!(fsm, meteo)

Surface energy balance solution for open and non-forest tiles.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions
"""
function ebalsrf!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}) where {Tf<:Real, Ti<:Integer}

  @unpack_constants(Tf)

  @unpack CANMOD, SNTRAN, SNSLID = fsm

  @unpack TILE, tthresh = fsm

  @unpack dt = fsm

  @unpack Nx, Ny = fsm

  @unpack trcn = fsm

  @unpack Sice, Tcan, Tsrf, Tveg = fsm

  @unpack fveg, tilefrac, glacierfrac = fsm

  @unpack SWsrf = fsm

  @unpack Ds1, Ts1, ks1 = fsm

  @unpack Esrf, Eveg, G, H, Hsrf, LE, LEsrf, LWsci, LWveg, Melt, Rnet, Rsrf = fsm

  @unpack KH, KWg, KHa, KHv, KWv, SWveg = fsm

  @unpack LW, Ps, Qa, Ta = meteo

  for j = 1:Ny
    for i = 1:Nx

      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        if ((CANMOD == 1 && fveg[i, j] == 0) || CANMOD == 0)

          Tveg[i, j] = Ta[i, j]
          Tcan[i, j] = Ta[i, j]

          # Saturation humidity and density of air
          Qs = qsat(Ps[i, j], Tsrf[i, j])
          Lh = Lv
          if (Tsrf[i, j] < Tm || Sice[1, i, j] > eps(Tf))
            Lh = Ls
          end
          D = Lh * Qs / (Rwat * Tsrf[i, j]^Tf(2))
          rho = Ps[i, j] / (Rair * Ta[i, j])

          # Explicit fluxes
          Esrf[i, j] = rho * KWg[i, j] * (Qs - Qa[i, j])
          G[i, j] = Tf(2) * ks1[i, j] * (Tsrf[i, j] - Ts1[i, j]) / Ds1[i, j]
          H[i, j] = cp * rho * KH[i, j] * (Tsrf[i, j] - Ta[i, j])
          LE[i, j] = Lh * Esrf[i, j]
          Melt[i, j] = Tf(0)
          Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LW[i, j] - sb * Tsrf[i, j]^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)

          # Surface energy balance increments without melt
          dTs = (Rnet[i, j] - G[i, j] - H[i, j] - LE[i, j]) / (Tf(4) * sb * Tsrf[i, j]^Tf(3) + Tf(2) * ks1[i, j] / Ds1[i, j] + rho * (cp * KH[i, j] + Lh * D * KWg[i, j]))
          dE = rho * KWg[i, j] * D * dTs
          dG = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
          dH = cp * rho * KH[i, j] * dTs
          dR = Tf(-4) * sb * Tsrf[i, j]^Tf(3) * dTs

          # Surface melting
          if (Tsrf[i, j] + dTs > Tm && Sice[1, i, j] > eps(Tf))
            Melt[i, j] = sum(@view Sice[:, i, j]) / dt
            dTs = (Rnet[i, j] - G[i, j] - H[i, j] - LE[i, j] - Lf * Melt[i, j]) / (Tf(4) * sb * Tsrf[i, j]^Tf(3) + Tf(2) * ks1[i, j] / Ds1[i, j] + rho * (cp * KH[i, j] + Ls * D * KWg[i, j]))
            dE = rho * KWg[i, j] * D * dTs
            dG = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
            dH = cp * rho * KH[i, j] * dTs
            dR = Tf(-4) * sb * Tsrf[i, j]^Tf(3) * dTs
            if (Tsrf[i, j] + dTs < Tm)
              Qs = qsat(Ps[i, j], Tm)
              Esrf[i, j] = rho * KWg[i, j] * (Qs - Qa[i, j])
              G[i, j] = Tf(2) * ks1[i, j] * (Tm - Ts1[i, j]) / Ds1[i, j]
              H[i, j] = cp * rho * KH[i, j] * (Tm - Ta[i, j])
              LE[i, j] = Ls * Esrf[i, j]
              Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LW[i, j] - sb * Tm^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
              Melt[i, j] = (Rnet[i, j] - H[i, j] - LE[i, j] - G[i, j]) / Lf
              Melt[i, j] = max(Melt[i, j], Tf(0.0))
              dE = Tf(0.0)
              dG = Tf(0.0)
              dH = Tf(0.0)
              dR = Tf(0.0)
              dTs = Tm - Tsrf[i, j]
            end
          end

          # In case of glacier without snow, cap Tsrf to 0°C
          # This adjustment:
          #     - assumes the glacier is an infinite heat reservoir.
          #     - does not conserve energy.
          # The excess energy would correspond to glacier melting, which we don't track.
          if (TILE == "glacier" || ((SNTRAN == 1 || SNSLID == 1) && glacierfrac[i,j] > eps(Tf)))
            if (Tsrf[i, j] + dTs > Tm && Sice[1, i, j] <= eps(Tf))
              Qs = qsat(Ps[i, j], Tm)
              Esrf[i, j] = rho * KWg[i, j] * (Qs - Qa[i, j])
              G[i, j] = Tf(2) * ks1[i, j] * (Tm - Ts1[i, j]) / Ds1[i, j]
              H[i, j] = cp * rho * KH[i, j] * (Tm - Ta[i, j])
              LE[i, j] = Ls * Esrf[i, j]
              Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LW[i, j] - sb * Tm^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
              dE = Tf(0.0)
              dG = Tf(0.0)
              dH = Tf(0.0)
              dR = Tf(0.0)
              dTs = Tm - Tsrf[i, j]
            end
          end

          # Update surface temperature and fluxes
          Esrf[i, j] = Esrf[i, j] + dE
          G[i, j] = G[i, j] + dG
          H[i, j] = H[i, j] + dH
          LE[i, j] = Lh * Esrf[i, j]
          Rnet[i, j] = Rnet[i, j] + dR
          Tsrf[i, j] = Tsrf[i, j] + dTs

          # Sublimation limited by amount of snow after melt
          Ssub = sum(@view Sice[:, i, j])
          Ssub -= Melt[i, j] * dt
          if (Ssub > eps(Tf) && Esrf[i, j] * dt > Ssub)
            Esrf[i, j] = Ssub / dt
            LE[i, j] = Ls * Esrf[i, j]
            H[i, j] = Rnet[i, j] - G[i, j] - LE[i, j] - Lf * Melt[i, j]
          end
          Hsrf[i, j] = H[i, j]
          LEsrf[i, j] = LE[i, j]
          Rsrf[i, j] = Rnet[i, j]

          # Ensure LWsci and LWveg exist as variable even in open runs
          LWsci[i, j] = LW[i, j]
          LWveg[i, j] = Tf(0.0)

          if (CANMOD == 0)
            # Add fluxes from canopy in zero-layer model
            Eveg[i, j] = Tf(0.0)
            if (fveg[i, j] > eps(Tf))
              Eveg[i, j] = -KWv[i, j] * Esrf[i, j] / (KHa[i, j] + KWv[i, j])
              H[i, j] = KHa[i, j] * H[i, j] / (KHa[i, j] + KHv[i, j])
              Lh = Ls
              if (Tveg[i, j] > Tm)
                Lh = Lv
              end
              LE[i, j] = LE[i, j] + Lh * Eveg[i, j]
              Rnet[i, j] = Rnet[i, j] + SWveg[i, j] + (Tf(1) - trcn[i, j]) * (LW[i, j] + sb * Tsrf[i, j]^Tf(4) - Tf(2) * sb * Tveg[i, j]^Tf(4))
            end
          end
        end

      end

    end
  end

end