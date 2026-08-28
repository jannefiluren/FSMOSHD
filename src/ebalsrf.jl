"""
    ebalsrf!(fsm, meteo)

Surface energy balance solution for open and non-forest tiles.

The per-cell physics lives in `ebalsrf_kernel!`, a KernelAbstractions kernel
launched over the whole grid: it runs partitioned across Julia threads on the
CPU and natively on GPU backends, with the same code and (on the CPU)
bit-identical results to the former plain loops.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions (read-only)
"""
function ebalsrf!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    @unpack TILE, tthresh = fsm

    @unpack dt = fsm

    @unpack Nx, Ny = fsm

    @unpack trcn = fsm

    @unpack Sice, Tcan, Tsrf, Tveg = fsm

    @unpack fveg, tilefrac = fsm

    @unpack SWsrf = fsm

    @unpack Ds1, Ts1, ks1 = fsm

    @unpack Esrf, Eveg, G, H, Hsrf, LE, LEsrf, LWsci, LWveg, Melt, Rnet, Rsrf = fsm

    @unpack KH, KWg, KHa, KHv, KWv, SWveg = fsm

    @unpack Qa, LWeff = fsm

    @unpack Ps, Ta = meteo

    # Strings cannot cross into kernels: resolve the tile test here
    glacier_tile = TILE == "glacier"

    backend = get_backend(Tsrf)
    kernel! = ebalsrf_kernel!(backend)
    kernel!(
        Tveg, Tcan, Tsrf, Esrf, Eveg, G, H, Hsrf, LE, LEsrf, LWsci, LWveg, Melt, Rnet, Rsrf,
        Sice, trcn, fveg, tilefrac, SWsrf, SWveg, Ds1, Ts1, ks1,
        KH, KWg, KHa, KHv, KWv, Qa, LWeff, Ps, Ta,
        dt, tthresh, glacier_tile;
        ndrange = (Int(Nx), Int(Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function ebalsrf_kernel!(
        Tveg, Tcan, Tsrf, Esrf, Eveg, G, H, Hsrf, LE, LEsrf, LWsci, LWveg, Melt, Rnet, Rsrf,
        Sice, trcn, fveg, tilefrac,
        SWsrf, SWveg, Ds1, Ts1, ks1,
        KH, KWg, KHa, KHv, KWv,
        Qa, LWeff, Ps, Ta,
        dt::Tf, tthresh::Tf, glacier_tile::Bool,
    ) where {Tf, Ti}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        if (fveg[i, j] == 0)

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
            Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LWeff[i, j] - sb * Tsrf[i, j]^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)

            # Surface energy balance increments without melt
            dTs = (Rnet[i, j] - G[i, j] - H[i, j] - LE[i, j]) / (Tf(4) * sb * Tsrf[i, j]^Tf(3) + Tf(2) * ks1[i, j] / Ds1[i, j] + rho * (cp * KH[i, j] + Lh * D * KWg[i, j]))
            dE = rho * KWg[i, j] * D * dTs
            dG = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
            dH = cp * rho * KH[i, j] * dTs
            dR = Tf(-4) * sb * Tsrf[i, j]^Tf(3) * dTs

            # Surface melting
            if (Tsrf[i, j] + dTs > Tm && Sice[1, i, j] > eps(Tf))
                Melt[i, j] = column_sum(Sice, i, j) / dt
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
                    Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LWeff[i, j] - sb * Tm^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
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
            if glacier_tile
                if (Tsrf[i, j] + dTs > Tm && Sice[1, i, j] <= eps(Tf))
                    Qs = qsat(Ps[i, j], Tm)
                    Esrf[i, j] = rho * KWg[i, j] * (Qs - Qa[i, j])
                    G[i, j] = Tf(2) * ks1[i, j] * (Tm - Ts1[i, j]) / Ds1[i, j]
                    H[i, j] = cp * rho * KH[i, j] * (Tm - Ta[i, j])
                    LE[i, j] = Ls * Esrf[i, j]
                    Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LWeff[i, j] - sb * Tm^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
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
            Ssub = column_sum(Sice, i, j)
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
            LWsci[i, j] = LWeff[i, j]
            LWveg[i, j] = Tf(0.0)

        end

    end
end
