"""
    ebalfor!(fsm, meteo)

Energy balance solution for forest tiles with 1-layer canopy model.

The per-cell physics lives in `ebalfor_kernel!`, a KernelAbstractions kernel
launched over the whole grid (see `ebalsrf!` for the pattern). Each cell
solves its own 4x4 linear system (`ludcmp!`) using kernel-local
`MMatrix`/`MVector` scratch, so the routine is thread-safe per cell.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions (read-only)
"""
function ebalfor!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    backend = get_backend(fsm.state.Tsrf)
    kernel! = ebalfor_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.landuse, fsm.params, meteo;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing

end

# inbounds = true keeps the kernel-local MMatrix/MVector scratch off the
# heap; see the note on soil_kernel! (a raw @inbounds block inside a @kernel
# body must not be used - it corrupts the KernelAbstractions CPU
# transformation)
@kernel inbounds = true function ebalfor_kernel!(
        state, diag, landuse, params::Parameters{Tf}, meteo,
    ) where {Tf}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    (; dt, tthresh) = params
    (; canh, fsky, trcn, tilefrac) = landuse
    (; Qcan, Sice, Tcan, Tsrf, Tveg) = state
    (; Esrf, Eveg, G, H, Hsrf, LE, LEsrf, LWsci, LWveg, Melt, Rnet, Rsrf,
       Ds1, KHa, KHg, KHv, KWg, KWv, ks1, SWsrf, SWveg, Ts1, Tveg0, Qa, LWeff) = diag
    (; Ps, Ta) = meteo

    # 1-layer canopy model
    if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest


        # Kernel-local scratch for the 4x4 energy balance system
        A_ebal = zero(MMatrix{4, 4, Tf})
        Acp_ebal = zero(MMatrix{4, 4, Tf})
        b_ebal = zero(MVector{4, Tf})
        x_ebal = zero(MVector{4, Tf})
        vv_ebal = zero(MVector{4, Tf})
        indx_ebal = zero(MVector{4, Int32})

        # Saturation humidity and density of air
        Qsrf = qsat(Ps[i, j], Tsrf[i, j])
        Lsrf = Ls
        if (Tsrf[i, j] > Tm)
            Lsrf = Lv
        end
        Dsrf = Lsrf * Qsrf / (Rwat * Tsrf[i, j]^Tf(2))
        Qveg = qsat(Ps[i, j], Tveg[i, j])
        Lveg = Ls
        if (Tveg[i, j] > Tm)
            Lveg = Lv
        end
        Dveg = Lveg * Qveg / (Rwat * Tveg[i, j]^Tf(2))
        rho = Ps[i, j] / (Rair * Ta[i, j])

        # Explicit fluxes
        E = rho * KHa[i, j] * (Qcan[i, j] - Qa[i, j])
        Esrf[i, j] = rho * KWg[i, j] * (Qsrf - Qcan[i, j])
        Eveg[i, j] = rho * KWv[i, j] * (Qveg - Qcan[i, j])
        G[i, j] = Tf(2) * ks1[i, j] * (Tsrf[i, j] - Ts1[i, j]) / Ds1[i, j]
        H[i, j] = rho * cp * KHa[i, j] * (Tcan[i, j] - Ta[i, j])
        Hsrf[i, j] = rho * cp * KHg[i, j] * (Tsrf[i, j] - Tcan[i, j])
        Hveg = rho * cp * KHv[i, j] * (Tveg[i, j] - Tcan[i, j])
        LE[i, j] = Lsrf * Esrf[i, j] + Lveg * Eveg[i, j]
        Melt[i, j] = Tf(0)
        Rsrf[i, j] = SWsrf[i, j] + trcn[i, j] * (fsky[i, j] * LWeff[i, j] + (Tf(1) - fsky[i, j]) * sb * Ta[i, j]^Tf(4)) - sb * Tsrf[i, j]^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)        # with near and distant canopy contributions
        # Rsrf[i,j] = SWsrf[i,j] + trcn[i,j]*LW[i,j] - sb*Tsrf[i,j]^4 + (1 - trcn[i,j])*sb*Tveg[i,j]^4  ! original formulations
        Rveg = SWveg[i, j] + (Tf(1) - trcn[i, j]) * (LWeff[i, j] + sb * Tsrf[i, j]^Tf(4) - Tf(2) * sb * Tveg[i, j]^Tf(4))

        # Surface energy balance increments without melt
        A_ebal[1, 1] = Tf(0)
        A_ebal[1, 2] = -(KHa[i, j] + KHv[i, j] + KHg[i, j])
        A_ebal[1, 3] = KHg[i, j]
        A_ebal[1, 4] = KHv[i, j]
        b_ebal[1] = (H[i, j] - Hveg - Hsrf[i, j]) / (rho * cp)
        A_ebal[2, 1] = -(KHa[i, j] + KWv[i, j] + KWg[i, j])
        A_ebal[2, 2] = Tf(0)
        A_ebal[2, 3] = Dsrf * KWg[i, j]
        A_ebal[2, 4] = Dveg * KWv[i, j]
        b_ebal[2] = (E - Eveg[i, j] - Esrf[i, j]) / rho
        A_ebal[3, 1] = -Lsrf * rho * KWg[i, j]
        A_ebal[3, 2] = -rho * cp * KHg[i, j]
        A_ebal[3, 3] = rho * (cp * KHg[i, j] + Lsrf * Dsrf * KWg[i, j]) + Tf(4) * sb * Tsrf[i, j]^Tf(3) + Tf(2) * ks1[i, j] / Ds1[i, j]
        A_ebal[3, 4] = -Tf(4) * (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(3)
        b_ebal[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j]
        A_ebal[4, 1] = -Lveg * rho * KWv[i, j]
        A_ebal[4, 2] = -rho * cp * KHv[i, j]
        A_ebal[4, 3] = -Tf(4) * (Tf(1) - trcn[i, j]) * sb * Tsrf[i, j]^Tf(3)
        A_ebal[4, 4] = canh[i, j] / dt + rho * (cp * KHv[i, j] + Lveg * Dveg * KWv[i, j]) + Tf(8) * (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(3)
        b_ebal[4] = Rveg - Hveg - Lveg * Eveg[i, j] - canh[i, j] * (Tveg[i, j] - Tveg0[i, j]) / dt
        ludcmp!(4, A_ebal, Acp_ebal, b_ebal, x_ebal, vv_ebal, indx_ebal)
        dQc = x_ebal[1]
        dTc = x_ebal[2]
        dTs = x_ebal[3]
        dTv = x_ebal[4]
        dEs = rho * KWg[i, j] * (Dsrf * dTs - dQc)
        dEv = rho * KWv[i, j] * (Dveg * dTv - dQc)
        dGs = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
        dHs = rho * cp * KHg[i, j] * (dTs - dTc)
        dHv = rho * cp * KHv[i, j] * (dTv - dTc)

        # Surface melting
        if (Tsrf[i, j] + dTs > Tm && Sice[1, i, j] > eps(Sice[1, i, j]))
            Melt[i, j] = column_sum(Sice, i, j) / dt
            b_ebal[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j] - Lf * Melt[i, j]
            ludcmp!(4, A_ebal, Acp_ebal, b_ebal, x_ebal, vv_ebal, indx_ebal)
            dQc = x_ebal[1]
            dTc = x_ebal[2]
            dTs = x_ebal[3]
            dTv = x_ebal[4]
            dEs = rho * KWg[i, j] * (Dsrf * dTs - dQc)
            dEv = rho * KWv[i, j] * (Dveg * dTv - dQc)
            dGs = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
            dHs = rho * cp * KHg[i, j] * (dTs - dTc)
            dHv = rho * cp * KHv[i, j] * (dTv - dTc)
            if (Tsrf[i, j] + dTs < Tm)
                Qsrf = qsat(Ps[i, j], Tm)
                Esrf[i, j] = rho * KWg[i, j] * (Qsrf - Qcan[i, j])
                G[i, j] = Tf(2) * ks1[i, j] * (Tm - Ts1[i, j]) / Ds1[i, j]
                Hsrf[i, j] = rho * cp * KHg[i, j] * (Tm - Tcan[i, j])
                Rsrf[i, j] = SWsrf[i, j] + trcn[i, j] * (fsky[i, j] * LWeff[i, j] + (Tf(1) - fsky[i, j]) * sb * Ta[i, j]^Tf(4)) - sb * Tm^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
                Rveg = SWveg[i, j] + (Tf(1) - trcn[i, j]) * (LWeff[i, j] + sb * Tm^Tf(4) - Tf(2) * sb * Tveg[i, j]^Tf(4))
                A_ebal[1, 3] = Tf(0)
                b_ebal[1] = (H[i, j] - Hveg - Hsrf[i, j]) / (rho * cp)
                A_ebal[2, 3] = Tf(0)
                b_ebal[2] = (E - Eveg[i, j] - Esrf[i, j]) / rho
                A_ebal[3, 3] = Tf(1)
                b_ebal[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j]
                A_ebal[4, 3] = Tf(0)
                b_ebal[4] = Rveg - Hveg - Lveg * Eveg[i, j] - canh[i, j] * (Tveg[i, j] - Tveg0[i, j]) / dt
                ludcmp!(4, A_ebal, Acp_ebal, b_ebal, x_ebal, vv_ebal, indx_ebal)
                dQc = x_ebal[1]
                dTc = x_ebal[2]
                Melt[i, j] = x_ebal[3] / Lf
                dTv = x_ebal[4]
                dTs = Tm - Tsrf[i, j]
                dEs = Tf(0)
                dEv = rho * KWv[i, j] * (Dveg * dTv - dQc)
                dGs = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
                dHs = Tf(0)
                dHv = rho * cp * KHv[i, j] * (dTv - dTc)
            end
        end

        # Update temperatures and fluxes
        Qcan[i, j] = Qcan[i, j] + dQc
        Tcan[i, j] = Tcan[i, j] + dTc
        Tsrf[i, j] = Tsrf[i, j] + dTs
        Tveg[i, j] = Tveg[i, j] + dTv
        Esrf[i, j] = Esrf[i, j] + dEs
        Eveg[i, j] = Eveg[i, j] + dEv
        Hsrf[i, j] = Hsrf[i, j] + dHs
        Hveg = Hveg + dHv
        G[i, j] = G[i, j] + dGs
        H[i, j] = Hsrf[i, j] + Hveg
        LE[i, j] = Lsrf * Esrf[i, j] + Lveg * Eveg[i, j]
        LEsrf[i, j] = Lsrf * Esrf[i, j]
        Rnet[i, j] = SWsrf[i, j] + SWveg[i, j] + LWeff[i, j] - trcn[i, j] * sb * Tsrf[i, j]^Tf(4) - (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
        LWsci[i, j] = trcn[i, j] * (fsky[i, j] * LWeff[i, j] + (Tf(1) - fsky[i, j]) * sb * Ta[i, j]^Tf(4)) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
        LWveg[i, j] = (Tf(1) - trcn[i, j]) * (LWeff[i, j] + sb * Tsrf[i, j]^Tf(4) - Tf(2) * sb * Tveg[i, j]^Tf(4))

        # Sublimation limited by amount of snow after melt
        Ssub = column_sum(Sice, i, j) - Melt[i, j] * dt
        if (Ssub > eps(Tf) && Esrf[i, j] * dt > Ssub)
            Esrf[i, j] = Ssub / dt
            LEsrf[i, j] = Ls * Esrf[i, j]
            Hsrf[i, j] = Rnet[i, j] - G[i, j] - LEsrf[i, j] - Lf * Melt[i, j]
        end

    end
end
