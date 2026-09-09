# Surface energy balance

"""
    energy_balance!(canopy, i, j, state, diag, surface, params, meteo, substrate)

Solve the energy balance at cell `(i, j)`, implemented for every `AbstractCanopy`:
`NoCanopy` solves the surface alone, `OneLayerCanopy` solves the joint surface and canopy
system. `substrate` caps the surface temperature over bare glacier ice and is
unused under a canopy.
"""
function energy_balance! end

"""
    surface_energy_balance!(fsm, meteo)

Surface energy balance solution, coupled to the canopy where there is one.

# Arguments
- `fsm::FSM`: Model state structure
- `meteo::MET`: Current meteorological conditions
"""
function surface_energy_balance!(fsm::FSM{Tf}, meteo::MET{Tf}) where {Tf <: Real}

    (; canopy, substrate) = fsm.physics

    backend = get_backend(fsm.state.Tsrf)
    kernel! = surface_energy_balance_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.params, meteo,
        canopy, substrate;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

# inbounds = true keeps the kernel-local scratch of the forest solver
# off the heap; raw @inbounds block inside a @kernel body corrupts the
# KernelAbstractions CPU transformation
@kernel inbounds = true function surface_energy_balance_kernel!(
        state, diag, surface, params::Parameters{Tf}, meteo,
        canopy::AbstractCanopy{Tf}, substrate::AbstractSubstrate{Tf},
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; tthresh) = params
    (; tilefrac) = surface

    if (tilefrac[i, j] >= tthresh)

        energy_balance!(canopy, i, j, state, diag, surface, params, meteo, substrate)

    end
end

# Open and non-forest tiles
@inline function energy_balance!(::NoCanopy{Tf}, i, j, state, diag, surface, params, meteo, substrate) where {Tf}

    @unpack_constants(Tf)

    (; dt) = params
    (; trcn) = surface
    (; Sice, Tcan, Tsrf, Tveg) = state
    (; Esrf, G, H, Hsrf, LE, LEsrf, LWsci, LWveg, Melt, Rnet, Rsrf, SWsrf, Ds1, Ts1, ks1, KH, KWg, Qa, LWeff) = diag
    (; Ps, Ta) = meteo

    # Reported as air temperature so that Tcan and Tveg are defined in open runs
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

    # Handle energy fluxes on bare glacier ice assuming glaciers are an
    # infinite heat reservoir. The capping of surface temperature at melting
    # point is not energy conserving
    if substrate isa IceSubstrate
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

    return nothing
end

# Forest tiles
@inline function energy_balance!(::OneLayerCanopy{Tf}, i, j, state, diag, surface, params, meteo, substrate) where {Tf}

    @unpack_constants(Tf)

    (; dt) = params
    (; canh, fsky, trcn) = surface
    (; Qcan, Sice, Tcan, Tsrf, Tveg) = state
    (; Esrf, Eveg, G, H, Hsrf, LE, LEsrf, LWsci, LWveg, Melt, Rnet, Rsrf, Ds1, KHa, KHg, KHv, KWg, KWv, ks1, SWsrf, SWveg, Ts1, Tveg0, Qa, LWeff) = diag
    (; Ps, Ta) = meteo

    # Kernel-local scratch variables
    A = zero(MMatrix{4, 4, Tf})
    Acp = zero(MMatrix{4, 4, Tf})
    b = zero(MVector{4, Tf})
    x = zero(MVector{4, Tf})
    vv = zero(MVector{4, Tf})
    indx = zero(MVector{4, Int32})

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
    Rveg = SWveg[i, j] + (Tf(1) - trcn[i, j]) * (LWeff[i, j] + sb * Tsrf[i, j]^Tf(4) - Tf(2) * sb * Tveg[i, j]^Tf(4))

    # Surface energy balance increments without melt
    A[1, 1] = Tf(0)
    A[1, 2] = -(KHa[i, j] + KHv[i, j] + KHg[i, j])
    A[1, 3] = KHg[i, j]
    A[1, 4] = KHv[i, j]
    b[1] = (H[i, j] - Hveg - Hsrf[i, j]) / (rho * cp)
    A[2, 1] = -(KHa[i, j] + KWv[i, j] + KWg[i, j])
    A[2, 2] = Tf(0)
    A[2, 3] = Dsrf * KWg[i, j]
    A[2, 4] = Dveg * KWv[i, j]
    b[2] = (E - Eveg[i, j] - Esrf[i, j]) / rho
    A[3, 1] = -Lsrf * rho * KWg[i, j]
    A[3, 2] = -rho * cp * KHg[i, j]
    A[3, 3] = rho * (cp * KHg[i, j] + Lsrf * Dsrf * KWg[i, j]) + Tf(4) * sb * Tsrf[i, j]^Tf(3) + Tf(2) * ks1[i, j] / Ds1[i, j]
    A[3, 4] = -Tf(4) * (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(3)
    b[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j]
    A[4, 1] = -Lveg * rho * KWv[i, j]
    A[4, 2] = -rho * cp * KHv[i, j]
    A[4, 3] = -Tf(4) * (Tf(1) - trcn[i, j]) * sb * Tsrf[i, j]^Tf(3)
    A[4, 4] = canh[i, j] / dt + rho * (cp * KHv[i, j] + Lveg * Dveg * KWv[i, j]) + Tf(8) * (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(3)
    b[4] = Rveg - Hveg - Lveg * Eveg[i, j] - canh[i, j] * (Tveg[i, j] - Tveg0[i, j]) / dt
    ludcmp!(4, A, Acp, b, x, vv, indx)
    dQc = x[1]
    dTc = x[2]
    dTs = x[3]
    dTv = x[4]
    dEs = rho * KWg[i, j] * (Dsrf * dTs - dQc)
    dEv = rho * KWv[i, j] * (Dveg * dTv - dQc)
    dGs = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
    dHs = rho * cp * KHg[i, j] * (dTs - dTc)
    dHv = rho * cp * KHv[i, j] * (dTv - dTc)

    # Surface melting
    if (Tsrf[i, j] + dTs > Tm && Sice[1, i, j] > eps(Sice[1, i, j]))
        Melt[i, j] = column_sum(Sice, i, j) / dt
        b[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j] - Lf * Melt[i, j]
        ludcmp!(4, A, Acp, b, x, vv, indx)
        dQc = x[1]
        dTc = x[2]
        dTs = x[3]
        dTv = x[4]
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
            A[1, 3] = Tf(0)
            b[1] = (H[i, j] - Hveg - Hsrf[i, j]) / (rho * cp)
            A[2, 3] = Tf(0)
            b[2] = (E - Eveg[i, j] - Esrf[i, j]) / rho
            A[3, 3] = Tf(1)
            b[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j]
            A[4, 3] = Tf(0)
            b[4] = Rveg - Hveg - Lveg * Eveg[i, j] - canh[i, j] * (Tveg[i, j] - Tveg0[i, j]) / dt
            ludcmp!(4, A, Acp, b, x, vv, indx)
            dQc = x[1]
            dTc = x[2]
            Melt[i, j] = x[3] / Lf
            dTv = x[4]
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

    return nothing
end
