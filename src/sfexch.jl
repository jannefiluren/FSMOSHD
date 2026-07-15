"""
    sfexch!(fsm, meteo)

Surface exchange coefficients and turbulent transfer calculations.

The per-cell physics lives in `sfexch_kernel!`, a KernelAbstractions kernel
launched over the whole grid (see `ebalsrf!` for the pattern).

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions (read-only)
"""
function sfexch!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    @unpack CANMOD, ZOFFST, EXCHNG, SNFRAC = fsm

    @unpack tthresh = fsm

    @unpack zT, zU = fsm

    @unpack Nx, Ny = fsm

    @unpack bstb, cden, cveg, gsnf, rchd, rchz, wcan, zsub, zgf, zgr, khcf, z0_snow = fsm

    @unpack VAI, z0sf = fsm

    @unpack Qcan, fsnow, Sice, Sveg, Tcan, Tsrf, Tveg, Ds = fsm

    @unpack fveg, fves, hcan, tilefrac = fsm

    @unpack KH, KHa, KHg, KHv, KWg, KWv, Usc = fsm

    @unpack gs1 = fsm

    @unpack Qa, Uaeff = fsm

    @unpack Ta, Ps = meteo

    backend = get_backend(Tsrf)
    kernel! = sfexch_kernel!(backend)
    kernel!(
        KH, KHa, KHg, KHv, KWg, KWv, Usc,
        z0_snow, z0sf, VAI, Qcan, fsnow, Sice, Sveg, Tcan, Tsrf, Tveg, Ds,
        fveg, fves, hcan, tilefrac, gs1, Qa, Uaeff, Ta, Ps,
        tthresh, zT, zU, bstb, cden, cveg, gsnf, rchd, rchz, wcan, zsub, zgf, zgr, khcf,
        CANMOD, ZOFFST, EXCHNG, SNFRAC;
        ndrange = (Int(Nx), Int(Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function sfexch_kernel!(
        KH, KHa, KHg, KHv, KWg, KWv, Usc,
        z0_snow, z0sf, VAI, Qcan, fsnow,
        Sice, Sveg, Tcan, Tsrf, Tveg, Ds,
        fveg, fves, hcan, tilefrac, gs1,
        Qa, Uaeff, Ta, Ps,
        tthresh::Tf, zT::Tf, zU::Tf, bstb::Tf, cden::Tf, cveg::Tf, gsnf::Tf,
        rchd::Tf, rchz::Tf, wcan::Tf, zsub::Tf, zgf::Tf, zgr::Tf, khcf::Tf,
        CANMOD::Ti, ZOFFST::Ti, EXCHNG::Ti, SNFRAC::Ti,
    ) where {Tf, Ti}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

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
        z0g = z0_snow[i, j]

        # BC, stabilize the tuning point runs by using a Ds threshold instead of fsnow.
        # TODO: test the impact for the grid points.
        if (SNFRAC == 3)
            sumtmp = column_sum(Ds, i, j)
            if (sumtmp <= Tf(0.05))
                z0g = z0sf[i, j]
            end
        else
            if (fsnow[i, j] <= eps(Tf))
                z0g = z0sf[i, j]
            end
        end

        # Additional roughness lengths and friction velocity
        if (EXCHNG == 2) # Forest - specific adjustment *GM
            # Open
            z0h = Tf(0.1) * z0g
            ustar = vkman * Uaeff[i, j] / log(zU / z0g)
            rgo = log(zT / z0h) / (vkman * ustar)

            # Forest
            if (fveg[i, j] > eps(Tf))
                z0g = (zgf + zgr * fveg[i, j]) * z0g
                z0h = Tf(0.1) * z0g
                dh = rchd * hcan[i, j]
                z0v = rchz * hcan[i, j]
                ustar = vkman * Uaeff[i, j] / log((zU1 - dh) / z0v)
                Uh = (ustar / vkman) * log((hcan[i, j] - dh) / z0v)
                KHh = vkman * ustar * (hcan[i, j] - dh)
                Usf = exp(wcan * (zsub / hcan[i, j] - Tf(1))) * Uh
            end
        else
            z0v = rchz * hcan[i, j]
            z0 = (z0v^fveg[i, j]) * (z0g^(Tf(1) - fveg[i, j]))
            z0h = Tf(0.1) * z0
            dh = fveg[i, j] * rchd * hcan[i, j]
            CD = (vkman / log((zU1 - dh) / z0))^Tf(2)
            ustar = sqrt(CD) * Uaeff[i, j]
        end
        Uso = Uaeff[i, j] * log(zsub / z0g) / log(zU / z0g)

        if (EXCHNG == 0)
            # No stability adjustment
            fh = Tf(1.0)
            Ric = Tf(0.0)
        end
        if (EXCHNG == 1)
            # Stability adjustment (Louis et al. 1982, quoted by Beljaars 1992)
            Tint = fveg[i, j] * Tveg[i, j] + (Tf(1) - fveg[i, j]) * Tsrf[i, j]
            RiB = grav * (Ta[i, j] - Tint) * (zU1 - dh)^Tf(2) / ((zT1 - dh) * Ta[i, j] * Uaeff[i, j]^Tf(2))
            if (RiB > Tf(0.2))
                RiB = Tf(0.2) # New maximum threshold for RiB
            end
            if (RiB > Tf(0))
                fh = Tf(1) / (Tf(1) + Tf(3) * bstb * RiB * sqrt(Tf(1) + bstb * RiB))
            else
                fh = Tf(1) - Tf(3) * bstb * RiB / (Tf(1) + Tf(3) * bstb^Tf(2) * CD * sqrt(-RiB * zU1 / z0))
            end
            Ric = grav * (Tcan[i, j] - Tsrf[i, j]) * hcan[i, j] / (Tcan[i, j] * ustar^Tf(2))
            Ric = max(min(Ric, Tf(10.0)), Tf(0.0))
        end
        # Note that currently, fh and Ric are not used in EXCHNG == 2, i.e. no stability correction

        # Eddy diffusivities
        if (fveg[i, j] == 0)
            KH[i, j] = fh * vkman * ustar / log(zT1 / z0h)
            Qs = qsat(Ps[i, j], Tsrf[i, j])
            if (Sice[1, i, j] > eps(Tf) || Qa[i, j] > Qs)
                KWg[i, j] = KH[i, j]
            else
                KWg[i, j] = gs1[i, j] * KH[i, j] / (gs1[i, j] + KH[i, j])
            end
            Usc[i, j] = Uso
        else
            if (EXCHNG == 2)
                rad = (log((zT1 - dh) / (hcan[i, j] - dh)) / (vkman * ustar) + hcan[i, j] * (exp(wcan * (Tf(1) - (z0v + dh) / hcan[i, j])) - Tf(1)) / (wcan * KHh)) / khcf
                KHa[i, j] = sqrt(fves[i, j]) / rad
                Usub = sqrt(fves[i, j]) * Usf + (Tf(1) - sqrt(fves[i, j])) * Uso
                Usub = max(Usub, Tf(0.1))
                rgd = Tf(1) / (vkman^Tf(2) * Usub) * log(zsub / z0h) * log(zsub / z0g)
                KHg[i, j] = Tf(1) / rgd
                Uc = exp(wcan * ((z0v + dh) / hcan[i, j] - Tf(1))) * Uh
                KHv[i, j] = VAI[i, j] * sqrt(Uc) / cveg
                Usc[i, j] = Usub
            else
                KHa[i, j] = fh * vkman * ustar / log((zT1 - dh) / z0)
                KHg[i, j] = vkman * ustar * ((Tf(1) - fveg[i, j]) * fh / log(z0 / z0h) + fveg[i, j] * cden / (Tf(1) + Tf(0.5) * Ric))
                KHv[i, j] = sqrt(ustar) * VAI[i, j] / cveg
            end
            Qs = qsat(Ps[i, j], Tsrf[i, j])
            if (Qcan[i, j] > Qs)
                KWg[i, j] = KHg[i, j]
            else
                KWg[i, j] = gs1[i, j] * KHg[i, j] / (gs1[i, j] + KHg[i, j])
            end
            Qs = qsat(Ps[i, j], Tveg[i, j])
            if (Sveg[i, j] > eps(Tf) || Qcan[i, j] > Qs)
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
