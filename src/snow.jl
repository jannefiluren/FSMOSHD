"""
    snow!(fsm, meteo, t)

Snow physics processes including heat conduction, melting, sublimation, hydraulics, and compaction.

The per-cell physics lives in `snow_kernel!`, a KernelAbstractions kernel
launched over the whole grid (see `ebalsrf!` for the pattern). Each cell
solves its own tridiagonal heat-conduction system using kernel-local
`MVector` scratch, so the routine is thread-safe per cell. The number of
snow layers is passed as `Val(Nsmax)` because the scratch size must be known
at compile time. Snow layering and the (CPU-only) transport routines are
called afterwards on the host, as before.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions (read-only)
- `t`: Current simulation time
"""
function snow!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}, t) where {Tf <: Real, Ti <: Integer}

    (; dt) = fsm.params
    (; Nsmax) = fsm.grid

    (; Gsoil, meltflux_out, Roff, Roff_bare, Roff_snow) = fsm.diag
    (; G, snowdepth0, Sice0) = fsm.diag
    (; fsnow) = fsm.state

    @unpack Rf = meteo

    Gsoil .= G
    Roff .= Tf(0)
    meltflux_out .= Tf(0)
    Roff_bare .= Rf .* dt .* (Tf(1) .- fsnow)
    Roff_snow .= Rf .* dt .* fsnow

    # Initialize arrays for snow layering
    snowdepth0 .= Tf(0)
    Sice0 .= Tf(0)

    # Points with existing snowpack
    backend = get_backend(fsm.state.Tsnow)
    kernel! = snow_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.landuse, fsm.grid, fsm.params, meteo,
        Val(Int(Nsmax));
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    # Accumulation of new snow, calculation of snow cover fraction and relayering
    snow_layering!(fsm, meteo, snowdepth0, Sice0, t)

    return nothing

end

# inbounds = true keeps the kernel-local MVector scratch off the heap; see
# the note on soil_kernel! (a raw @inbounds block inside a @kernel body must
# not be used - it corrupts the KernelAbstractions CPU transformation)
@kernel inbounds = true function snow_kernel!(
        state, diag, landuse, grid, params::Parameters{Tf, Ti}, meteo,
        ::Val{Nsmax},
    ) where {Tf, Ti, Nsmax}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    (; dt, tthresh, Wirr, rho0, rhob, rhoc, rhof, rhos_min, rcld, rmlt,
       snda, trho, eta0, eta1, a_eta, b_eta, c_eta, rhos_max,
       HYDROL, DENSTY, SNFRAC, FSNRHO, Tsnow_min) = params
    (; Dzsoil) = grid
    (; dem, tilefrac) = landuse
    (; Tsnow, Ds, Sice, Sliq, histowet, rgrn, Nsnow, fsnow, Tsoil, Tsrf) = state
    (; Sbsrf, Roff_bare, Roff_snow, Roff, meltflux_out, Gsoil, Sice0, snowdepth0,
       unload, ksnow, ksoil, G, Melt, Esrf, Uaeff, Sfeff) = diag
    (; Ta) = meteo

    if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        # Kernel-local scratch (one set per grid cell)
        csnow = zero(MVector{Nsmax, Tf})
        Gs = zero(MVector{Nsmax, Tf})
        a = zero(MVector{Nsmax, Tf})
        bsnow = zero(MVector{Nsmax, Tf})
        c = zero(MVector{Nsmax, Tf})
        rhs = zero(MVector{Nsmax, Tf})
        dTssnow = zero(MVector{Nsmax, Tf})
        gammasnow = zero(MVector{Nsmax, Tf})

        Sbsrf[i, j] = Tf(0)

        # Handle snow unloading above freezing
        if (Ta[i, j] >= Tm)
            # Unloading on bare ground fraction is added to runoff
            Roff_bare[i, j] = Roff_bare[i, j] + unload[i, j] * (Tf(1) - fsnow[i, j])
            # Unloading on snow covered fraction is later added to snow liquid water (see hydraulics)
            Roff_snow[i, j] = Roff_snow[i, j] + unload[i, j] * fsnow[i, j]
        end

        if (fsnow[i, j] > eps(Tf)) # This condition should be equivalent to Nsnow[i,j] > 0

            # Except for point case, apply a minimum threshold of 0.1 to fsnow
            # to avoid 'long tails' in SWE due to slowing down depletion rates
            if (SNFRAC == 3)
                fsnow_thres = fsnow[i, j]
            else
                fsnow_thres = min(fsnow[i, j] + Tf(0.25), Tf(1.0))
            end

            # Heat conduction
            for k in 1:Nsnow[i, j]
                csnow[k] = (Sice[k, i, j] * hcap_ice + Sliq[k, i, j] * hcap_wat) / fsnow[i, j]
            end
            if (Nsnow[i, j] == 1)
                Gs[1] = Tf(2) / (Ds[1, i, j] / ksnow[1, i, j] + Dzsoil[1] / ksoil[1, i, j])
                dTssnow[1] = (G[i, j] + Gs[1] * (Tsoil[1, i, j] - Tsnow[1, i, j])) * dt / (csnow[1] + Gs[1] * dt)
            else
                for k in 1:(Nsnow[i, j] - 1)
                    Gs[k] = Tf(2) / (Ds[k, i, j] / ksnow[k, i, j] + Ds[k + 1, i, j] / ksnow[k + 1, i, j])
                end
                a[1] = Tf(0.0)
                bsnow[1] = csnow[1] + Gs[1] * dt
                c[1] = -Gs[1] * dt
                rhs[1] = (G[i, j] - Gs[1] * (Tsnow[1, i, j] - Tsnow[2, i, j])) * dt
                for k in 2:(Nsnow[i, j] - 1)
                    a[k] = c[k - 1]
                    bsnow[k] = csnow[k] + (Gs[k - 1] + Gs[k]) * dt
                    c[k] = -Gs[k] * dt
                    rhs[k] = Gs[k - 1] * (Tsnow[k - 1, i, j] - Tsnow[k, i, j]) * dt + Gs[k] * (Tsnow[k + 1, i, j] - Tsnow[k, i, j]) * dt
                end
                k = Nsnow[i, j]
                Gs[k] = Tf(2) / (Ds[k, i, j] / ksnow[k, i, j] + Dzsoil[1] / ksoil[1, i, j])
                a[k] = c[k - 1]
                bsnow[k] = csnow[k] + (Gs[k - 1] + Gs[k]) * dt
                c[k] = Tf(0)
                rhs[k] = Gs[k - 1] * (Tsnow[k - 1, i, j] - Tsnow[k, i, j]) * dt + Gs[k] * (Tsoil[1, i, j] - Tsnow[k, i, j]) * dt
                tridiag!(dTssnow, Nsnow[i, j], gammasnow, Nsmax, a, bsnow, c, rhs)
            end
            for k in 1:Nsnow[i, j]
                Tsnow[k, i, j] = Tsnow[k, i, j] + dTssnow[k]
                Tsnow[k, i, j] = max(Tsnow[k, i, j], Tsnow_min)
            end
            k = Nsnow[i, j]
            Gsoil[i, j] = Gs[k] * (Tsnow[k, i, j] - Tsoil[1, i, j])

            # Convert melting ice to liquid water
            dSice = Melt[i, j] * fsnow_thres * dt

            meltflux_out[i, j] = dSice

            for k in 1:Nsnow[i, j]
                coldcont = csnow[k] * (Tm - Tsnow[k, i, j])
                if (coldcont < Tf(0))
                    dSice = dSice - fsnow[i, j] * coldcont / Lf
                    Tsnow[k, i, j] = Tm
                end
                if (dSice > eps(Tf))
                    if (dSice > Sice[k, i, j])  # Layer melts completely
                        dSice = dSice - Sice[k, i, j]
                        Ds[k, i, j] = Tf(0)
                        Sliq[k, i, j] = Sliq[k, i, j] + Sice[k, i, j]
                        Sice[k, i, j] = Tf(0)
                    else                       # Layer melts partially
                        Ds[k, i, j] = (Tf(1) - dSice / Sice[k, i, j]) * Ds[k, i, j]
                        Sice[k, i, j] = Sice[k, i, j] - dSice
                        Sliq[k, i, j] = Sliq[k, i, j] + dSice
                        dSice = Tf(0.0)                # Melt exhausted
                    end
                end
            end

            # Remove snow by sublimation
            dSice = max(Esrf[i, j] * fsnow_thres, Tf(0.0)) * dt
            if (dSice > eps(Tf))
                for k in 1:Nsnow[i, j]
                    if (dSice > Sice[k, i, j])  # Layer sublimates completely
                        dSice = dSice - Sice[k, i, j]
                        Ds[k, i, j] = Tf(0)
                        Sbsrf[i, j] = Sbsrf[i, j] + Sice[k, i, j]
                        Sice[k, i, j] = Tf(0)
                    else                       # Layer sublimates partially
                        Ds[k, i, j] = (Tf(1) - dSice / Sice[k, i, j]) * Ds[k, i, j]
                        Sice[k, i, j] = Sice[k, i, j] - dSice
                        Sbsrf[i, j] = Sbsrf[i, j] + dSice
                        dSice = Tf(0.0)                # Sublimation exhausted
                    end
                end
            end

            # Snow hydraulics
            if (HYDROL == 0)
                # Free-draining snow
                meltflux_out[i, j] = Tf(0)
                for k in 1:Nsnow[i, j]
                    Roff_snow[i, j] = Roff_snow[i, j] + Sliq[k, i, j]
                    meltflux_out[i, j] = meltflux_out[i, j] + Sliq[k, i, j]
                    Sliq[k, i, j] = Tf(0)
                end
            elseif (HYDROL == 1)
                # Bucket storage
                for k in 1:Nsnow[i, j]
                    phi = Tf(0.0)
                    if (Ds[k, i, j] > eps(Tf))
                        phi = Tf(1) - Sice[k, i, j] / (rho_ice * Ds[k, i, j] * fsnow[i, j])
                    end
                    SliqMax = fsnow[i, j] * rho_wat * Ds[k, i, j] * phi * Wirr
                    Sliq[k, i, j] = Sliq[k, i, j] + Roff_snow[i, j]
                    Roff_snow[i, j] = Tf(0)
                    if (Sliq[k, i, j] > SliqMax)       # Liquid capacity exceeded
                        Roff_snow[i, j] = Sliq[k, i, j] - SliqMax   # so drainage to next layer
                        Sliq[k, i, j] = SliqMax
                        histowet[k, i, j] = Tf(1.0)
                    end
                    # csnow needs to be updated after changing Sliq and Sice
                    csnow[k] = (Sice[k, i, j] * hcap_ice + Sliq[k, i, j] * hcap_wat) / fsnow[i, j]
                    coldcont = csnow[k] * (Tm - Tsnow[k, i, j])
                    if (coldcont > Tf(0))       # Liquid can freeze
                        dSice = min(Sliq[k, i, j], fsnow[i, j] * coldcont / Lf)
                        Sliq[k, i, j] = Sliq[k, i, j] - dSice
                        Sice[k, i, j] = Sice[k, i, j] + dSice
                        meltflux_out[i, j] = meltflux_out[i, j] - dSice
                        Tsnow[k, i, j] = Tsnow[k, i, j] + Lf * dSice / csnow[k] / fsnow[i, j]
                    end
                end

                if (meltflux_out[i, j] < Tf(0))
                    meltflux_out[i, j] = Tf(0)
                end

            else  # HYDROL == 2
                # Density-dependent bucket storage
                for k in 1:Nsnow[i, j]
                    SliqCap = Tf(0.0)
                    if (Ds[k, i, j] > eps(Tf))
                        rhos = Sice[k, i, j] / Ds[k, i, j] / fsnow[i, j]
                        SliqCap = Tf(0.03) + Tf(0.07) * (Tf(1) - rhos / Tf(200))
                        SliqCap = max(SliqCap, Tf(0.03))
                    end
                    SliqMax = SliqCap * Sice[k, i, j]
                    Sliq[k, i, j] = Sliq[k, i, j] + Roff_snow[i, j]
                    Roff_snow[i, j] = Tf(0)
                    if (Sliq[k, i, j] > SliqMax)       # Liquid capacity exceeded
                        Roff_snow[i, j] = Sliq[k, i, j] - SliqMax   # so drainage to next layer
                        Sliq[k, i, j] = SliqMax
                        histowet[k, i, j] = Tf(1.0)
                    end
                    # csnow needs to be updated after changing Sliq and Sice
                    csnow[k] = (Sice[k, i, j] * hcap_ice + Sliq[k, i, j] * hcap_wat) / fsnow[i, j]
                    coldcont = csnow[k] * (Tm - Tsnow[k, i, j])
                    if (coldcont > eps(Tf))       # Liquid can freeze
                        dSice = min(Sliq[k, i, j], fsnow[i, j] * coldcont / Lf)
                        Sliq[k, i, j] = Sliq[k, i, j] - dSice
                        Sice[k, i, j] = Sice[k, i, j] + dSice
                        # to account for refreezing of melt
                        meltflux_out[i, j] = meltflux_out[i, j] - dSice
                        Tsnow[k, i, j] = Tsnow[k, i, j] + Lf * dSice / csnow[k] / fsnow[i, j]
                    end
                end

                if (meltflux_out[i, j] < Tf(0))
                    meltflux_out[i, j] = Tf(0)
                end

            end

            # Snow compaction
            if (DENSTY == 0)
                # Fixed snow density
                for k in 1:Nsnow[i, j]
                    Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rho0 / fsnow[i, j]
                end
            elseif (DENSTY == 1)
                # Snow compaction with age
                for k in 1:Nsnow[i, j]
                    if (Ds[k, i, j] > eps(Tf))
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
                mass = Tf(0.0)
                for k in 1:Nsnow[i, j]
                    mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
                    if (Ds[k, i, j] > eps(Tf))
                        rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
                        rhos = rhos + (rhos * grav * mass * dt / (eta0 * exp(-(Tsnow[k, i, j] - Tm) / Tf(12.4) + rhos / Tf(55.6))) + dt * rhos * snda * exp((Tsnow[k, i, j] - Tm) / Tf(23.8) - max(rhos - Tf(150), Tf(0.0)) / Tf(21.7)))
                        rhos = min(rhos, rhos_max)
                        Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rhos / fsnow[i, j]
                    end
                    mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
                end
            else # DENSTY == 3
                # Snow compaction by overburden, dependent on liquid water content (Crocus B92)
                mass = Tf(0.0)
                for k in 1:Nsnow[i, j]
                    mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
                    if (Ds[k, i, j] > eps(Tf))
                        rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
                        f1 = Tf(1) / (Tf(1) + Tf(600) * Sliq[k, i, j] / (rho_wat * Ds[k, i, j] * fsnow[i, j]))
                        f2 = Tf(1.0)
                        eta = f1 * f2 * eta1 * (rhos / c_eta) * exp(a_eta * (Tm - Tsnow[k, i, j]) + b_eta * rhos)
                        rhos = rhos + rhos * grav * mass * dt / eta
                        rhos = min(rhos, rhos_max)
                        Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rhos / fsnow[i, j]
                    end
                    mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
                end
            end

            # Snow grain growth --> *GM for now, this code feature is not functional because the state variable rgrn is not tracked (no bin output)
            for k in 1:Nsnow[i, j]
                ggr = Tf(2.0e-13)
                if (Tsnow[k, i, j] < Tm)
                    if (rgrn[k, i, j] < Tf(1.5e-4))
                        ggr = Tf(2.0e-14)
                    else
                        ggr = Tf(7.3e-8) * exp(Tf(-4600) / Tsnow[k, i, j])
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

        # Add snowfall and frost to new snow with fresh snow density and grain size
        Esnow = Tf(0.0)
        if (Esrf[i, j] < Tf(0) && Tsrf[i, j] < Tm)
            Esnow = fsnow[i, j] * Esrf[i, j]
            Sbsrf[i, j] = Esnow * dt
        end
        dSice = (Sfeff[i, j] - Esnow) * dt  # Think about how to scale for fsnow...

        # Catch to round infinitesimally small new snow amounts.
        # The small amounts were due to EnKF-assimilated daily precip being downscaled to hourly
        # based on the hourly (COSMO) to daily (COSMO) mass ratio.  Since min daily EnKF was 1mm,
        # when 1mm was downscaled to hourly this mass could become
        # quite small.  When this occurred on bare ground, Tsnow calculation would blow up,
        # place NaN in Tsnow and snow would never again melt through the season
        # Do we still need this Catch in FSM??
        if (Nsnow[i, j] <= 1 && dSice < Tf(0.001) && Sice[1, i, j] < Tf(0.001))
            dSice = Tf(trunc(Ti, dSice * Tf(1000) + Tf(0.5))) / Tf(1000)    # TODO verify against original code
        end

        rhonew = fresh_snow_density(FSNRHO, rho0, rhob, rhoc, rhof, rhos_min, Ta[i, j], Uaeff[i, j], dem[i, j])

        Sice0[i, j] = dSice
        snowdepth0[i, j] = dSice / rhonew
        # Add canopy unloading to new snow with bulk snow density
        rhos = rhof
        if (Ta[i, j] < Tm)  # only if it's cold enough
            mass = column_sum(Sice, i, j) + column_sum(Sliq, i, j)
            snowdepth = column_sum(Ds, i, j) * fsnow[i, j]
            if (snowdepth > eps(Tf))
                rhos = mass / snowdepth
            end
            Sice0[i, j] = Sice0[i, j] + unload[i, j]
            snowdepth0[i, j] = snowdepth0[i, j] + unload[i, j] / rhos
        end

    end
end
