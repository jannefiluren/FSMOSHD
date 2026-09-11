"""
    snow!(fsm, meteo, t)

Snow physics processes including heat conduction, melting, sublimation, hydraulics, and compaction.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions (read-only)
- `t`: Current simulation time
"""
function snow!(fsm::FSM{Tf}, meteo::MET{Tf}, t) where {Tf <: Real}

    (; Nsmax) = fsm.grid

    # Resolve dates here as they cannot cross into kernels
    update_hist = 4.5 < hour(t) < 5.5

    backend = get_backend(fsm.state.Tsnow)
    kernel! = snow_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.grid, fsm.params, meteo,
        fsm.physics.new_snow_density, fsm.physics.compaction, fsm.physics.hydrology, fsm.physics.snow_fraction,
        fsm.physics.layering, update_hist, Val(Int(Nsmax));
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing

end

@kernel inbounds = true function snow_kernel!(
        state, diag, surface, grid, params::Parameters{Tf}, meteo,
        new_snow_density::AbstractFreshSnowDensity{Tf}, compaction::AbstractCompaction{Tf},
        hydrology::AbstractHydrology{Tf}, snow_fraction::AbstractSnowFraction{Tf},
        layering::AbstractLayering{Tf}, update_hist::Bool, ::Val{Nsmax},
    ) where {Tf, Nsmax}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    (; dt, tthresh, rho0, rhob, rhoc, rhof, rhos_min, Tsnow_min) = params
    (; Dzsoil) = grid
    (; dem, tilefrac) = surface
    (; Tsnow, Ds, Sice, Sliq, Nsnow, fsnow, Tsoil, Tsrf) = state
    (;
        Sbsrf, Roff_bare, Roff_snow, Roff, meltflux_out, Gsoil, Sice0, snowdepth0,
        unload, ksnow, ksoil, G, Melt, Esrf, Uaeff, Sfeff,
    ) = diag
    (; Rf, Ta) = meteo

    if (tilefrac[i, j] >= tthresh)

        # Kernel-local scratch
        csnow = zero(MVector{Nsmax, Tf})
        Gs = zero(MVector{Nsmax, Tf})
        a = zero(MVector{Nsmax, Tf})
        b = zero(MVector{Nsmax, Tf})
        c = zero(MVector{Nsmax, Tf})
        rhs = zero(MVector{Nsmax, Tf})
        dTs = zero(MVector{Nsmax, Tf})

        # Accumulators for this step
        Sbsrf[i, j] = Tf(0)
        Gsoil[i, j] = G[i, j]
        Roff[i, j] = Tf(0)
        meltflux_out[i, j] = Tf(0)
        Roff_bare[i, j] = Rf[i, j] * dt * (Tf(1) - fsnow[i, j])
        Roff_snow[i, j] = Rf[i, j] * dt * fsnow[i, j]
        snowdepth0[i, j] = Tf(0)
        Sice0[i, j] = Tf(0)

        # Handle snow unloading above freezing
        if (Ta[i, j] >= Tm)
            # Unloading on bare ground fraction is added to runoff
            Roff_bare[i, j] = Roff_bare[i, j] + unload[i, j] * (Tf(1) - fsnow[i, j])
            # Unloading on snow-covered fraction is added to snow liquid water later
            Roff_snow[i, j] = Roff_snow[i, j] + unload[i, j] * fsnow[i, j]
        end

        if (fsnow[i, j] > eps(Tf)) # This condition should be equivalent to Nsnow[i,j] > 0

            fsnow_thres = melt_snow_fraction(snow_fraction, i, j, state)

            # Heat conduction
            for k in 1:Nsnow[i, j]
                csnow[k] = (Sice[k, i, j] * hcap_ice + Sliq[k, i, j] * hcap_wat) / fsnow[i, j]
            end
            if (Nsnow[i, j] == 1)
                Gs[1] = Tf(2) / (Ds[1, i, j] / ksnow[1, i, j] + Dzsoil[1] / ksoil[1, i, j])
                dTs[1] = (G[i, j] + Gs[1] * (Tsoil[1, i, j] - Tsnow[1, i, j])) * dt / (csnow[1] + Gs[1] * dt)
            else
                for k in 1:(Nsnow[i, j] - 1)
                    Gs[k] = Tf(2) / (Ds[k, i, j] / ksnow[k, i, j] + Ds[k + 1, i, j] / ksnow[k + 1, i, j])
                end
                a[1] = Tf(0.0)
                b[1] = csnow[1] + Gs[1] * dt
                c[1] = -Gs[1] * dt
                rhs[1] = (G[i, j] - Gs[1] * (Tsnow[1, i, j] - Tsnow[2, i, j])) * dt
                for k in 2:(Nsnow[i, j] - 1)
                    a[k] = c[k - 1]
                    b[k] = csnow[k] + (Gs[k - 1] + Gs[k]) * dt
                    c[k] = -Gs[k] * dt
                    rhs[k] = Gs[k - 1] * (Tsnow[k - 1, i, j] - Tsnow[k, i, j]) * dt + Gs[k] * (Tsnow[k + 1, i, j] - Tsnow[k, i, j]) * dt
                end
                k = Nsnow[i, j]
                Gs[k] = Tf(2) / (Ds[k, i, j] / ksnow[k, i, j] + Dzsoil[1] / ksoil[1, i, j])
                a[k] = c[k - 1]
                b[k] = csnow[k] + (Gs[k - 1] + Gs[k]) * dt
                c[k] = Tf(0)
                rhs[k] = Gs[k - 1] * (Tsnow[k - 1, i, j] - Tsnow[k, i, j]) * dt + Gs[k] * (Tsoil[1, i, j] - Tsnow[k, i, j]) * dt
                tridiag!(dTs, Nsnow[i, j], a, b, c, rhs)
            end
            for k in 1:Nsnow[i, j]
                Tsnow[k, i, j] = Tsnow[k, i, j] + dTs[k]
                Tsnow[k, i, j] = max(Tsnow[k, i, j], Tsnow_min)
            end
            k = Nsnow[i, j]
            Gsoil[i, j] = Gs[k] * (Tsnow[k, i, j] - Tsoil[1, i, j])

            # Convert melting ice to liquid water
            dSice = Melt[i, j] * fsnow_thres * dt

            meltflux_out[i, j] = dSice

            # Update layers due to refreezing and melt
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
                    else                        # Layer melts partially
                        Ds[k, i, j] = (Tf(1) - dSice / Sice[k, i, j]) * Ds[k, i, j]
                        Sice[k, i, j] = Sice[k, i, j] - dSice
                        Sliq[k, i, j] = Sliq[k, i, j] + dSice
                        dSice = Tf(0.0)         # Melt exhausted
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
                    else                        # Layer sublimates partially
                        Ds[k, i, j] = (Tf(1) - dSice / Sice[k, i, j]) * Ds[k, i, j]
                        Sice[k, i, j] = Sice[k, i, j] - dSice
                        Sbsrf[i, j] = Sbsrf[i, j] + dSice
                        dSice = Tf(0.0)         # Sublimation exhausted
                    end
                end
            end

            snow_hydrology!(hydrology, i, j, state, diag, params)

            compact_snow!(compaction, i, j, state, params)

        end  # Existing snowpack

        # Cap meltflux_out at snow runoff: it ignores liquid retention during percolation
        if (meltflux_out[i, j] > Roff_snow[i, j])
            meltflux_out[i, j] = Roff_snow[i, j]
        end

        # Total runoff
        Roff[i, j] = Roff_snow[i, j] + Roff_bare[i, j]

        # Add snowfall and frost to new snow with fresh snow density
        Esnow = Tf(0.0)
        if (Esrf[i, j] < Tf(0) && Tsrf[i, j] < Tm)
            Esnow = fsnow[i, j] * Esrf[i, j]
            Sbsrf[i, j] = Esnow * dt
        end
        dSice = (Sfeff[i, j] - Esnow) * dt

        # Catch to round infinitesimally small new snow amounts
        if (Nsnow[i, j] <= 1 && dSice < Tf(0.001) && Sice[1, i, j] < Tf(0.001))
            dSice = Tf(trunc(Int, dSice * Tf(1000) + Tf(0.5))) / Tf(1000)
        end

        rhonew = fresh_snow_density(new_snow_density, rho0, rhob, rhoc, rhof, rhos_min, Ta[i, j], Uaeff[i, j], dem[i, j])

        Sice0[i, j] = dSice
        snowdepth0[i, j] = dSice / rhonew

        # Add canopy unloading to new snow with bulk snow density
        rhos = rhof
        if (Ta[i, j] < Tm)
            mass = column_sum(Sice, i, j) + column_sum(Sliq, i, j)
            snowdepth = column_sum(Ds, i, j) * fsnow[i, j]
            if (snowdepth > eps(Tf))
                rhos = mass / snowdepth
            end
            Sice0[i, j] = Sice0[i, j] + unload[i, j]
            snowdepth0[i, j] = snowdepth0[i, j] + unload[i, j] / rhos
        end

        # Accumulation of new snow, snow cover fraction and relayering
        snow_layering!(
            layering, snow_fraction, i, j, state, diag, surface, grid, params, meteo,
            update_hist, Val(Nsmax)
        )

    end
end
