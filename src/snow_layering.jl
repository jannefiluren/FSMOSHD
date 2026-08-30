"""
    snow_layering!(fsm, meteo, snowdepth0, Sice0, t)

Accumulation of new snow, snow cover fraction update and relayering.

The per-cell physics lives in `snow_layering_kernel!`, a KernelAbstractions
kernel launched over the whole grid (see `ebalsrf!` for the pattern). All
layer scratch is kernel-local `MVector`s sized via `Val(Nsmax)`, so the
routine is thread-safe per cell. The snow cover fraction update is the
device function [`snowcoverfraction_point!`](@ref); the "6:00 am" history
test on `t` is resolved on the host.
"""
function snow_layering!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}, snowdepth0, Sice0, t) where {Tf <: Real, Ti <: Integer}

    (; Nsmax) = fsm.grid
    (; Ds0) = fsm.diag

    # Initialize Ds0
    Ds0 .= Tf(0)

    # Dates cannot cross into kernels: update history of SWE and hs only if
    # they correspond to 6:00am values
    update_hist = 4.5 < hour(t) < 5.5

    backend = get_backend(fsm.state.Tsnow)
    kernel! = snow_layering_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.landuse, fsm.grid, fsm.params, meteo,
        snowdepth0, Sice0, update_hist, fsm.physics.LAYERING, fsm.physics.SNFRAC, Val(Int(Nsmax));
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

# inbounds = true keeps the kernel-local MVector scratch off the heap; see
# the note on soil_kernel! (a raw @inbounds block inside a @kernel body must
# not be used - it corrupts the KernelAbstractions CPU transformation)
@kernel inbounds = true function snow_layering_kernel!(
        state, diag, landuse, grid, params::Parameters{Tf}, meteo,
        snowdepth0, Sice0, update_hist::Bool,
        LAYERING::AbstractLayering{Tf}, SNFRAC::AbstractSnowFraction{Tf}, ::Val{Nsmax},
    ) where {Tf, Nsmax}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    (; tthresh, hfsn, Tsnow_min) = params
    (; tilefrac) = landuse
    (; Ds, Sice, Sliq, Tsnow, histowet, Nsnow, fsnow) = state
    (; Ds0) = diag
    (; Ta) = meteo

    if tilefrac[i, j] >= tthresh

        # Decrease Nsnow if necessary (e.g. after melting)
        while Nsnow[i, j] > 0 && Ds[1, i, j] < eps(Tf)
            if Nsnow[i, j] > 1
                for n in 1:(Nsnow[i, j] - 1)
                    Ds[n, i, j] = Ds[n + 1, i, j]
                    Sice[n, i, j] = Sice[n + 1, i, j]
                    Sliq[n, i, j] = Sliq[n + 1, i, j]
                    Tsnow[n, i, j] = Tsnow[n + 1, i, j]
                    histowet[n, i, j] = histowet[n + 1, i, j]
                end
            end
            Ds[Nsnow[i, j], i, j] = 0
            Sice[Nsnow[i, j], i, j] = 0
            Sliq[Nsnow[i, j], i, j] = 0
            Tsnow[Nsnow[i, j], i, j] = Tm
            histowet[Nsnow[i, j], i, j] = Tf(0)
            Nsnow[i, j] = Nsnow[i, j] - 1
        end

        if LAYERING isa OriginalLayering
            Sice[1, i, j] = Sice[1, i, j] + Sice0[i, j]
        end
        snowdepth = column_sum(Ds, i, j) * fsnow[i, j] + snowdepth0[i, j]

        # Store previous snow cover fraction
        fold = fsnow[i, j]
        # Updated Fractional Snow-Covered Area
        SWEtmp = column_sum(Sice, i, j) + column_sum(Sliq, i, j)
        if LAYERING isa DensityLayering
            SWEtmp = SWEtmp + Sice0[i, j]
        end

        snowcoverfraction_point!(
            SNFRAC, state, landuse, snowdepth, SWEtmp, hfsn, i, j, update_hist
        )

        # Rescale Ds with new snow cover fraction
        if fsnow[i, j] > eps(Tf)
            Ds0[i, j] = snowdepth0[i, j] / fsnow[i, j]
            # Update surface layer thickness based on new fsnow
            if LAYERING isa OriginalLayering
                Ds[1, i, j] = Ds[1, i, j] * fold / fsnow[i, j] + Ds0[i, j]
            else
                Ds[1, i, j] = Ds[1, i, j] * fold / fsnow[i, j]
            end
        else
            Nsnow[i, j] = 0
            for k in 1:Nsmax
                Ds[k, i, j] = 0
                Sice[k, i, j] = 0
                Sliq[k, i, j] = 0
                Tsnow[k, i, j] = Tm
                histowet[k, i, j] = Tf(0)
            end
        end
        if Nsnow[i, j] > 1
            for k in 2:Nsnow[i, j]
                Ds[k, i, j] = Ds[k, i, j] * fold / fsnow[i, j]
            end
        end

        # New snow temperature
        Tsnow0 = min(Ta[i, j], Tm)
        Tsnow0 = max(Tsnow0, Tsnow_min)

        relayer_snow!(LAYERING, i, j, state, diag, grid, params, snowdepth, Tsnow0, Val(Nsmax))

    end
end
