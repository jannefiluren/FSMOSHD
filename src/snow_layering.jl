"""
    snow_layering!(layering, snowfraction, i, j, state, diag, landuse, grid, params, meteo, update_hist, ::Val{Nsmax})

Accumulation of new snow, snow cover fraction update and relayering at cell `(i, j)`,
after the melt, sublimation and compaction of the same step. The snow cover fraction
update is [`snowcoverfraction_point!`](@ref); `update_hist` refreshes the 14-day
history state and is resolved by the caller, since `Dates` cannot run in a kernel.
"""
# @propagate_inbounds: this is the link that carries the kernel's inbounds context
# down to relayer_snow!, whose MVector scratch would otherwise go to the heap
Base.@propagate_inbounds function snow_layering!(
        LAYERING::AbstractLayering{Tf}, SNFRAC::AbstractSnowFraction{Tf},
        i, j, state, diag, landuse, grid, params, meteo, update_hist::Bool, ::Val{Nsmax},
    ) where {Tf, Nsmax}

    @unpack_constants(Tf)

    (; hfsn, Tsnow_min) = params
    (; Ds, Sice, Sliq, Tsnow, histowet, Nsnow, fsnow) = state
    (; Ds0, Sice0, snowdepth0) = diag
    (; Ta) = meteo

    Ds0[i, j] = Tf(0)

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

    return nothing
end
