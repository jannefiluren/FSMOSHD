# Snow relayering parameterizations

struct OriginalLayering{Tf} <: AbstractLayering{Tf} end
struct DensityLayering{Tf} <: AbstractLayering{Tf} end

OriginalLayering{Tf}(grid::Grid; kwargs...) where {Tf} = OriginalLayering{Tf}()
DensityLayering{Tf}(grid::Grid; kwargs...) where {Tf} = DensityLayering{Tf}()

"""
    relayer_snow!(scheme, i, j, state, diag, grid, params, snowdepth, Tsnow0, ::Val{Nsmax})

Rebuild the snow layer structure at cell `(i, j)` in place after mass changes:
repartition `Ds`/`Sice`/`Sliq`/`Tsnow`/`histowet`/`Nsnow` into the scheme's
target layering, conserving mass and internal energy. `snowdepth` and `Tsnow0`
are the (scalar) total snow depth and fresh-snow temperature computed by the
caller. Every `AbstractLayering` implements it.
"""
function relayer_snow! end

# Both methods are @propagate_inbounds, not @inline: without the kernel's inbounds
# context the bounds-check paths capture their MVector scratch onto the heap, once
# per grid cell

# Original layering routine
Base.@propagate_inbounds function relayer_snow!(::OriginalLayering{Tf}, i, j, state, diag, grid, params, snowdepth, Tsnow0, ::Val{Nsmax}) where {Tf, Nsmax}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, Tsnow, Nsnow, fsnow) = state
    (; Sice0) = diag
    (; Dzsnow) = grid
    (; Tsnow_min) = params

    # Kernel-local scratch
    csnow = zero(MVector{Nsmax, Tf})
    D = zero(MVector{Nsmax, Tf})
    E = zero(MVector{Nsmax, Tf})
    S = zero(MVector{Nsmax, Tf})
    U = zero(MVector{Nsmax, Tf})
    W = zero(MVector{Nsmax, Tf})

    # New snowpack
    if (Nsnow[i, j] == 0 && Sice[1, i, j] > eps(Tf))
        Nsnow[i, j] = 1
        Tsnow[1, i, j] = Tsnow0
    end

    # Store state of old layers
    for si in 1:Nsmax
        D[si] = Ds[si, i, j]
        S[si] = Sice[si, i, j]
        W[si] = Sliq[si, i, j]
    end
    if (fsnow[i, j] > eps(Tf))
        csnow[1] = (Sice[1, i, j] * hcap_ice + Sliq[1, i, j] * hcap_wat) / fsnow[i, j]
        E[1] = csnow[1] * (Tsnow[1, i, j] - Tm) + (Sice0[i, j] * hcap_ice / fsnow[i, j]) * (Tsnow0 - Tsnow[1, i, j]) # Adjustment given that csnow[1] already includes the new snow
    else
        fill!(E, Tf(0))
    end
    if (Nsnow[i, j] > 1)
        for k in 2:Nsnow[i, j]
            csnow[k] = (Sice[k, i, j] * hcap_ice + Sliq[k, i, j] * hcap_wat) / fsnow[i, j]
            E[k] = csnow[k] * (Tsnow[k, i, j] - Tm)
        end
    end
    Nold = Nsnow[i, j]

    # Initialise new layers
    for k in 1:Nsmax
        Ds[k, i, j] = Tf(0)
        Sice[k, i, j] = Tf(0)
        Sliq[k, i, j] = Tf(0)
        Tsnow[k, i, j] = Tm
    end
    fill!(U, Tf(0))
    Nsnow[i, j] = 0

    if (fsnow[i, j] > eps(Tf))  # Existing or new snowpack

        # Re-assign and count snow layers
        dnew = snowdepth / fsnow[i, j]
        Ds[1, i, j] = dnew
        if (Ds[1, i, j] > Dzsnow[1])
            for k in 1:Nsmax
                Ds[k, i, j] = Dzsnow[k]
                dnew = dnew - Dzsnow[k]
                if (dnew <= Dzsnow[k] || k == Nsmax)
                    Ds[k, i, j] = Ds[k, i, j] + dnew
                    break
                end
            end
        end
        Nsnow[i, j] = 0
        for si in 1:Nsmax
            if Ds[si, i, j] > Tf(0)
                Nsnow[i, j] += 1
            end
        end

        # Fill new layers from the top downwards
        knew = 1
        dnew = Ds[1, i, j]
        for kold in 1:Nold
            while true
                if (D[kold] < dnew)
                    # All snow from old layer partially fills new layer
                    Sice[knew, i, j] = Sice[knew, i, j] + S[kold]
                    Sliq[knew, i, j] = Sliq[knew, i, j] + W[kold]
                    U[knew] = U[knew] + E[kold]
                    dnew = dnew - D[kold]
                    break
                else
                    # Some snow from old layer fills new layer
                    wt = dnew / D[kold]
                    Sice[knew, i, j] = Sice[knew, i, j] + wt * S[kold]
                    Sliq[knew, i, j] = Sliq[knew, i, j] + wt * W[kold]
                    U[knew] = U[knew] + wt * E[kold]
                    D[kold] = (1 - wt) * D[kold]
                    E[kold] = (1 - wt) * E[kold]
                    S[kold] = (1 - wt) * S[kold]
                    W[kold] = (1 - wt) * W[kold]
                    knew = knew + 1
                    if (knew > Nsnow[i, j])
                        break
                    end
                    dnew = Ds[knew, i, j]
                end
            end
        end

        # Diagnose snow layer temperatures
        for k in 1:Nsnow[i, j]
            csnow[k] = (Sice[k, i, j] * hcap_ice + Sliq[k, i, j] * hcap_wat) / fsnow[i, j]
            Tsnow[k, i, j] = Tm + U[k] / csnow[k]
            Tsnow[k, i, j] = max(Tsnow[k, i, j], Tsnow_min)
        end
    end # Existing or new snowpack

    return nothing
end

# Density dependent snowpack layering
Base.@propagate_inbounds function relayer_snow!(::DensityLayering{Tf}, i, j, state, diag, grid, params, snowdepth, Tsnow0, ::Val{Nsmax}) where {Tf, Nsmax}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, Tsnow, histowet, Nsnow, fsnow) = state
    (; Ds0, Sice0) = diag
    (; rho0, Ds_surflay, Ds_min) = params

    # Kernel-local scratch
    rho = zero(MVector{Nsmax + 1, Tf})
    diff_rho = zero(MVector{Nsmax, Tf})
    csnow_loc = zero(MVector{Nsmax + 1, Tf})
    Sice_loc = zero(MVector{Nsmax + 1, Tf})
    Sliq_loc = zero(MVector{Nsmax + 1, Tf})
    Ds_loc = zero(MVector{Nsmax + 1, Tf})
    histowet_loc = zero(MVector{Nsmax + 1, Tf})
    U_loc = zero(MVector{Nsmax + 1, Tf})
    Tsnow_loc = zero(MVector{Nsmax + 1, Tf})

    # Compute total snow thickness (including new snow if present)
    if Ds0[i, j] > eps(Tf)
        snowthickness = Ds0[i, j] + column_sum(Ds, i, j)
    else
        snowthickness = column_sum(Ds, i, j)
    end
    fill!(rho, rho0)

    # Step 0: Save state variables in local variables that can be up to Nsmax+1
    for k in 1:Nsmax
        Sice_loc[k] = Sice[k, i, j]
        Sliq_loc[k] = Sliq[k, i, j]
        Ds_loc[k] = Ds[k, i, j]
        histowet_loc[k] = histowet[k, i, j]
        Tsnow_loc[k] = Tsnow[k, i, j]
    end
    Sice_loc[Nsmax + 1] = 0
    Sliq_loc[Nsmax + 1] = 0
    Ds_loc[Nsmax + 1] = 0
    histowet_loc[Nsmax + 1] = 0
    Tsnow_loc[Nsmax + 1] = Tm
    Nsnow_loc = Nsnow[i, j]
    if fsnow[i, j] > eps(Tf)
        for k in 1:Nsmax
            csnow_loc[k] = (Sice_loc[k] * hcap_ice + Sliq_loc[k] * hcap_wat) / fsnow[i, j]
            U_loc[k] = csnow_loc[k] * (Tsnow_loc[k] - Tm)
        end
        U_loc[Nsmax + 1] = 0
    else
        fill!(U_loc, Tf(0))
    end

    # Step 1: If there is fresh snow, add the top fresh snow layer and shift layer numbers
    if Ds0[i, j] > eps(Tf)
        Nsnow_loc = Nsnow_loc + 1
        if Nsnow_loc > 1
            for k in 1:(Nsnow_loc - 1)
                Ds_loc[Nsnow_loc - k + 1] = Ds_loc[Nsnow_loc - k]
                Sice_loc[Nsnow_loc - k + 1] = Sice_loc[Nsnow_loc - k]
                Sliq_loc[Nsnow_loc - k + 1] = Sliq_loc[Nsnow_loc - k]
                U_loc[Nsnow_loc - k + 1] = U_loc[Nsnow_loc - k]
                histowet_loc[Nsnow_loc - k + 1] = histowet_loc[Nsnow_loc - k]
            end
        end
        Ds_loc[1] = Ds0[i, j]             # Set new top layer thickness
        Sice_loc[1] = Sice0[i, j]         # Set new top layer ice content
        Sliq_loc[1] = 0                   # No liquid water in new snow
        Tsnow_loc[1] = Tsnow0             # Set new snow temperature
        if fsnow[i, j] > eps(Tf)
            csnow_loc[1] = (Sice_loc[1] * hcap_ice + Sliq_loc[1] * hcap_wat) / fsnow[i, j]
            U_loc[1] = csnow_loc[1] * (Tsnow_loc[1] - Tm)
        else
            U_loc[1] = 0
        end
        histowet_loc[1] = 0               # New snow has never been wet
    end

    # Step 2: Initialise new layers for the case of no snow
    if snowthickness < eps(Tf)
        fill!(Ds_loc, Tf(0))
        fill!(Sice_loc, Tf(0))
        fill!(Sliq_loc, Tf(0))
        Nsnow_loc = zero(Nsnow_loc)
        fill!(U_loc, Tf(0))
        fill!(histowet_loc, Tf(0))
    end
    fill!(Tsnow_loc, Tm)

    if snowthickness >= eps(Tf)

        # Step 3: Restrict surface fine snow layering to Ds_surflay if this thickness is exceeded by the fresh snow addition
        if Nsnow_loc > 1
            Dtemp_surflay = Tf(0)
            k_surflay = 0
            for k in 1:(Nsnow_loc - 1)
                Dtemp_surflay = Dtemp_surflay + Ds_loc[k]
                if Dtemp_surflay > Ds_surflay
                    k_surflay = k
                    break
                end
            end
            if k_surflay > 0
                Ds_excess = Dtemp_surflay - Ds_surflay
                Ds_old = Ds_loc[k_surflay]
                Sice_old = Sice_loc[k_surflay]
                Sliq_old = Sliq_loc[k_surflay]
                U_old = U_loc[k_surflay]
                Ds_loc[k_surflay] = Ds_old - Ds_excess
                Sice_loc[k_surflay] = Sice_old * (Ds_loc[k_surflay] / Ds_old)
                Sliq_loc[k_surflay] = Sliq_old * (Ds_loc[k_surflay] / Ds_old)
                U_loc[k_surflay] = U_old * (Ds_loc[k_surflay] / Ds_old)
                # histowet_loc[k_surflay] unchanged
                Ds_loc[k_surflay + 1] = snowthickness - Ds_surflay
                # Merge the layers below the restricted surface layers
                # into one bottom layer (sums as plain loops so that
                # this is kernel-safe; same operation order as the
                # former vectorized expressions)
                wsum = zero(Tf)
                msum = zero(Tf)
                sicesum = zero(Tf)
                sliqsum = zero(Tf)
                usum = zero(Tf)
                for k in (k_surflay + 1):Nsnow_loc
                    wsum += (Sice_loc[k] + Sliq_loc[k]) * histowet_loc[k]
                    msum += Sice_loc[k] + Sliq_loc[k]
                    sicesum += Sice_loc[k]
                    sliqsum += Sliq_loc[k]
                    usum += U_loc[k]
                end
                histowet_loc[k_surflay + 1] = (
                    wsum +
                        (Sice_old + Sliq_old) * (1 - Ds_loc[k_surflay] / Ds_old) * histowet_loc[k_surflay]
                ) /
                    (
                    msum +
                        (Sice_old + Sliq_old) * (1 - Ds_loc[k_surflay] / Ds_old)
                )
                Sice_loc[k_surflay + 1] = sicesum + Sice_old * (1 - Ds_loc[k_surflay] / Ds_old)
                Sliq_loc[k_surflay + 1] = sliqsum + Sliq_old * (1 - Ds_loc[k_surflay] / Ds_old)
                U_loc[k_surflay + 1] = usum + U_old * (1 - Ds_loc[k_surflay] / Ds_old)
                if Nsnow_loc > k_surflay + 1
                    for k in (k_surflay + 2):Nsnow_loc
                        Ds_loc[k] = 0
                        Sice_loc[k] = 0
                        Sliq_loc[k] = 0
                        U_loc[k] = 0
                        histowet_loc[k] = 0
                    end
                    Nsnow_loc = k_surflay + 1
                end
            end
        end

        # Step 4: If one layer is too thin, merge it with the neighbouring layer of closest density
        while Nsnow_loc > 1
            # Find the thinnest layer
            kmin = first_argmin(Ds_loc, Nsnow_loc)
            if !(Ds_loc[kmin] < Ds_min)
                break
            end
            if kmin == 1
                # The thinnest layer is the top one
                # Merge top two layers
                Ds_loc[1] = Ds_loc[1] + Ds_loc[2]
                histowet_loc[1] = ((Sice_loc[1] + Sliq_loc[1]) * histowet_loc[1] + (Sice_loc[2] + Sliq_loc[2]) * histowet_loc[2]) /
                    (Sice_loc[1] + Sliq_loc[1] + Sice_loc[2] + Sliq_loc[2])
                Sice_loc[1] = Sice_loc[1] + Sice_loc[2]
                Sliq_loc[1] = Sliq_loc[1] + Sliq_loc[2]
                U_loc[1] = U_loc[1] + U_loc[2]
                if Nsnow_loc > 2
                    for k in 2:(Nsnow_loc - 1)
                        Ds_loc[k] = Ds_loc[k + 1]
                        Sice_loc[k] = Sice_loc[k + 1]
                        Sliq_loc[k] = Sliq_loc[k + 1]
                        U_loc[k] = U_loc[k + 1]
                        histowet_loc[k] = histowet_loc[k + 1]
                    end
                end
                Ds_loc[Nsnow_loc] = 0
                Sice_loc[Nsnow_loc] = 0
                Sliq_loc[Nsnow_loc] = 0
                U_loc[Nsnow_loc] = 0
                histowet_loc[Nsnow_loc] = 0
                Nsnow_loc = Nsnow_loc - 1
            elseif kmin == Nsnow_loc
                # The thinnest layer is the bottom one
                # Merge bottom two layers
                Ds_loc[Nsnow_loc - 1] = Ds_loc[Nsnow_loc - 1] + Ds_loc[Nsnow_loc]
                Ds_loc[Nsnow_loc] = 0
                histowet_loc[Nsnow_loc - 1] = (
                    (Sice_loc[Nsnow_loc - 1] + Sliq_loc[Nsnow_loc - 1]) * histowet_loc[Nsnow_loc - 1] +
                        (Sice_loc[Nsnow_loc] + Sliq_loc[Nsnow_loc]) * histowet_loc[Nsnow_loc]
                ) /
                    (
                    Sice_loc[Nsnow_loc - 1] + Sliq_loc[Nsnow_loc - 1] +
                        Sice_loc[Nsnow_loc] + Sliq_loc[Nsnow_loc]
                )
                histowet_loc[Nsnow_loc] = 0
                Sice_loc[Nsnow_loc - 1] = Sice_loc[Nsnow_loc - 1] + Sice_loc[Nsnow_loc]
                Sice_loc[Nsnow_loc] = 0
                Sliq_loc[Nsnow_loc - 1] = Sliq_loc[Nsnow_loc - 1] + Sliq_loc[Nsnow_loc]
                Sliq_loc[Nsnow_loc] = 0
                U_loc[Nsnow_loc - 1] = U_loc[Nsnow_loc - 1] + U_loc[Nsnow_loc]
                U_loc[Nsnow_loc] = 0
                Nsnow_loc = Nsnow_loc - 1
            else
                # The thinnest layer is in the middle
                # Merge with neighbour of closest density (either above or below)
                kup = kmin - 1
                kdown = kmin + 1
                rho_kup = (Sice_loc[kup] + Sliq_loc[kup]) / Ds_loc[kup] / fsnow[i, j]
                rho_kdown = (Sice_loc[kdown] + Sliq_loc[kdown]) / Ds_loc[kdown] / fsnow[i, j]
                rho_kmin = (Sice_loc[kmin] + Sliq_loc[kmin]) / Ds_loc[kmin] / fsnow[i, j]
                if abs(rho_kmin - rho_kup) < abs(rho_kmin - rho_kdown)
                    # Layer with closest density is up
                    # Merge with upper neighbour
                    Ds_loc[kup] = Ds_loc[kup] + Ds_loc[kmin]
                    histowet_loc[kup] = (
                        (Sice_loc[kup] + Sliq_loc[kup]) * histowet_loc[kup] +
                            (Sice_loc[kmin] + Sliq_loc[kmin]) * histowet_loc[kmin]
                    ) /
                        (
                        Sice_loc[kup] + Sliq_loc[kup] +
                            Sice_loc[kmin] + Sliq_loc[kmin]
                    )
                    Sice_loc[kup] = Sice_loc[kup] + Sice_loc[kmin]
                    Sliq_loc[kup] = Sliq_loc[kup] + Sliq_loc[kmin]
                    U_loc[kup] = U_loc[kup] + U_loc[kmin]
                    for k in kmin:(Nsnow_loc - 1)
                        Ds_loc[k] = Ds_loc[k + 1]
                        Sice_loc[k] = Sice_loc[k + 1]
                        Sliq_loc[k] = Sliq_loc[k + 1]
                        U_loc[k] = U_loc[k + 1]
                        histowet_loc[k] = histowet_loc[k + 1]
                    end
                    Ds_loc[Nsnow_loc] = 0
                    Sice_loc[Nsnow_loc] = 0
                    Sliq_loc[Nsnow_loc] = 0
                    U_loc[Nsnow_loc] = 0
                    histowet_loc[Nsnow_loc] = 0
                    Nsnow_loc = Nsnow_loc - 1
                else
                    # Layer with closest density is down
                    # Merge with lower neighbour
                    Ds_loc[kmin] = Ds_loc[kmin] + Ds_loc[kdown]
                    histowet_loc[kmin] = (
                        (Sice_loc[kmin] + Sliq_loc[kmin]) * histowet_loc[kmin] +
                            (Sice_loc[kdown] + Sliq_loc[kdown]) * histowet_loc[kdown]
                    ) /
                        (
                        Sice_loc[kmin] + Sliq_loc[kmin] +
                            Sice_loc[kdown] + Sliq_loc[kdown]
                    )
                    Sice_loc[kmin] = Sice_loc[kmin] + Sice_loc[kdown]
                    Sliq_loc[kmin] = Sliq_loc[kmin] + Sliq_loc[kdown]
                    U_loc[kmin] = U_loc[kmin] + U_loc[kdown]
                    if kdown < Nsnow_loc
                        for k in kdown:(Nsnow_loc - 1)
                            Ds_loc[k] = Ds_loc[k + 1]
                            Sice_loc[k] = Sice_loc[k + 1]
                            Sliq_loc[k] = Sliq_loc[k + 1]
                            U_loc[k] = U_loc[k + 1]
                            histowet_loc[k] = histowet_loc[k + 1]
                        end
                    end
                    Ds_loc[Nsnow_loc] = 0
                    Sice_loc[Nsnow_loc] = 0
                    Sliq_loc[Nsnow_loc] = 0
                    U_loc[Nsnow_loc] = 0
                    histowet_loc[Nsnow_loc] = 0
                    Nsnow_loc = Nsnow_loc - 1
                end
            end
        end

        # Step 5: If too many layers, merge the two neighbours with closest density
        while Nsnow_loc > Nsmax
            # Compute the density of each layer
            for k in 1:Nsnow_loc
                rho[k] = (Sice_loc[k] + Sliq_loc[k]) / Ds_loc[k] / fsnow[i, j]
            end
            for k in 1:(Nsnow_loc - 1)
                # Compute the density difference between each layer and its bottom neighbour
                diff_rho[k] = abs(rho[k] - rho[k + 1])
            end
            # Find neighbours with smallest density difference
            kmerge = first_argmin(diff_rho, Nsnow_loc - 1)
            Ds_loc[kmerge] = Ds_loc[kmerge] + Ds_loc[kmerge + 1]
            histowet_loc[kmerge] = (
                (Sice_loc[kmerge] + Sliq_loc[kmerge]) * histowet_loc[kmerge] +
                    (Sice_loc[kmerge + 1] + Sliq_loc[kmerge + 1]) * histowet_loc[kmerge + 1]
            ) /
                (
                Sice_loc[kmerge] + Sliq_loc[kmerge] +
                    Sice_loc[kmerge + 1] + Sliq_loc[kmerge + 1]
            )
            Sice_loc[kmerge] = Sice_loc[kmerge] + Sice_loc[kmerge + 1]
            Sliq_loc[kmerge] = Sliq_loc[kmerge] + Sliq_loc[kmerge + 1]
            U_loc[kmerge] = U_loc[kmerge] + U_loc[kmerge + 1]
            if kmerge + 1 < Nsnow_loc
                for k in (kmerge + 1):(Nsnow_loc - 1)
                    Ds_loc[k] = Ds_loc[k + 1]
                    Sice_loc[k] = Sice_loc[k + 1]
                    Sliq_loc[k] = Sliq_loc[k + 1]
                    U_loc[k] = U_loc[k + 1]
                    histowet_loc[k] = histowet_loc[k + 1]
                end
            end
            Ds_loc[Nsnow_loc] = 0
            Sice_loc[Nsnow_loc] = 0
            Sliq_loc[Nsnow_loc] = 0
            U_loc[Nsnow_loc] = 0
            histowet_loc[Nsnow_loc] = 0
            Nsnow_loc = Nsnow_loc - 1
        end

        # Step 6: If more layers could be used, split the thickest ones
        if Nsnow_loc > 0
            while Nsnow_loc < Nsmax
                if Nsnow_loc == 1
                    # Only one layer
                    kmax = 1
                    # If the thickest layer is too thin to split, stop
                    if Ds_loc[kmax] < Tf(2.0) * Ds_min
                        break
                    end
                    if Ds_loc[kmax] / Tf(2.0) <= Ds_surflay
                        # We can split the bottom layer in two
                        # and the new surface layer will be thinner than Ds_surflay.
                        Nsnow_loc = Nsnow_loc + 1
                        Ds_loc[Nsnow_loc] = Ds_loc[kmax] / Tf(2.0)
                        Ds_loc[kmax] = Ds_loc[kmax] / Tf(2.0)
                        Sice_loc[Nsnow_loc] = Sice_loc[kmax] / Tf(2.0)
                        Sice_loc[kmax] = Sice_loc[kmax] / Tf(2.0)
                        Sliq_loc[Nsnow_loc] = Sliq_loc[kmax] / Tf(2.0)
                        Sliq_loc[kmax] = Sliq_loc[kmax] / Tf(2.0)
                        histowet_loc[Nsnow_loc] = histowet_loc[kmax]
                        histowet_loc[kmax] = histowet_loc[kmax]
                        U_loc[Nsnow_loc] = U_loc[kmax] / Tf(2.0)
                        U_loc[kmax] = U_loc[kmax] / Tf(2.0)
                    else
                        # The bottom layer is more than twice Ds_surflay thick.
                        # We can only remove Ds_surflay from it.
                        Nsnow_loc = Nsnow_loc + 1
                        wt = Ds_surflay / Ds_loc[kmax] # Ratio of layer taken away
                        Ds_loc[Nsnow_loc] = (Tf(1.0) - wt) * Ds_loc[kmax]
                        Ds_loc[kmax] = wt * Ds_loc[kmax]
                        Sice_loc[Nsnow_loc] = (Tf(1.0) - wt) * Sice_loc[kmax]
                        Sice_loc[kmax] = wt * Sice_loc[kmax]
                        Sliq_loc[Nsnow_loc] = (Tf(1.0) - wt) * Sliq_loc[kmax]
                        Sliq_loc[kmax] = wt * Sliq_loc[kmax]
                        histowet_loc[Nsnow_loc] = histowet_loc[kmax]
                        histowet_loc[kmax] = histowet_loc[kmax]
                        U_loc[Nsnow_loc] = (Tf(1.0) - wt) * U_loc[kmax]
                        U_loc[kmax] = wt * U_loc[kmax]
                    end
                else
                    # More than one layer
                    # Find the thickest layer
                    kmax = first_argmax(Ds_loc, Nsnow_loc)
                    # If the thickest layer is too thin to split, stop
                    if Ds_loc[kmax] < Tf(2.0) * Ds_min
                        break
                    end
                    Dtemp_surflay = Tf(0)
                    for k in 1:(Nsnow_loc - 1)
                        Dtemp_surflay += Ds_loc[k]
                    end
                    if (kmax == Nsnow_loc) && (Ds_surflay - Dtemp_surflay > Ds_min)
                        # The thickest layer is the bottom one
                        # AND we can add at least Ds_min to the surface layers before reaching the max Ds_surflay
                        # There are surface layers on top of the thickest (Nsnow_loc >1)
                        # Calculate the thickness of the surface layers
                        if Ds_loc[kmax] - (Ds_surflay - Dtemp_surflay) > Ds_min
                            # We can remove (Ds_surflay - Dtemp_surflay) from the bottom layer to a new layer on top of it.
                            # The bottom layer will still be thicker than Ds_min.
                            Nsnow_loc = Nsnow_loc + 1
                            wt = (Ds_surflay - Dtemp_surflay) / Ds_loc[kmax] # Ratio of layer taken away
                            Ds_loc[Nsnow_loc] = (Tf(1.0) - wt) * Ds_loc[kmax]
                            Ds_loc[kmax] = wt * Ds_loc[kmax]
                            Sice_loc[Nsnow_loc] = (Tf(1.0) - wt) * Sice_loc[kmax]
                            Sice_loc[kmax] = wt * Sice_loc[kmax]
                            Sliq_loc[Nsnow_loc] = (Tf(1.0) - wt) * Sliq_loc[kmax]
                            Sliq_loc[kmax] = wt * Sliq_loc[kmax]
                            histowet_loc[Nsnow_loc] = histowet_loc[kmax]
                            histowet_loc[kmax] = histowet_loc[kmax]
                            U_loc[Nsnow_loc] = (Tf(1.0) - wt) * U_loc[kmax]
                            U_loc[kmax] = wt * U_loc[kmax]
                        else
                            # If we fill the surface layers to the max Ds_surflay, the bottom layer will be thinner than Ds_min.
                            # Ds_loc[kmax] > 2 * Ds_min anyway
                            # We then remove Ds_min from the bottom layer to a new layer on top of it.
                            Nsnow_loc = Nsnow_loc + 1
                            wt = Ds_min / Ds_loc[kmax] # Ratio of layer taken away
                            Ds_loc[Nsnow_loc] = (Tf(1.0) - wt) * Ds_loc[kmax]
                            Ds_loc[kmax] = wt * Ds_loc[kmax]
                            Sice_loc[Nsnow_loc] = (Tf(1.0) - wt) * Sice_loc[kmax]
                            Sice_loc[kmax] = wt * Sice_loc[kmax]
                            Sliq_loc[Nsnow_loc] = (Tf(1.0) - wt) * Sliq_loc[kmax]
                            Sliq_loc[kmax] = wt * Sliq_loc[kmax]
                            histowet_loc[Nsnow_loc] = histowet_loc[kmax]
                            histowet_loc[kmax] = histowet_loc[kmax]
                            U_loc[Nsnow_loc] = (Tf(1.0) - wt) * U_loc[kmax]
                            U_loc[kmax] = wt * U_loc[kmax]
                        end
                    else
                        # The thickest layer is not the last one.
                        # OR it is, but surface layers are already too full.
                        # In case of the second condition, we need to recalculate kmax excluding bottom layer.
                        kmax = first_argmax(Ds_loc, Nsnow_loc - 1)
                        # If the thickest layer is too thin to split, stop
                        if Ds_loc[kmax] < Tf(2.0) * Ds_min
                            break
                        end
                        # So we split an internal layer and shift layers down to make space
                        Nsnow_loc = Nsnow_loc + 1
                        for k in Nsnow_loc:-1:(kmax + 2)
                            Ds_loc[k] = Ds_loc[k - 1]
                            Sice_loc[k] = Sice_loc[k - 1]
                            Sliq_loc[k] = Sliq_loc[k - 1]
                            histowet_loc[k] = histowet_loc[k - 1]
                            U_loc[k] = U_loc[k - 1]
                        end
                        Ds_loc[kmax + 1] = Ds_loc[kmax] / Tf(2.0)
                        Ds_loc[kmax] = Ds_loc[kmax] / Tf(2.0)
                        Sice_loc[kmax + 1] = Sice_loc[kmax] / Tf(2.0)
                        Sice_loc[kmax] = Sice_loc[kmax] / Tf(2.0)
                        Sliq_loc[kmax + 1] = Sliq_loc[kmax] / Tf(2.0)
                        Sliq_loc[kmax] = Sliq_loc[kmax] / Tf(2.0)
                        histowet_loc[kmax + 1] = histowet_loc[kmax]
                        histowet_loc[kmax] = histowet_loc[kmax]
                        U_loc[kmax + 1] = U_loc[kmax] / Tf(2.0)
                        U_loc[kmax] = U_loc[kmax] / Tf(2.0)
                    end
                end
            end
        end

        # Step 7: Diagnose snow layer temperatures
        for k in 1:Nsnow_loc
            if fsnow[i, j] > eps(Tf)
                csnow_loc[k] = (Sice_loc[k] * hcap_ice + Sliq_loc[k] * hcap_wat) / fsnow[i, j]
                Tsnow_loc[k] = Tm + U_loc[k] / csnow_loc[k]
            else
                Tsnow_loc[k] = Tm
            end
            # Bring back histowet to the [0,1] range if it is out by epsilon
            histowet_loc[k] = min(max(histowet_loc[k], Tf(0.0)), Tf(1.0))
        end

    end

    # Step 8: Copy local variables to state variables
    for k in 1:Nsmax
        Sice[k, i, j] = Sice_loc[k]
        Sliq[k, i, j] = Sliq_loc[k]
        Ds[k, i, j] = Ds_loc[k]
        histowet[k, i, j] = histowet_loc[k]
        Tsnow[k, i, j] = Tsnow_loc[k]
    end
    Nsnow[i, j] = Nsnow_loc

    return nothing
end

# High-level layering step: accumulate new snow, update the snow-cover fraction and relayer.
# Called per-cell by snow_kernel! (processes/snow.jl) and by the transport relayer! pass.
"""
    snow_layering!(layering, snowfraction, i, j, state, diag, surface, grid, params, meteo, update_hist, ::Val{Nsmax})

Accumulation of new snow, snow cover fraction update and relayering at cell `(i, j)`,
after the melt, sublimation and compaction of the same step. The snow cover fraction
update is [`snowcoverfraction_point!`](@ref); `update_hist` refreshes the 14-day
history state and is resolved by the caller, since `Dates` cannot run in a kernel.
"""
# @propagate_inbounds: this is the link that carries the kernel's inbounds context
# down to relayer_snow!, whose MVector scratch would otherwise go to the heap
Base.@propagate_inbounds function snow_layering!(
        layering::AbstractLayering{Tf}, snow_fraction::AbstractSnowFraction{Tf},
        i, j, state, diag, surface, grid, params, meteo, update_hist::Bool, ::Val{Nsmax},
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

    if layering isa OriginalLayering
        Sice[1, i, j] = Sice[1, i, j] + Sice0[i, j]
    end
    snowdepth = column_sum(Ds, i, j) * fsnow[i, j] + snowdepth0[i, j]

    # Store previous snow cover fraction
    fold = fsnow[i, j]
    # Updated Fractional Snow-Covered Area
    SWEtmp = column_sum(Sice, i, j) + column_sum(Sliq, i, j)
    if layering isa DensityLayering
        SWEtmp = SWEtmp + Sice0[i, j]
    end

    snowcoverfraction_point!(
        snow_fraction, state, surface, snowdepth, SWEtmp, hfsn, i, j, update_hist
    )

    # Rescale Ds with new snow cover fraction
    if fsnow[i, j] > eps(Tf)
        Ds0[i, j] = snowdepth0[i, j] / fsnow[i, j]
        # Update surface layer thickness based on new fsnow
        if layering isa OriginalLayering
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

    relayer_snow!(layering, i, j, state, diag, grid, params, snowdepth, Tsnow0, Val(Nsmax))

    return nothing
end
