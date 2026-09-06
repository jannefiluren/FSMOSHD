"""
    layer_sum(A, i, j)

Sum of the layer dimension of a `(Nsmax, Nx, Ny)` array at pixel (i, j).
Equivalent to `sum(A[:, i, j])` but without the overhead of the generic
reduction machinery, which matters in the snow slide hot loop.
"""
@inline function layer_sum(A::Array{Tf, 3}, i::Integer, j::Integer) where {Tf}
    s = zero(Tf)
    @inbounds for k in axes(A, 1)
        s += A[k, i, j]
    end
    return s
end

"""
    swe_from_hs(fsm, hs, i, j)

Compute the SWE (kg/m^2) contained in the top `hs` meters of the snowpack at
pixel (i, j).

Pure Julia translation of the Fortran routine `SWE_FROM_HS` (deps/SWE_FROM_HS.F90).

# Arguments
- `fsm::FSM`: Model state structure
- `hs::Real`: Snow depth at the top of the snowpack (m), averaged over the grid cell
- `i::Integer`, `j::Integer`: Grid indices
"""
function swe_from_hs(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti}, hs::Tf, i::Integer, j::Integer
    ) where {Tf <: Real, Ti <: Integer}

    (; Nsnow, fsnow, Sice, Sliq, Ds) = fsm.state
    (; rhos_min, rhos_max, rho_snow) = w

    # Epsilon Ds tolerance to avoid instabilities
    eps_Ds = Tf(1.0e-4)

    rho_avg = rho_snow
    swe = zero(Tf)
    Ds_tmp = zero(Tf)

    @inbounds if fsnow[i, j] > eps(Tf)

        dDs = hs / fsnow[i, j]
        snowthickness = layer_sum(Ds, i, j)

        if dDs > eps(Tf) && dDs <= snowthickness + eps_Ds

            if dDs >= snowthickness

                # Fix computing approximations to avoid instabilities
                swe = sum(Sice[k, i, j] + Sliq[k, i, j] for k in axes(Sice, 1))

            else

                # Normal case
                k = 1
                swe_tmp = zero(Tf)

                while k <= Nsnow[i, j] && Ds_tmp < dDs
                    swe_layer = Sice[k, i, j] + Sliq[k, i, j]
                    if dDs - Ds_tmp > Ds[k, i, j]
                        Ds_tmp += Ds[k, i, j]
                        swe_tmp += swe_layer
                    else
                        weight = (dDs - Ds_tmp) / Ds[k, i, j]
                        swe_tmp += weight * swe_layer
                        Ds_tmp = dDs
                    end
                    k += 1
                end

                swe = swe_tmp

            end

            rho_avg = swe / hs

        elseif dDs < -eps(Tf)

            error("swe_from_hs: dDs < 0 at ($i, $j): dDs = $dDs")

        elseif dDs > snowthickness + eps_Ds

            error(
                "swe_from_hs: dDs > snowthickness at ($i, $j): " *
                    "dDs = $dDs, snowthickness = $snowthickness"
            )

        end

    end

    if (rho_avg < rhos_min - Tf(0.5) || rho_avg > rhos_max + Tf(0.5)) && Ds_tmp > Tf(0.001)
        error("swe_from_hs: invalid density $rho_avg at ($i, $j)")
    end

    return swe
end

"""
    snow_ablation!(fsm, dhs, dswe, i, j, Tm)

Erode a snow depth `dhs` (m) corresponding to mass `dswe` (kg/m^2) at the top
of the snowpack at pixel (i, j), reducing the number of layers if necessary.

Pure Julia translation of the Fortran routine `SNOW_ABLATION` (deps/SNOW_ABLATION.F90).
Mutates `Sice`, `Sliq`, `Ds`, `histowet`, `Tsnow` and `Nsnow` in `fsm`.

# Arguments
- `fsm::FSM`: Model state structure
- `dhs::Real`: Snow depth to erode (m), averaged over the grid cell
- `dswe::Real`: SWE to erode (kg/m^2)
- `i::Integer`, `j::Integer`: Grid indices
- `Tm::Real`: Melting temperature (K), assigned to emptied layers
"""
function snow_ablation!(
        fsm::FSM{Tf, Ti}, dhs::Tf, dswe::Tf, i::Integer, j::Integer, Tm::Tf
    ) where {Tf <: Real, Ti <: Integer}

    (; Nsmax) = fsm.grid
    (; Ds_min) = fsm.params
    (; Sice, Sliq, Ds, histowet, Nsnow, fsnow, Tsnow) = fsm.state

    if fsnow[i, j] > eps(Tf)

        k = 0
        swe_tmp = zero(Tf)
        Ds_tmp = zero(Tf)
        swe_layer = zero(Tf)
        dDs = dhs / fsnow[i, j]

        if dDs > eps(Tf) && dswe > eps(Tf)

            while k <= Nsnow[i, j] && Ds_tmp < dDs
                k += 1
                swe_layer = Sice[k, i, j] + Sliq[k, i, j]
                Ds_tmp += Ds[k, i, j]
                swe_tmp += swe_layer
            end

            Nsnow[i, j] = Nsnow[i, j] - k + 1

            if Nsnow[i, j] >= 1
                Ds[1, i, j] = Ds_tmp - dDs
                swe_layer_new = swe_tmp - dswe
                if swe_layer > eps(Tf)
                    Sice[1, i, j] = Sice[k, i, j] * (swe_layer_new / swe_layer)
                    Sliq[1, i, j] = Sliq[k, i, j] * (swe_layer_new / swe_layer)
                else
                    Sice[1, i, j] = swe_layer_new
                    Sliq[1, i, j] = zero(Tf)
                end
                histowet[1, i, j] = histowet[k, i, j]
                Tsnow[1, i, j] = Tsnow[k, i, j]
            end

            @inbounds if Nsnow[i, j] >= 2
                for l in 2:Nsnow[i, j]
                    Ds[l, i, j] = Ds[l + k - 1, i, j]
                    Sice[l, i, j] = Sice[l + k - 1, i, j]
                    Sliq[l, i, j] = Sliq[l + k - 1, i, j]
                    histowet[l, i, j] = histowet[l + k - 1, i, j]
                    Tsnow[l, i, j] = Tsnow[l + k - 1, i, j]
                end
            end

            @inbounds for l in (Nsnow[i, j] + 1):Nsmax
                Ds[l, i, j] = zero(Tf)
                Sice[l, i, j] = zero(Tf)
                Sliq[l, i, j] = zero(Tf)
                histowet[l, i, j] = zero(Tf)
                Tsnow[l, i, j] = Tm
            end

            # If the top layer gets too thin, aggregate it with the next one, if it exists
            @inbounds if Nsnow[i, j] >= 2 && Ds[1, i, j] < Ds_min
                Ds[1, i, j] = Ds[1, i, j] + Ds[2, i, j]
                Sice[1, i, j] = Sice[1, i, j] + Sice[2, i, j]
                Sliq[1, i, j] = Sliq[1, i, j] + Sliq[2, i, j]
                histowet[1, i, j] = histowet[2, i, j]
                Tsnow[1, i, j] = Tsnow[2, i, j]
                Nsnow[i, j] = Nsnow[i, j] - 1
                if Nsnow[i, j] >= 2
                    for l in 2:Nsnow[i, j]
                        Ds[l, i, j] = Ds[l + 1, i, j]
                        Sice[l, i, j] = Sice[l + 1, i, j]
                        Sliq[l, i, j] = Sliq[l + 1, i, j]
                        histowet[l, i, j] = histowet[l + 1, i, j]
                        Tsnow[l, i, j] = Tsnow[l + 1, i, j]
                    end
                end
                for l in (Nsnow[i, j] + 1):Nsmax
                    Ds[l, i, j] = zero(Tf)
                    Sice[l, i, j] = zero(Tf)
                    Sliq[l, i, j] = zero(Tf)
                    histowet[l, i, j] = zero(Tf)
                    Tsnow[l, i, j] = Tm
                end
            end

        else

            @warn "snow_ablation!: dDs = $dDs, dswe = $dswe"

        end

    end

    return nothing
end

"""
    remove_slide_snow!(fsm, snowdepth0, Sice0, snowdepth_available, i, j, Tm)

Remove the snow depth `snowdepth_available` from pixel (i, j), taking snow
first from the fresh avalanche deposit (`snowdepth0`, `Sice0`) and then from
the snowpack itself. Returns the corresponding SWE (kg/m^2) made available for
transport.

Helper for [`snowslide_julia!`](@ref), corresponding to the "move first snow
coming from fresh avalanche deposit" block of the Fortran routine `SNOWSLIDE`
(deps/SNOWSLIDE.F90).
"""
function remove_slide_snow!(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti}, snowdepth0::Matrix{Tf}, Sice0::Matrix{Tf},
        snowdepth_available::Tf, i::Integer, j::Integer, Tm::Tf
    ) where {Tf <: Real, Ti <: Integer}

    @inbounds if snowdepth0[i, j] - snowdepth_available > eps(Tf)

        wt = snowdepth_available / snowdepth0[i, j]
        snowdepth0[i, j] = snowdepth0[i, j] - snowdepth_available
        swe_available = wt * Sice0[i, j]
        Sice0[i, j] = Sice0[i, j] - swe_available

    elseif snowdepth_available - snowdepth0[i, j] > eps(Tf)

        snowdepth_available2 = snowdepth_available - snowdepth0[i, j]
        snowdepth0[i, j] = zero(Tf)
        swe_available = Sice0[i, j]
        Sice0[i, j] = zero(Tf)

        # Compute the mass of snow available for transport
        swe_available2 = swe_from_hs(fsm, w, snowdepth_available2, i, j)

        snow_ablation!(fsm, snowdepth_available2, swe_available2, i, j, Tm)

        swe_available = swe_available + swe_available2

    else # snowdepth_available == snowdepth0[i, j]

        snowdepth0[i, j] = zero(Tf)
        swe_available = Sice0[i, j]
        Sice0[i, j] = zero(Tf)

    end

    return swe_available
end

"""
    snowslide_julia!(fsm, snowdepth0, Sice0, dSWE_slide)

Lateral redistribution of snow through gravity using the SnowSlide model.

Pure Julia implementation of the Bernhardt and Schulz (2010) SnowSlide model,
translated from the Fortran routine `SNOWSLIDE` (deps/SNOWSLIDE.F90) and
producing the same results as [`snowslide!`](@ref).
Reference: Quéno et al. (2024)

Note on orientation (as in the Fortran code): y is the W->E axis, while x is
the S->N axis, i.e. South of (i,j) is (i-1,j) and West of (i,j) is (i,j-1).

# Arguments
- `fsm::FSM`: Model state structure
- `snowdepth0::Matrix`: Snow depth of deposited snow (m) - modified in-place
- `Sice0::Matrix`: Ice content of deposited snow (kg/m²) - modified in-place
- `dSWE_slide::Matrix`: SWE change due to snow slides (kg/m²) - output
"""
function snowslide_julia!(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti}, snowdepth0::Matrix{Tf},
        Sice0::Matrix{Tf}, dSWE_slide::Matrix{Tf}
    ) where {Tf <: Real, Ti <: Integer}

    (; Nx, Ny) = fsm.grid
    (; fsnow, Ds) = fsm.state
    (; dem) = fsm.landuse
    (; dSWE_tot_slide, index_sorted_dem, slope, Shd, forestfrac) = w
    (; dyn_ratio, trig_ratio, rho_deposit, slope_min, Shd_min) = w
    (; tiled_trans_run, snow_depo, Shd_corr) = w

    # Constants
    Tm = Tf(273.15)     # TODO get from constants function...

    # Neighbour offsets in the order S, N, W, E, SW, SE, NW, NE
    offsets = ((-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1))

    # Boolean to identify pixels where snow is deposited
    fill!(snow_depo, false)

    # Threshold the snow holding depth with Shd_min
    @. Shd_corr = max(Shd, Shd_min)

    # Treating grid points from the highest to the lowest
    @inbounds for n in 1:(Nx * Ny)

        i = index_sorted_dem[n, 1]
        j = index_sorted_dem[n, 2]

        # Start slide processes only if slope higher than the defined minimum.
        # If it is a pixel receiving avalanche snow, no slope threshold.
        # No sliding if more than 50% forest
        if (slope[i, j] >= slope_min || snow_depo[i, j]) && forestfrac[i, j] < Tf(0.5)

            # Update snowdepth in case snow has been transported to this pixel earlier in the loop
            snowdepth_updated = layer_sum(Ds, i, j) * fsnow[i, j] + snowdepth0[i, j]

            # Local elevation accounting for updated snowdepth
            elev = dem[i, j] + snowdepth_updated

            # If an avalanche is occurring, the snow holding depth is reduced
            # to mimic the dynamic effect
            if snow_depo[i, j]
                Shd_corr[i, j] = max(Shd_corr[i, j] * dyn_ratio, Shd_min)
            end

            # Compute the depth of snow available for avalanche transport
            snowdepth_available = max(zero(Tf), snowdepth_updated - Shd_corr[i, j])

            # Only the open part of the pixel can slide
            # In case of tiled run, the weighting is done later when combining tiles.
            if !tiled_trans_run
                snowdepth_available = snowdepth_available * (one(Tf) - forestfrac[i, j])
            end

            # In absence of current avalanche, an hysteretic ratio of the snow holding depth
            # must be overcome to trigger a new avalanche.
            snow_excess_trigger = !snow_depo[i, j] &&
                snowdepth_updated - Shd_corr[i, j] * trig_ratio > eps(Tf)
            # If an avalanche is arriving at the pixel, the available snowdepth over holding
            # threshold must be overcome to continue the avalanche.
            snow_excess_depo = snow_depo[i, j] && snowdepth_available > eps(Tf)

            # There is a snow excess at the pixel.
            if snow_excess_trigger || snow_excess_depo

                if i == 1 || i == Nx || j == 1 || j == Ny
                    # Case 1: edge pixels. Excess snow is dumped out of the domain
                    # to avoid accumulation artefacts. It is not routed to domain pixels.

                    swe_available = remove_slide_snow!(
                        fsm, w, snowdepth0, Sice0, snowdepth_available, i, j, Tm
                    )

                    dSWE_slide[i, j] = dSWE_slide[i, j] - swe_available
                    dSWE_tot_slide[i, j] = dSWE_tot_slide[i, j] - swe_available

                else
                    # Case 2: inner domain pixels. Excess snow is routed to lower pixels if there are.
                    # Mass transfers are weighted by elevation (snowdepth included) differences.
                    # The forest fraction reduces the weight of the pixel. No snow is
                    # transferred to fully forested pixels.

                    wts = map(offsets) do (di, dj)
                        @inbounds elev_nb = dem[i + di, j + dj] +
                            layer_sum(Ds, i + di, j + dj) * fsnow[i + di, j + dj] +
                            snowdepth0[i + di, j + dj]
                        @inbounds return max(zero(Tf), elev - elev_nb) *
                            (one(Tf) - forestfrac[i + di, j + dj])
                    end

                    delev_tot = sum(wts)
                    # if delev_tot == 0 then it's a "sink", no avalanche possible.

                    if delev_tot > eps(Tf)

                        swe_available = remove_slide_snow!(
                            fsm, w, snowdepth0, Sice0, snowdepth_available, i, j, Tm
                        )

                        dSWE_slide[i, j] = dSWE_slide[i, j] - swe_available
                        dSWE_tot_slide[i, j] = dSWE_tot_slide[i, j] - swe_available

                        # Transport to neighbour pixels, weighted by elevation difference
                        # Higher pixels have a weight of 0
                        for ((di, dj), wk) in zip(offsets, wts)
                            wk = wk / delev_tot
                            if wk > eps(Tf)
                                dswe = wk * swe_available
                                # The lower pixel only receives snow from the open part of the
                                # higher pixel. If not tiled run, snow transfer has already been
                                # weighted. If tiled run, the SWE transfer is concentrated on the
                                # open part of the pixel, hence division by open fraction to
                                # cancel counter-weighting later when combining tiles.
                                if tiled_trans_run
                                    dswe = dswe * (one(Tf) - forestfrac[i, j]) /
                                        (one(Tf) - forestfrac[i + di, j + dj])
                                end
                                Sice0[i + di, j + dj] = Sice0[i + di, j + dj] + dswe
                                snowdepth0[i + di, j + dj] = snowdepth0[i + di, j + dj] + dswe / rho_deposit
                                snow_depo[i + di, j + dj] = true
                                dSWE_slide[i + di, j + dj] = dSWE_slide[i + di, j + dj] + dswe
                                dSWE_tot_slide[i + di, j + dj] = dSWE_tot_slide[i + di, j + dj] + dswe
                            end
                        end

                    end

                end

            end

        end

    end

    return nothing
end
