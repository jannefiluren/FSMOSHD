"""
    snowcoverfraction_point!(fsnow, swehist, swemin, swemax, snowdepthhist,
                             snowdepthmin, snowdepthmax, slopemu, xi, Ld,
                             snowdepth, SWEtmp, i, j, SNFRAC, hfsn, update_hist)

Snow cover fraction calculation for one grid cell, using multiple
parameterizations selected by `SNFRAC`.

Device-safe per-cell function called from the `snow_layering!` kernel (and
from the [`snowcoverfraction!`](@ref) host wrapper). The 14-day SWE and snow
depth history buffers are function-local `MVector`s; `update_hist` selects
whether the history state is refreshed this time step (the caller resolves
the "6:00 am" test, since `Dates` cannot run inside kernels).

# Arguments
- `fsnow`, `swehist`, `swemin`, `swemax`, `snowdepthhist`, `snowdepthmin`,
  `snowdepthmax`: State arrays from `FSM` (modified in-place)
- `slopemu`, `xi`, `Ld`: Terrain property arrays from `FSM`
- `snowdepth::Real`: Current snow depth (m)
- `SWEtmp::Real`: Current snow water equivalent (kg/m²)
- `i::Integer`, `j::Integer`: Grid indices
- `SNFRAC`: Snow cover fraction configuration
- `hfsn::Real`: Snowcover fraction depth scale (m)
- `update_hist::Bool`: Refresh the 14-day history state (true at 6:00 am)
"""
@inline function snowcoverfraction_point!(
        fsnow, swehist, swemin, swemax, snowdepthhist, snowdepthmin, snowdepthmax,
        slopemu, xi, Ld,
        snowdepth::Tf, SWEtmp::Tf, i::Integer, j::Integer,
        SNFRAC::Integer, hfsn::Tf, update_hist::Bool
    ) where {Tf <: Real}

    # @inbounds so that the bounds-check error paths do not capture the
    # local MVector history buffers (which would force them onto the heap,
    # allocating once per grid cell); tests run with --check-bounds=yes,
    # which overrides this
    @inbounds if SNFRAC == 0

        # calculate topo terms needed for standard deviation of snow depth (done)
        sd_snowdepth1 = exp(Tf(-1) / (Ld[i, j] / xi[i, j])^Tf(2))
        sd_snowdepth3 = slopemu[i, j]^Tf(0.309)

        # merge current SWEtmp with SWEtmp history from past 14 days
        SWEbuffer = MVector{15, Tf}(undef)
        snowdepthbuffer = MVector{15, Tf}(undef)
        SWEbuffer[1] = SWEtmp
        snowdepthbuffer[1] = snowdepth
        @inbounds for k in 1:14
            SWEbuffer[k + 1] = swehist[k, i, j]
            snowdepthbuffer[k + 1] = snowdepthhist[k, i, j]
        end

        # calculate snowdepthmin_buffer, snowdepthmax_buffer, snowdepthmin_recent
        # find indices of global min and max in SWEbuffer
        iabsmax = first_argmax(SWEbuffer, 15)
        iabsmin = first_argmin(SWEbuffer, 15)

        # find index of recent min in SWEbuffer
        # calculate diff vector of SWEBuffer
        ifinal = 1
        for iloop in 1:14
            ifinal = iloop
            diffSWEbuffer = SWEbuffer[iloop + 1] - SWEbuffer[iloop]
            if (diffSWEbuffer > Tf(0.5))
                break
            else
                ifinal = iloop + 1
            end
        end
        irecentmin = first_argmin(SWEbuffer, ifinal)

        # use indices to determine snowdepth amounts
        snowdepthmin_buffer = snowdepthbuffer[iabsmin]
        snowdepthmax_buffer = snowdepthbuffer[iabsmax]
        snowdepthmin_recent = snowdepthbuffer[irecentmin]

        # Compute storage of new snow on old snow in snowdepthbuffer
        dsnowdepth = snowdepth - snowdepthmin_buffer
        if (dsnowdepth < eps(Tf))
            dsnowdepth = Tf(0)
        end

        # compute dswemax in SWEbuffer
        dsnowdepthmax = snowdepthmax_buffer - snowdepthmin_buffer
        if (dsnowdepthmax < eps(Tf))
            dsnowdepthmax = Tf(0)
        end

        # don't accept dsnowdepthmax to be larger then dsnowdepth, otherwise larger fnsnow values
        # todo: think about doing this for snowdepthmin and snowdepthmax as well, and swemin and swemax
        if (dsnowdepthmax < dsnowdepth)
            dsnowdepthmax = dsnowdepth
        end

        # Compute storage of recent new snow on old snow in SWEbuffer (done)
        dsnowdepth_recent = snowdepth - snowdepthmin_recent
        if (dsnowdepth_recent < eps(Tf))
            dsnowdepth_recent = Tf(0)
        end

        # state variables interpeting the whole SWEtmp history, not only the past 14 days in the buffer
        # Set swemax and swemin equal to zero if no snow, same with corresponding snow depth values
        if (SWEtmp < eps(Tf))
            swemax[i, j] = Tf(0)
            swemin[i, j] = Tf(0)
        end
        if (snowdepth < eps(Tf))
            snowdepthmax[i, j] = Tf(0)
            snowdepthmin[i, j] = Tf(0)
        end

        # Set swemax and swemin equal to SWEtmp if maximum, store also snowdepthmax and snowdepthmin of those time steps
        if (SWEtmp >= swemax[i, j])
            swemax[i, j] = SWEtmp
            swemin[i, j] = SWEtmp
        end

        # BC: same as with the dsnowdepth, it is possible that snowdepth >snowdepthmax because the position of the max is determined
        #   based on SWE values.
        if (snowdepth >= snowdepthmax[i, j])
            snowdepthmax[i, j] = snowdepth
            snowdepthmin[i, j] = snowdepth
        end

        # Set swemin equal SWEtmp if smaller than swemin, same with corresponding snow depth value
        if (SWEtmp < swemax[i, j] && SWEtmp < swemin[i, j])
            swemin[i, j] = SWEtmp
        end
        if (snowdepth < snowdepthmax[i, j] && snowdepth < snowdepthmin[i, j])
            snowdepthmin[i, j] = snowdepth
        end

        ### calculating SCF
        # Initial guess of snow covered fraction
        fsnow_season = Tf(0)

        ####### seasonal scf, inserting snow depth in formulas of Helbig et al. and Egli and Jonas
        # calculate standard deviation (done)
        sd_snowdepth2 = snowdepthmax[i, j]^Tf(0.549)
        sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3
        # set completely flat pixels to values determined by Luca (instead of 1 or 0)
        if (!(slopemu[i, j] > eps(Tf)))
            sd_snowdepth0 = snowdepthmax[i, j]^Tf(0.84)
        end
        # calculate snow covered fraction
        if (snowdepthmax[i, j] > eps(Tf))
            fsnow_season = tanh(Tf(1.3) * snowdepthmin[i, j] / sd_snowdepth0)
        end

        # calculate cv
        coeff_vari = sd_snowdepth0 / snowdepthmax[i, j]

        ## scf based on dswe of last 14 days
        # calculate standard deviation of dhs, taking Luca's formula (flat field approximation)
        fsnow_nsnow = Tf(0)

        sd_snowdepth0_dhs = dsnowdepthmax^Tf(0.84)
        # calculate snow covered fraction of nsnow
        if (dsnowdepthmax > eps(Tf))
            fsnow_nsnow = tanh(dsnowdepth^Tf(0.14) + dsnowdepth / Tf(0.13))
        end
        #######

        ####### scf based on dswe_recent since last minimum
        # calculate standard deviation of dsnowdepth_recent, taking Luca's formula (flat field approximation)
        fsnow_nsnow_recent = Tf(0)

        sd_snowdepth0_dhs_recent = dsnowdepth_recent^Tf(0.84)
        # calculate snow covered fraction of nsnow with recent dswe, converting SWEtmp into snow depth
        if (dsnowdepth_recent > eps(Tf))
            fsnow_nsnow_recent = tanh(dsnowdepth_recent^Tf(0.14) + dsnowdepth_recent / Tf(0.13))
        end

        # take maximum between the two new snow scf, similar to taking the maximum of all three regimes at the end (done)
        fsnow_nsnow = max(fsnow_nsnow, fsnow_nsnow_recent)

        # RESET PART OF THE CODE IS TEMPORARILY REMOVED - SOLUTION TO BE FOUND TO AVOID INSTABILITIES
        #    !! recalculate scf_season if new snow has melted after a snow fall to account for a higher CV
        #    ! If new snow has melted away, update parameters swemin and swemax of
        #    ! "seasonal snow" so that the scf trajectory continues along last
        #    ! scf-value given by scf_nsnow, added that it is really (SWE yesterday > SWE current + threshold of 2) melting, including the threshold for more than spurious differences
        #    if (fsnow_nsnow .NE. 0 .and. fsnow_nsnow < fsnow_season .and. SWEbuffer(2) > SWEtmp + 2) then
        #      swemin(i,j)       = SWEtmp
        #      snowdepthmin(i,j) = snowdepth
        #      rhomax = swemax(i,j)/snowdepthmax(i,j) ! rhomax should remain constant with time, i.e the modelleded density at timestep of swemax
        #      snowdepthmax(i,j) = 1.3 * snowdepthmin(i,j) / (coeff_vari * atanh(fsnow_season))
        #      swemax(i,j) = rhomax * snowdepthmax(i,j)
        #      ! re-calculate standard deviation with new snowdepthmax
        #      sd_snowdepth2 = snowdepthmax(i,j)**0.549
        #      sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3
        #      ! set completely flat pixels to values determined by Luca (instead of 1 or 0)
        #      if (.not.(slopemu(i,j) > epsilon(0.0))) then
        #        sd_snowdepth0 = snowdepthmax(i,j)**0.84
        #      end if
        #      ! calculate snow covered fraction
        #      fsnow_season = tanh(1.3 * snowdepthmin(i,j) / !sd_snowdepth0)
        #      fsnow_nsnow = 0
        #    end if

        # Use the largest of the two fsnow estimates
        fsnow[i, j] = max(fsnow_season, fsnow_nsnow)

        # BC update history of SWE and hs only if they correspond to 6:00am values
        if update_hist
            @inbounds for k in 1:14
                swehist[k, i, j] = SWEbuffer[k]
                snowdepthhist[k, i, j] = snowdepthbuffer[k]
            end
        end

        fsnow[i, j] = max(fsnow[i, j], Tf(0.01))

    elseif SNFRAC == 1
        # HelbigHS
        sd_snowdepth2 = snowdepth^Tf(0.549)
        sd_snowdepth1 = exp(Tf(-1) / (Ld[i, j] / xi[i, j])^Tf(2))
        sd_snowdepth3 = slopemu[i, j]^Tf(0.309)
        sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3

        fsnow[i, j] = tanh(Tf(1.3) * snowdepth / sd_snowdepth0)

    elseif SNFRAC == 2
        # HelbigHS0
        if snowdepth == Tf(0)
            snowdepthmax[i, j] = Tf(0.0)
        end

        if snowdepth > snowdepthmax[i, j]
            snowdepthmax[i, j] = snowdepth
        end

        sd_snowdepth2 = snowdepthmax[i, j]^Tf(0.549)
        sd_snowdepth1 = exp(Tf(-1) / (Ld[i, j] / xi[i, j])^Tf(2))
        sd_snowdepth3 = slopemu[i, j]^Tf(0.309)
        sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3

        fsnow[i, j] = tanh(Tf(1.3) * snowdepth / sd_snowdepth0)

    elseif SNFRAC == 3
        # Point model
        fsnow[i, j] = snowdepth > eps(Tf) ? Tf(1.0) : Tf(0.0)

    else
        # tanh model / original FSM
        fsnow[i, j] = tanh(snowdepth / hfsn)
    end

    # Final adjustments
    if snowdepth < eps(Tf)
        fsnow[i, j] = Tf(0.0)
    else
        fsnow[i, j] = min(fsnow[i, j], Tf(1.0))
    end

    return nothing
end

"""
    snowcoverfraction!(fsm, snowdepth, SWEtmp, t, i, j, SWEbuffer, snowdepthbuffer, diffSWEbuffer)

Snow cover fraction calculation for one grid cell (host convenience wrapper
around [`snowcoverfraction_point!`](@ref), kept for API compatibility).

The buffer arguments are accepted but ignored: the history buffers are now
function-local (they were always pure workspace).

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `snowdepth::Real`: Current snow depth (m)
- `SWEtmp::Real`: Current snow water equivalent (kg/m²)
- `t::DateTime`: Current simulation time
- `i::Int, j::Int`: Grid indices
"""
function snowcoverfraction!(fsm::FSM{Tf, Ti}, snowdepth::Tf, SWEtmp::Tf, t::DateTime, i::Int, j::Int, SWEbuffer::AbstractArray{Tf}, snowdepthbuffer::AbstractArray{Tf}, diffSWEbuffer::AbstractArray{Tf}) where {Tf <: Real, Ti <: Integer}

    @unpack SNFRAC, hfsn = fsm
    @unpack fsnow, swehist, swemin, swemax = fsm
    @unpack snowdepthhist, snowdepthmin, snowdepthmax = fsm
    @unpack slopemu, xi, Ld = fsm

    # update history of SWE and hs only if they correspond to 6:00am values
    update_hist = 4.5 < hour(t) < 5.5

    snowcoverfraction_point!(
        fsnow, swehist, swemin, swemax, snowdepthhist, snowdepthmin, snowdepthmax,
        slopemu, xi, Ld,
        snowdepth, SWEtmp, i, j, SNFRAC, hfsn, update_hist
    )

    return nothing
end
