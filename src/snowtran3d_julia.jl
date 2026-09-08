"""
    hs_from_swe(fsm, swe, i, j)

Compute the snow depth (m) corresponding to the top `swe` kg/m^2 of the
snowpack at pixel (i, j).

Pure Julia translation of the Fortran routine `HS_FROM_SWE` (deps/HS_FROM_SWE.F90).
"""
function hs_from_swe(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti}, swe::Tf, i::Integer, j::Integer
    ) where {Tf <: Real, Ti <: Integer}

    (; Nsnow, fsnow, Sice, Sliq, Ds) = fsm.state
    (; rhos_min, rhos_max, rho_snow) = w

    # Epsilon SWE tolerance to avoid instabilities
    eps_swe = Tf(1.0e-4)

    rho_avg = rho_snow
    hs = zero(Tf)
    Ds_tmp = zero(Tf)

    swe_tot = sum(Sice[k, i, j] + Sliq[k, i, j] for k in axes(Sice, 1))

    @inbounds if swe > eps(Tf) && swe <= swe_tot + eps_swe

        if swe >= swe_tot

            # Fix computing approximations to avoid instabilities
            Ds_tmp = layer_sum(Ds, i, j)

        else

            # Normal case
            k = 1
            swe_tmp = zero(Tf)

            while k <= Nsnow[i, j] && swe_tmp < swe
                swe_layer = Sice[k, i, j] + Sliq[k, i, j]
                if swe - swe_tmp > swe_layer
                    Ds_tmp += Ds[k, i, j]
                    swe_tmp += swe_layer
                else
                    weight = (swe - swe_tmp) / swe_layer
                    Ds_tmp += weight * Ds[k, i, j]
                    swe_tmp = swe
                end
                k += 1
            end

        end

        hs = Ds_tmp * fsnow[i, j]

        if hs > eps(Tf)
            rho_avg = swe / hs
        elseif hs < -eps(Tf)
            error("hs_from_swe: hs < 0 at ($i, $j): hs = $hs")
        end

    elseif swe < -eps(Tf)

        error("hs_from_swe: swe < 0 at ($i, $j): swe = $swe")

    elseif swe > swe_tot + eps_swe

        error("hs_from_swe: swe > swe_tot at ($i, $j): swe = $swe, swe_tot = $swe_tot")

    end

    if (rho_avg < rhos_min - Tf(0.5) || rho_avg > rhos_max + Tf(0.5)) && Ds_tmp > Tf(0.001)
        error("hs_from_swe: invalid density $rho_avg at ($i, $j)")
    end

    return hs
end

"""
    compute_soft_snow!(fsm, w)

Determine the soft (movable) snow thickness `w.Ds_soft` from the layer wetting
history: layers from the top are soft while they are dry (`Sliq == 0`) and have
never been wet (`histowet < 0.5`).

Translation of the Fortran subroutine `compute_soft_snow` (deps/SNOWTRAN3D.F90).
"""
function compute_soft_snow!(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti}
    ) where {Tf <: Real, Ti <: Integer}

    (; Nx, Ny) = fsm.grid
    (; Ds, Nsnow, histowet, Sliq) = fsm.state
    (; Ds_soft) = w

    @inbounds for j in 1:Ny, i in 1:Nx
        k = 1
        while k <= Nsnow[i, j] && Sliq[k, i, j] < eps(Tf) && histowet[k, i, j] < Tf(0.5)
            Ds_soft[i, j] += Ds[k, i, j]
            k += 1
        end
    end

    return nothing
end

"""
    update_soft_snow!(fsm, w)

Update `w.Ds_soft` by removing layers that are too dense to be moved by the
friction velocity `w.Utau`.

Translation of the Fortran subroutine `update_soft_snow` (deps/SNOWTRAN3D.F90).
"""
function update_soft_snow!(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti}
    ) where {Tf <: Real, Ti <: Integer}

    (; Nx, Ny) = fsm.grid
    (; Ds, Nsnow, fsnow, Sice, Sliq) = fsm.state
    (; Utau, Ds_soft) = w

    @inbounds for j in 1:Ny, i in 1:Nx

        Ds_soft_old = Ds_soft[i, j]
        Ds_soft[i, j] = zero(Tf)
        hard_layer_reached = false
        k = 1

        while k <= Nsnow[i, j] && Ds_soft_old - Ds_soft[i, j] > eps(Tf) && !hard_layer_reached

            # Calculate snow layer density
            rho_layer = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]

            # Calculate the snow threshold friction velocity
            if rho_layer <= Tf(300)
                Utau_t_layer = Tf(0.1) * exp(Tf(0.003) * rho_layer)
            else
                Utau_t_layer = Tf(0.005) * exp(Tf(0.013) * rho_layer)
            end

            if Utau[i, j] - Utau_t_layer > eps(Tf)
                Ds_soft[i, j] += Ds[k, i, j]
            else
                hard_layer_reached = true
            end

            k += 1

        end

    end

    return nothing
end

"""
    surface_snow!(fsm, w)

Compute the threshold friction velocity `w.Utau_t` from the density of the top
snow layer.

Translation of the Fortran subroutine `surface_snow` (deps/SNOWTRAN3D.F90).
"""
function surface_snow!(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti}
    ) where {Tf <: Real, Ti <: Integer}

    (; Nx, Ny) = fsm.grid
    (; Nsnow, fsnow, Sice, Sliq, Ds) = fsm.state
    (; rho_snow, Utau_t) = w

    @inbounds for j in 1:Ny, i in 1:Nx

        rho_surf_snow = rho_snow
        if Nsnow[i, j] > 0 && Ds[1, i, j] > eps(Tf)
            rho_surf_snow = (Sice[1, i, j] + Sliq[1, i, j]) / Ds[1, i, j] / fsnow[i, j]
        end

        # Calculate the snow threshold friction velocity
        if rho_surf_snow <= Tf(300)
            Utau_t[i, j] = Tf(0.1) * exp(Tf(0.003) * rho_surf_snow)
        else
            Utau_t[i, j] = Tf(0.005) * exp(Tf(0.013) * rho_surf_snow)
        end

    end

    return nothing
end

"""
    getdirection!(w)

Build the wind direction index arrays: for every row (column), find the runs of
grid cells where the u (v) wind component keeps the same sign. Column 1 of each
index array holds the number of runs; columns 2k, 2k+1 hold the start and end
indices of run k.

Translation of the Fortran subroutine `getdirection` (deps/SNOWTRAN3D.F90).
"""
function getdirection!(w::SnowTransport{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    (; Nx, Ny, index_ue, index_uw, index_vn, index_vs, uwind, vwind) = w

    # Sweep looking for WESTERLY winds, looking for positive numbers
    @inbounds for i in 1:Nx
        sign1 = uwind[i, 1] <= zero(Tf) ? -one(Tf) : one(Tf)
        npairs = 0
        if sign1 > 0
            npairs = 1
            index_uw[i, 2] = 1
        end
        for j in 2:Ny
            sign1 = uwind[i, j - 1] <= zero(Tf) ? -one(Tf) : one(Tf)
            sign2 = uwind[i, j] <= zero(Tf) ? -one(Tf) : one(Tf)
            if sign2 != sign1
                if sign2 > 0
                    npairs += 1
                    index_uw[i, npairs * 2] = j
                else
                    index_uw[i, npairs * 2 + 1] = j - 1
                end
            end
        end
        sign1 = uwind[i, Ny] <= zero(Tf) ? -one(Tf) : one(Tf)
        if sign1 > 0
            index_uw[i, npairs * 2 + 1] = Ny
        end
        index_uw[i, 1] = npairs
    end

    # Sweep looking for EASTERLY winds, looking for negative numbers
    @inbounds for i in 1:Nx
        sign1 = uwind[i, 1] <= zero(Tf) ? -one(Tf) : one(Tf)
        npairs = 0
        if sign1 < 0
            npairs = 1
            index_ue[i, 2] = 1
        end
        for j in 2:Ny
            sign1 = uwind[i, j - 1] <= zero(Tf) ? -one(Tf) : one(Tf)
            sign2 = uwind[i, j] <= zero(Tf) ? -one(Tf) : one(Tf)
            if sign2 != sign1
                if sign2 < 0
                    npairs += 1
                    index_ue[i, npairs * 2] = j
                else
                    index_ue[i, npairs * 2 + 1] = j - 1
                end
            end
        end
        sign1 = uwind[i, Ny] <= zero(Tf) ? -one(Tf) : one(Tf)
        if sign1 < 0
            index_ue[i, npairs * 2 + 1] = Ny
        end
        index_ue[i, 1] = npairs
    end

    # Sweep looking for SOUTHERLY winds, looking for positive numbers
    @inbounds for j in 1:Ny
        sign1 = vwind[1, j] <= zero(Tf) ? -one(Tf) : one(Tf)
        npairs = 0
        if sign1 > 0
            npairs = 1
            index_vs[j, 2] = 1
        end
        for i in 2:Nx
            sign1 = vwind[i - 1, j] <= zero(Tf) ? -one(Tf) : one(Tf)
            sign2 = vwind[i, j] <= zero(Tf) ? -one(Tf) : one(Tf)
            if sign2 != sign1
                if sign2 > 0
                    npairs += 1
                    index_vs[j, npairs * 2] = i
                else
                    index_vs[j, npairs * 2 + 1] = i - 1
                end
            end
        end
        sign1 = vwind[Nx, j] <= zero(Tf) ? -one(Tf) : one(Tf)
        if sign1 > 0
            index_vs[j, npairs * 2 + 1] = Nx
        end
        index_vs[j, 1] = npairs
    end

    # Sweep looking for NORTHERLY winds, looking for negative numbers
    @inbounds for j in 1:Ny
        sign1 = vwind[1, j] <= zero(Tf) ? -one(Tf) : one(Tf)
        npairs = 0
        if sign1 < 0
            npairs = 1
            index_vn[j, 2] = 1
        end
        for i in 2:Nx
            sign1 = vwind[i - 1, j] <= zero(Tf) ? -one(Tf) : one(Tf)
            sign2 = vwind[i, j] <= zero(Tf) ? -one(Tf) : one(Tf)
            if sign2 != sign1
                if sign2 < 0
                    npairs += 1
                    index_vn[j, npairs * 2] = i
                else
                    index_vn[j, npairs * 2 + 1] = i - 1
                end
            end
        end
        sign1 = vwind[Nx, j] <= zero(Tf) ? -one(Tf) : one(Tf)
        if sign1 < 0
            index_vn[j, npairs * 2 + 1] = Nx
        end
        index_vn[j, 1] = npairs
    end

    return nothing
end

"""
    solve1(zU, guess, windtmp, threshold_flag)

Newton iterations for the friction velocity under saltation, following
Liston and Sturm (1998) eq. 5 (threshold_flag = 1) or the high-wind relation
z0 = 0.00734 u* - 0.0022 (threshold_flag = 2). Returns the solution for Utau.

Translation of the Fortran subroutine `solve1` (deps/SNOWTRAN3D.F90).
"""
function solve1(zU::Tf, guess::Tf, windtmp::Tf, threshold_flag::Tf) where {Tf <: Real}

    @unpack_constants(Tf)

    # Coefficient 0.12 in Liston and Sturm (1998) eq. 5 p. 500
    C_z = Tf(0.12)

    tol = Tf(1.0e-3)
    maxiter = 20
    old = guess
    xnew = old

    if threshold_flag == Tf(1)

        for _ in 1:maxiter
            denom = log(zU) - log(C_z / (Tf(2) * grav)) - Tf(2) * log(old)
            fprime = -Tf(1) + Tf(2) / old * windtmp * vkman / (denom * denom)
            funct = -old + windtmp * vkman / denom
            xnew = old - funct / fprime
            if abs(xnew - old) < tol
                return xnew
            end
            old = xnew
        end

    elseif threshold_flag == Tf(2)

        old = Tf(0.6)
        for _ in 1:maxiter
            denom = log(zU) - log(Tf(0.00734) * old - Tf(0.0022))
            fprime = -Tf(1) + windtmp * vkman * Tf(0.00734) / (Tf(0.00734) * old - Tf(0.0022)) / (denom * denom)
            funct = -old + windtmp * vkman / denom
            xnew = old - funct / fprime
            if abs(xnew - old) < tol
                return xnew
            end
            old = xnew
        end

    end

    @warn "solve1: max iteration exceeded when solving for Utau, Utau = $old"

    return xnew
end

"""
    solve_utau!(fsm, w, met)

Solve for the friction velocity `w.Utau`, surface roughness length `w.z_0` and
saltation-layer height `w.h_star`. Returns the blowing snow flag (1 if snow is
saltating anywhere in the domain, 0 otherwise).

Translation of the Fortran subroutine `solveUtau` (deps/SNOWTRAN3D.F90).
"""
function solve_utau!(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti}, met::MET{Tf, Ti}
    ) where {Tf <: Real, Ti <: Integer}

    (; Nx, Ny) = fsm.grid
    (; zU) = fsm.params
    (; z0_snow) = fsm.surface
    Ua_eff = fsm.diag.Uaeff
    (; vegsnowd_xy, Utau, Utau_t, z_0, h_star, snowthickness, veg_z0, Ds_soft) = w
    @unpack_constants(Tf)

    # Coefficients in Liston and Sturm (1998), eq. 5 and eq. 14
    C_z = Tf(0.12)
    h_const = Tf(1.6)

    bs_flag = zero(Tf)

    # Build the Utau array
    guess = Tf(0.1)
    @inbounds for j in 1:Ny, i in 1:Nx

        # Determine whether snow is saltating (this influences how Utau
        # and z_0 are computed)
        if snowthickness[i, j] <= vegsnowd_xy[i, j]

            # Saltation will not occur
            sfrac = snowthickness[i, j] / max(vegsnowd_xy[i, j], veg_z0[i, j])
            z_0[i, j] = sfrac * z0_snow[i, j] + (Tf(1) - sfrac) * veg_z0[i, j]
            z_0_tmp = min(Tf(0.25) * zU, z_0[i, j])
            Utau[i, j] = Ua_eff[i, j] * vkman / log(zU / z_0_tmp)
            h_star[i, j] = z_0[i, j] * h_const / C_z

        elseif Ds_soft[i, j] <= eps(Tf)

            # Saltation will not occur
            z_0[i, j] = z0_snow[i, j]
            Utau[i, j] = Ua_eff[i, j] * vkman / log(zU / z_0[i, j])
            h_star[i, j] = z_0[i, j]

        else

            # Saltation may occur. Test for that possibility by assuming that
            # saltation is present, solving for Utau and z_0, and comparing
            # whether Utau exceeds Utau_t.

            # To help insure that the iteration converges, set the minimum
            # wind speed to be 1.0 m/s, and the maximum wind speed to be
            # 30 m/s at 10-m height.
            windtmp = max(Tf(1), Ua_eff[i, j])
            wind_max = Tf(30) * log(zU / z0_snow[i, j]) / log(Tf(10) / z0_snow[i, j])
            windtmp = min(windtmp, wind_max)

            # For u* over 0.6, use the relation z0 = 0.00734 u* - 0.0022,
            # instead of Equation (5) in Liston and Sturm (1998).
            threshold = Tf(0.6) / vkman * log(zU / Tf(0.0022))
            threshold_flag = windtmp <= threshold ? Tf(1) : Tf(2)

            Utautmp = solve1(zU, guess, windtmp, threshold_flag)

            if Utautmp > Utau_t[i, j]

                # We have saltation
                Utau[i, j] = Utautmp
                z_0[i, j] = C_z * Utau[i, j]^2 / (Tf(2) * grav)
                h_star[i, j] = h_const * Utau[i, j]^2 / (Tf(2) * grav)
                bs_flag = one(Tf)

            else

                # We do not have saltation, but the vegetation is covered by
                # snow. Make sure Utau does not exceed Utau_t.
                z_0[i, j] = z0_snow[i, j]
                Utau[i, j] = Ua_eff[i, j] * vkman / log(zU / z_0[i, j])
                Utau[i, j] = min(Utau[i, j], Utau_t[i, j])
                h_star[i, j] = z_0[i, j] * h_const / C_z

            end

        end

    end

    return bs_flag
end

"""
    getsublim(zRH, z, RH, Ta, Utau, z_0, Utau_t, flag)

Sublimation loss rate coefficients at height `z`, following the appendix of
Liston and Sturm (1998). Returns `(V_susp, V_salt)`; `flag = 1` computes the
suspension-layer coefficient, `flag = 0` the saltation-layer coefficient.

Translation of the Fortran subroutine `getsublim` (deps/SNOWTRAN3D.F90).
"""
function getsublim(
        zRH::Tf, z::Tf, RH::Tf, Ta::Tf, Utau::Tf, z_0::Tf, Utau_t::Tf, flag::Tf
    ) where {Tf <: Real}

    @unpack_constants(Tf)

    # Constants from the Fortran CONSTANTS module not present in get_constants
    Runi = Tf(8313)         # Universal gas constant (J/kmol/K)
    xM = Tf(18.01)          # Molecular weight of water (kg/kmol)
    visc_air = Tf(13.0e-6)  # Kinematic viscosity of air (m^2/s)

    D = Tf(2.06e-5) * (Ta / Tf(273))^Tf(1.75)
    rho_sat = Tf(0.622) / (Rair * Ta) * Tf(610.78) * exp(Tf(21.875) * (Ta - Tm) / (Ta - Tf(7.66)))

    # Assume that the rh varies according to a modification to
    # Pomeroy's humidity variation with height equation.
    rh_offset = Tf(1) - Tf(0.027) * log(zRH)
    sigma = (Tf(0.01) * RH - Tf(1)) * (rh_offset + Tf(0.027) * log(z))
    sigma = min(Tf(0), sigma)
    sigma = max(Tf(-1), sigma)

    alpha = Tf(4.08) + Tf(12.6) * z
    rbar_r = Tf(4.6e-5) * z^Tf(-0.258)
    xmbar = Tf(4) / Tf(3) * pi * rho_ice * rbar_r^3 * (Tf(1) + Tf(3) / alpha + Tf(2) / alpha^2)
    rbar = ((Tf(3) * xmbar) / (Tf(4) * pi * rho_ice))^Tf(0.33)
    u_z = Utau / vkman * log(z / z_0)
    x_r = Tf(0.005) * u_z^Tf(1.36)
    wbar = Tf(1.1e7) * rbar^Tf(1.8)

    V_susp = zero(Tf)
    V_salt = zero(Tf)

    if flag == Tf(1)

        # Compute the sublimation loss rate coefficient for the suspension layer
        V_r = wbar + Tf(3) * x_r * cos(pi / Tf(4))
        xN_r = Tf(2) * rbar * V_r / visc_air
        xNu = Tf(1.79) + Tf(0.606) * xN_r^Tf(0.5)
        xSh = xNu
        tmp1 = (Ls * xM) / (Runi * Ta) - Tf(1)
        tmp2 = hcon_air * Ta * xNu
        top = Tf(2) * pi * rbar * sigma
        bottom = Ls / tmp2 * tmp1 + Tf(1) / (D * rho_sat * xSh)
        V_susp = (top / bottom) / xmbar
        V_salt = zero(Tf)

    elseif flag == Tf(0)

        # Compute the sublimation loss rate coefficient for the saltation layer
        V_rsalt = Tf(0.68) * Utau + Tf(2.3) * Utau_t
        xN_r = Tf(2) * rbar * V_rsalt / visc_air
        xNu = Tf(1.79) + Tf(0.606) * xN_r^Tf(0.5)
        xSh = xNu
        tmp1 = (Ls * xM) / (Runi * Ta) - Tf(1)
        tmp2 = hcon_air * Ta * xNu
        top = Tf(2) * pi * rbar * sigma
        bottom = Ls / tmp2 * tmp1 + Tf(1) / (D * rho_sat * xSh)
        V_salt = (top / bottom) / xmbar
        V_susp = zero(Tf)

    end

    return V_susp, V_salt
end

"""
    suspension!(fsm, w, met)

Compute the suspension flux `w.Qsusp` (vertical quadrature of the suspended
snow concentration following Kind, 1992) and the sublimation flux `w.Qsubl`.

Translation of the Fortran subroutine `suspension` (deps/SNOWTRAN3D.F90).
"""
function suspension!(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti}, met::MET{Tf, Ti}
    ) where {Tf <: Real, Ti <: Integer}

    (; Nx, Ny) = fsm.grid
    (; zRH) = fsm.params
    (; Ta, RH) = met
    (; Utau, Utau_t, z_0, h_star, uwind, vwind) = w
    (; Qsalt, conc_salt, Qsusp, Qsusp_u, Qsusp_v, Qsubl) = w
    @unpack_constants(Tf)

    # Constants from the Fortran CONSTANTS_SNOWTRAN3D module
    Up_const = Tf(2.8)      # Constant coefficient for calculation of U_p
    dz_susp = Tf(0.2)      # dz in suspension layer (m)
    ztop_susp = Tf(2.0)     # Height of the top of suspension layer (m)
    fall_vel = Tf(0.3)      # Particle-settling velocity (m/s)
    Ur_const = Tf(0.5)      # Constant coefficient beta in L&S 1998, eq. 13

    @inbounds for j in 1:Ny, i in 1:Nx

        if Qsalt[i, j] > eps(Tf)

            U_r = Utau[i, j] / vkman * log(h_star[i, j] / z_0[i, j])
            phistar_Cr = Utau[i, j] / U_r * Ur_const
            prd = phistar_Cr * Utau[i, j] / fall_vel
            U_p = Up_const * Utau_t[i, j]

            # Compute the concentration in the saltation layer (kg/m^3)
            conc_salt[i, j] = Qsalt[i, j] / (h_star[i, j] * U_p)

            nzsteps = trunc(Int, (ztop_susp - h_star[i, j]) / dz_susp)

            Qsusp[i, j] = zero(Tf)
            Qsubl[i, j] = zero(Tf)

            for iz in 1:nzsteps
                z = h_star[i, j] + Tf(0.5) * dz_susp + Tf(iz - 1) * dz_susp

                # Compute the concentration of the suspended snow at height z
                conc = conc_salt[i, j] * (
                    (prd + Tf(1)) *
                        (z / h_star[i, j])^((-fall_vel) / (vkman * Utau[i, j])) - prd
                )
                conc = max(conc, Tf(0))

                # Only do the integration if the concentration is non-zero
                if conc > eps(Tf)

                    # Compute the sublimation due to suspension
                    V_susp, _ = getsublim(
                        zRH, z, RH[i, j], Ta[i, j], Utau[i, j],
                        z_0[i, j], Utau_t[i, j], Tf(1)
                    )

                    # Perform the quadrature (summation), without the constants
                    if z == z_0[i, j]
                        z = Tf(1.2) * z_0[i, j]
                    end
                    Qsusp[i, j] += conc * log(z / z_0[i, j]) * dz_susp
                    Qsubl[i, j] += conc * V_susp * dz_susp

                end

            end

            # Finish the quadratures; include the constants for Qsusp
            Qsusp[i, j] = Utau[i, j] / vkman * Qsusp[i, j]

            # Include the sublimation contribution due to saltation
            z = h_star[i, j] / Tf(2)
            _, V_salt = getsublim(
                zRH, z, RH[i, j], Ta[i, j], Utau[i, j],
                z_0[i, j], Utau_t[i, j], Tf(0)
            )
            Qsubl[i, j] += V_salt * conc_salt[i, j] * h_star[i, j]

        else
            conc_salt[i, j] = zero(Tf)
            Qsusp[i, j] = zero(Tf)
            Qsubl[i, j] = zero(Tf)
        end

    end

    # Separate the east-west and the north-south suspended transport
    # components; the vector sum should equal Qsusp.
    @inbounds for j in 1:Ny, i in 1:Nx
        Qsusp_u[i, j] = Qsusp[i, j] * abs(uwind[i, j]) / sqrt(uwind[i, j]^2 + vwind[i, j]^2)
        Qsusp_v[i, j] = Qsusp[i, j] * abs(vwind[i, j]) / sqrt(uwind[i, j]^2 + vwind[i, j]^2)
    end

    return nothing
end

"""
    saltation!(fsm, w, delta_WE, delta_SN)

Compute the saltation flux `w.Qsalt` and its components `w.Qsalt_u`/`w.Qsalt_v`
with four directional upwind sweeps (westerly, easterly, southerly, northerly)
following Liston and Sturm (1998) eq. 9.

Translation of the Fortran subroutine `saltation` (deps/SNOWTRAN3D.F90).
"""
function saltation!(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti}, delta_WE::Tf, delta_SN::Tf
    ) where {Tf <: Real, Ti <: Integer}

    (; Nx, Ny) = fsm.grid
    (; vegsnowd_xy, forestfrac) = w
    (; Utau, Utau_t, snowthickness, uwind, vwind, Ds_soft) = w
    (; Qsalt, Qsalt_u, Qsalt_v, Qsalt_max, Qsalt_maxu, Qsalt_maxv) = w
    (; index_ue, index_uw, index_vn, index_vs) = w
    @unpack_constants(Tf)

    # Constants from the Fortran CONSTANTS_SNOWTRAN3D and PARAM_SNOWTRAN3D modules
    rho_air = Tf(1.275)              # Density of air (kg/m^3)
    fetch = Tf(500)                  # Equilibrium fetch distance (m)
    xmu = Tf(3)                      # Scaling constant for non-equilibrium saltation transport
    flag_boundary_condition = false  # Boundary condition flag
    blowby = Tf(0.75)                # Fraction of the saltation flux transferred downwind

    # Compute the maximum possible saltation flux, assuming that
    # an abundance of snow is available at the surface.
    @inbounds for j in 1:Ny, i in 1:Nx

        # For a given wind speed, find Qsalt_max
        Qsalt_max[i, j] = Tf(0.68) * rho_air / grav *
            Utau_t[i, j] / Utau[i, j] * (Utau[i, j]^2 - Utau_t[i, j]^2)
        Qsalt_max[i, j] = max(Qsalt_max[i, j], Tf(0))

        # Now weight the max saltation flux for the u and v wind
        # components, where the vector sum should equal Qsalt_max.
        Qsalt_maxu[i, j] = Qsalt_max[i, j] * abs(uwind[i, j]) / sqrt(uwind[i, j]^2 + vwind[i, j]^2)
        Qsalt_maxv[i, j] = Qsalt_max[i, j] * abs(vwind[i, j]) / sqrt(uwind[i, j]^2 + vwind[i, j]^2)

    end

    # Define an upwind boundary condition (flag_boundary_condition = false: zero
    # incoming flux, flag_boundary_condition = true: steady-state maximum incoming flux)
    if !flag_boundary_condition
        fill!(Qsalt_u, zero(Tf))
        fill!(Qsalt_v, zero(Tf))
    else
        copyto!(Qsalt_u, Qsalt_maxu)
        copyto!(Qsalt_v, Qsalt_maxv)
    end

    # Define the scaling coefficients for Eqn. 9 in L&S 1998. Don't
    # let them be greater than 1.0 or you will make more snow than
    # there was before.
    scale_EW = min(Tf(1), xmu * delta_WE / fetch)
    scale_NS = min(Tf(1), xmu * delta_SN / fetch)

    # Consider WESTERLY winds
    @inbounds for i in 1:Nx
        for k in 1:index_uw[i, 1]
            jstart = index_uw[i, k * 2] + 1
            jend = index_uw[i, k * 2 + 1]
            for j in jstart:jend
                dUtau = Utau[i, j] - Utau[i, j - 1]
                if dUtau >= eps(Tf)
                    Qsalt_u[i, j] = Qsalt_u[i, j - 1] + scale_EW * (Qsalt_maxu[i, j] - Qsalt_u[i, j - 1])
                else
                    if Qsalt_u[i, j - 1] < Qsalt_maxu[i, j]
                        Qsalt_u[i, j] = Qsalt_u[i, j - 1]
                    else
                        Qsalt_u[i, j] = max(blowby * Qsalt_u[i, j - 1], Qsalt_maxu[i, j])
                    end
                end
                # No flux from fully-forested cells.
                # For safety, because very high snow holding depth should anyway prevent it.
                if forestfrac[i, j] > Tf(0.9)
                    Qsalt_u[i, j] = zero(Tf)
                end
                # Downwind forest cover reduces the flow.
                if j + 1 <= Ny
                    if forestfrac[i, j + 1] > Tf(0.9)
                        Qsalt_u[i, j] = zero(Tf)
                    else
                        Qsalt_u[i, j] = Qsalt_u[i, j] * (one(Tf) - forestfrac[i, j + 1])
                    end
                end
            end
        end
    end

    # Consider EASTERLY winds
    @inbounds for i in 1:Nx
        for k in 1:index_ue[i, 1]
            jend = index_ue[i, k * 2]
            jstart = index_ue[i, k * 2 + 1] - 1
            for j in jstart:-1:jend
                dUtau = Utau[i, j] - Utau[i, j + 1]
                if dUtau >= eps(Tf)
                    Qsalt_u[i, j] = Qsalt_u[i, j + 1] + scale_EW * (Qsalt_maxu[i, j] - Qsalt_u[i, j + 1])
                else
                    if Qsalt_u[i, j + 1] < Qsalt_maxu[i, j]
                        Qsalt_u[i, j] = Qsalt_u[i, j + 1]
                    else
                        Qsalt_u[i, j] = max(blowby * Qsalt_u[i, j + 1], Qsalt_maxu[i, j])
                    end
                end
                # No flux from fully-forested cells.
                # For safety, because very high snow holding depth should anyway prevent it.
                if forestfrac[i, j] > Tf(0.9)
                    Qsalt_u[i, j] = zero(Tf)
                end
                # Downwind forest cover reduces the flow.
                if j - 1 >= 1
                    if forestfrac[i, j - 1] > Tf(0.9)
                        Qsalt_u[i, j] = zero(Tf)
                    else
                        Qsalt_u[i, j] = Qsalt_u[i, j] * (one(Tf) - forestfrac[i, j - 1])
                    end
                end
            end
        end
    end

    # Consider SOUTHERLY winds
    @inbounds for j in 1:Ny
        for k in 1:index_vs[j, 1]
            istart = index_vs[j, k * 2] + 1
            iend = index_vs[j, k * 2 + 1]
            for i in istart:iend
                dUtau = Utau[i, j] - Utau[i - 1, j]
                if dUtau >= eps(Tf)
                    Qsalt_v[i, j] = Qsalt_v[i - 1, j] + scale_NS * (Qsalt_maxv[i, j] - Qsalt_v[i - 1, j])
                else
                    if Qsalt_v[i - 1, j] < Qsalt_maxv[i, j]
                        Qsalt_v[i, j] = Qsalt_v[i - 1, j]
                    else
                        Qsalt_v[i, j] = max(blowby * Qsalt_v[i - 1, j], Qsalt_maxv[i, j])
                    end
                end
                # No flux from fully-forested cells.
                # For safety, because very high snow holding depth should anyway prevent it.
                if forestfrac[i, j] > Tf(0.9)
                    Qsalt_v[i, j] = zero(Tf)
                end
                # Downwind forest cover reduces the flow.
                if i + 1 <= Nx
                    if forestfrac[i + 1, j] > Tf(0.9)
                        Qsalt_v[i, j] = zero(Tf)
                    else
                        Qsalt_v[i, j] = Qsalt_v[i, j] * (one(Tf) - forestfrac[i + 1, j])
                    end
                end
            end
        end
    end

    # Consider NORTHERLY winds
    @inbounds for j in 1:Ny
        for k in 1:index_vn[j, 1]
            iend = index_vn[j, k * 2]
            istart = index_vn[j, k * 2 + 1] - 1
            for i in istart:-1:iend
                dUtau = Utau[i, j] - Utau[i + 1, j]
                if dUtau >= eps(Tf)
                    Qsalt_v[i, j] = Qsalt_v[i + 1, j] + scale_NS * (Qsalt_maxv[i, j] - Qsalt_v[i + 1, j])
                else
                    if Qsalt_v[i + 1, j] < Qsalt_maxv[i, j]
                        Qsalt_v[i, j] = Qsalt_v[i + 1, j]
                    else
                        Qsalt_v[i, j] = max(blowby * Qsalt_v[i + 1, j], Qsalt_maxv[i, j])
                    end
                end
                # No flux from fully-forested cells.
                # For safety, because very high snow holding depth should anyway prevent it.
                if forestfrac[i, j] > Tf(0.9)
                    Qsalt_v[i, j] = zero(Tf)
                end
                # Downwind forest cover reduces the flow.
                if i - 1 >= 1
                    if forestfrac[i - 1, j] > Tf(0.9)
                        Qsalt_v[i, j] = zero(Tf)
                    else
                        Qsalt_v[i, j] = Qsalt_v[i, j] * (one(Tf) - forestfrac[i - 1, j])
                    end
                end
            end
        end
    end

    # Combine the u and v components to yield the total saltation flux
    # at each grid cell.
    @inbounds for j in 1:Ny, i in 1:Nx
        Qsalt[i, j] = Qsalt_u[i, j] + Qsalt_v[i, j]
    end

    # Adjust Qsalt to account for the availability of snow for transport
    @inbounds for j in 1:Ny, i in 1:Nx
        if snowthickness[i, j] <= vegsnowd_xy[i, j]
            Qsalt[i, j] = zero(Tf)
            Qsalt_u[i, j] = zero(Tf)
            Qsalt_v[i, j] = zero(Tf)
        end
        if Ds_soft[i, j] <= eps(Tf)
            Qsalt[i, j] = zero(Tf)
            Qsalt_u[i, j] = zero(Tf)
            Qsalt_v[i, j] = zero(Tf)
        end
    end

    return nothing
end

"""
    getnewdepth_point!(fsm, w, Qs, dSWE_c, dSWE_c_loss, dSWE_c_gain, dh_c, dh_c_loss,
                       i, j, iu, ju, delta)

Per-pixel body shared by the four directional sweeps of [`getnewdepth!`](@ref):
compute the SWE loss/gain of pixel (i, j) for one flux component `Qs` with
upwind neighbour (iu, ju), including the adjustment for the case where there is
not enough erodible snow on the ground.
"""
@inline function getnewdepth_point!(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti},
        Qs::Matrix{Tf}, dSWE_c::Matrix{Tf}, dSWE_c_loss::Matrix{Tf},
        dSWE_c_gain::Matrix{Tf}, dh_c::Matrix{Tf}, dh_c_loss::Matrix{Tf},
        i::Integer, j::Integer, iu::Integer, ju::Integer, delta::Tf
    ) where {Tf <: Real, Ti <: Integer}

    (; dt) = fsm.params
    (; fsnow, Sice, Sliq) = fsm.state
    (; vegsnowd_xy, forestfrac, tiled_trans_run, rho_snow, snowthickness, Ds_soft) = w

    @inbounds begin

        swe_loc = sum(Sice[k, i, j] + Sliq[k, i, j] for k in axes(Sice, 1))
        dSWE_c_loss[i, j] = dt * Qs[i, j] * fsnow[i, j] / delta
        dSWE_c_gain[i, j] = dt * Qs[iu, ju] * fsnow[iu, ju] / delta
        # The SWE loss due to wind only occurs on the open part of the pixel.
        # The SWE gained by neighbouring is then weighted by the open fraction of the pixel.
        # If it is a tiled run or open run, the transferred SWE is then weighted differently.
        if tiled_trans_run
            # In the open tile, SWE loss stays the same since it will be weighted later
            # when combining the tiles. SWE gain from neighbouring pixel is first weighted
            # by the open fraction of origin pixel, then divided by the open fraction of
            # the target pixel in the open tile, since it will be counter-weighted later
            # when combining the tiles.
            if forestfrac[i, j] <= Tf(0.9)
                dSWE_c_gain[i, j] = dSWE_c_gain[i, j] * (one(Tf) - forestfrac[iu, ju]) /
                    (one(Tf) - forestfrac[i, j])
                # else: the flux is already zero
            end
        else
            dSWE_c_loss[i, j] = dSWE_c_loss[i, j] * (one(Tf) - forestfrac[i, j])
            dSWE_c_gain[i, j] = dSWE_c_gain[i, j] * (one(Tf) - forestfrac[iu, ju])
        end
        # No need to adjust Qs here because if thresholded, it will be done
        # in the next if block anyway
        dSWE_c_loss[i, j] = min(dSWE_c_loss[i, j], swe_loc)
        if dSWE_c_loss[i, j] > eps(Tf)
            dh_c_loss[i, j] = hs_from_swe(fsm, w, dSWE_c_loss[i, j], i, j)
        else
            dSWE_c_loss[i, j] = zero(Tf)
            dh_c_loss[i, j] = zero(Tf)
        end
        dh_c_gain = dSWE_c_gain[i, j] / rho_snow
        dh_c[i, j] = dh_c_gain - dh_c_loss[i, j]
        dSWE_c[i, j] = dSWE_c_gain[i, j] - dSWE_c_loss[i, j]

        # Make adjustments for the case where there is no snow available
        # on the ground (or captured within the vegetation) to be eroded.
        Ds_hard = snowthickness[i, j] - Ds_soft[i, j]
        snowdmin = max(vegsnowd_xy[i, j], Ds_hard)
        if snowthickness[i, j] > snowdmin && fsnow[i, j] > eps(Tf)
            if snowthickness[i, j] - dh_c_loss[i, j] / fsnow[i, j] <= snowdmin
                dh_c_loss[i, j] = (snowthickness[i, j] - snowdmin) * fsnow[i, j]
                if dh_c_loss[i, j] > eps(Tf)
                    dSWE_c_loss[i, j] = swe_from_hs(fsm, w, dh_c_loss[i, j], i, j)
                    # Same open tile weighting as before.
                    if !tiled_trans_run
                        dSWE_c_loss[i, j] = dSWE_c_loss[i, j] * (one(Tf) - forestfrac[i, j])
                    end
                else
                    dh_c_loss[i, j] = zero(Tf)
                    dSWE_c_loss[i, j] = zero(Tf)
                end
                dh_c[i, j] = dh_c_gain - dh_c_loss[i, j]
                dSWE_c[i, j] = dSWE_c_gain[i, j] - dSWE_c_loss[i, j]
                Qs[i, j] = Qs[iu, ju] - dSWE_c[i, j] * delta / dt / fsnow[i, j]
            end
        else
            Qs[i, j] = zero(Tf)
            dh_c[i, j] = dh_c_gain
            dSWE_c_loss[i, j] = zero(Tf)
            dh_c_loss[i, j] = zero(Tf)
            dSWE_c[i, j] = dSWE_c_gain[i, j]
        end

    end

    return nothing
end

"""
    getnewdepth!(fsm, w, Qs_u, Qs_v, dSWE_s, snowdepth0, Sice0, delta_WE, delta_SN, Tm)

Convert the flux components `Qs_u`/`Qs_v` (either saltation or suspension) into
SWE and depth changes per pixel: erode snow from the snowpack (via
[`snow_ablation!`](@ref)) and deposit transported snow into `Sice0`/`snowdepth0`.

Translation of the Fortran subroutine `getnewdepth` (deps/SNOWTRAN3D.F90).
"""
function getnewdepth!(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti},
        Qs_u::Matrix{Tf}, Qs_v::Matrix{Tf}, dSWE_s::Matrix{Tf},
        snowdepth0::Matrix{Tf}, Sice0::Matrix{Tf},
        delta_WE::Tf, delta_SN::Tf, Tm::Tf
    ) where {Tf <: Real, Ti <: Integer}

    (; Nx, Ny) = fsm.grid
    (; Ds) = fsm.state
    (; rho_snow, snowthickness, Ds_soft) = w
    (; index_ue, index_uw, index_vn, index_vs) = w
    (; dh_s_u, dh_s_v, dSWE_s_u, dSWE_s_v) = w
    (; dSWE_s_u_loss, dSWE_s_v_loss, dh_s_u_loss, dh_s_v_loss) = w
    (; dSWE_s_u_gain, dSWE_s_v_gain) = w

    fill!(dh_s_u, zero(Tf))
    fill!(dh_s_v, zero(Tf))
    fill!(dSWE_s_u, zero(Tf))
    fill!(dSWE_s_v, zero(Tf))
    fill!(dSWE_s_u_gain, zero(Tf))
    fill!(dSWE_s_v_gain, zero(Tf))
    fill!(dSWE_s_u_loss, zero(Tf))
    fill!(dSWE_s_v_loss, zero(Tf))
    fill!(dh_s_u_loss, zero(Tf))
    fill!(dh_s_v_loss, zero(Tf))

    # Consider WESTERLY winds
    @inbounds for i in 1:Nx
        for k in 1:index_uw[i, 1]
            jstart = index_uw[i, k * 2] + 1
            jend = index_uw[i, k * 2 + 1]
            for j in jstart:jend
                getnewdepth_point!(
                    fsm, w, Qs_u, dSWE_s_u, dSWE_s_u_loss, dSWE_s_u_gain,
                    dh_s_u, dh_s_u_loss, i, j, i, j - 1, delta_WE
                )
            end
        end
    end

    # Consider EASTERLY winds
    @inbounds for i in 1:Nx
        for k in 1:index_ue[i, 1]
            jend = index_ue[i, k * 2]
            jstart = index_ue[i, k * 2 + 1] - 1
            for j in jstart:-1:jend
                getnewdepth_point!(
                    fsm, w, Qs_u, dSWE_s_u, dSWE_s_u_loss, dSWE_s_u_gain,
                    dh_s_u, dh_s_u_loss, i, j, i, j + 1, delta_WE
                )
            end
        end
    end

    # Consider SOUTHERLY winds
    @inbounds for j in 1:Ny
        for k in 1:index_vs[j, 1]
            istart = index_vs[j, k * 2] + 1
            iend = index_vs[j, k * 2 + 1]
            for i in istart:iend
                getnewdepth_point!(
                    fsm, w, Qs_v, dSWE_s_v, dSWE_s_v_loss, dSWE_s_v_gain,
                    dh_s_v, dh_s_v_loss, i, j, i - 1, j, delta_SN
                )
            end
        end
    end

    # Consider NORTHERLY winds
    @inbounds for j in 1:Ny
        for k in 1:index_vn[j, 1]
            iend = index_vn[j, k * 2]
            istart = index_vn[j, k * 2 + 1] - 1
            for i in istart:-1:iend
                getnewdepth_point!(
                    fsm, w, Qs_v, dSWE_s_v, dSWE_s_v_loss, dSWE_s_v_gain,
                    dh_s_v, dh_s_v_loss, i, j, i + 1, j, delta_SN
                )
            end
        end
    end

    # Update the snow depth changes due to transport from the east and west,
    # and north and south.
    eps_w = Tf(1.0e-6)
    @inbounds for j in 1:Ny, i in 1:Nx

        weight_u = abs(dh_s_u[i, j]) / (abs(dh_s_u[i, j]) + abs(dh_s_v[i, j]) + eps_w)
        weight_v = abs(dh_s_v[i, j]) / (abs(dh_s_u[i, j]) + abs(dh_s_v[i, j]) + eps_w)

        dSWE_s_u[i, j] = weight_u * dSWE_s_u[i, j]
        dSWE_s_v[i, j] = weight_v * dSWE_s_v[i, j]
        dSWE_s_u_gain[i, j] = weight_u * dSWE_s_u_gain[i, j]
        dSWE_s_v_gain[i, j] = weight_v * dSWE_s_v_gain[i, j]
        dSWE_s_u_loss[i, j] = weight_u * dSWE_s_u_loss[i, j]
        dSWE_s_v_loss[i, j] = weight_v * dSWE_s_v_loss[i, j]

        dSWE_s[i, j] = dSWE_s_u[i, j] + dSWE_s_v[i, j]
        dSWE_s_gain = dSWE_s_u_gain[i, j] + dSWE_s_v_gain[i, j]
        dSWE_s_loss = dSWE_s_u_loss[i, j] + dSWE_s_v_loss[i, j]
        dh_s_gain = dSWE_s_gain / rho_snow
        if dSWE_s_loss > eps(Tf)
            dh_s_loss = hs_from_swe(fsm, w, dSWE_s_loss, i, j)
        else
            dSWE_s_loss = zero(Tf)
            dh_s_loss = zero(Tf)
        end

        Ds_soft[i, j] = max(Ds_soft[i, j] - dh_s_loss, zero(Tf))

        if dSWE_s_loss > eps(Tf) && dh_s_loss > eps(Tf)
            snow_ablation!(fsm, dh_s_loss, dSWE_s_loss, i, j, Tm)
        end

        # Net mass gain for this grid cell at this time step
        if dSWE_s_gain > eps(Tf) && dh_s_gain > eps(Tf)
            # Add to the existing top layer
            Sice0[i, j] += dSWE_s_gain
            snowdepth0[i, j] += dh_s_gain
        end

        # Update the snow layer thicknesses
        snowthickness[i, j] = layer_sum(Ds, i, j)

    end

    return nothing
end

"""
    accum!(fsm, w, snowdepth0, Sice0, dSWE_salt, dSWE_susp, dSWE_subl,
           bs_flag, delta_WE, delta_SN, Tm)

Compute the new snow depth due to accumulation from saltation and suspension,
and the mass loss due to sublimation, then update the cumulated transport
arrays `dSWE_tot_*` in `fsm`.

Translation of the Fortran subroutine `accum` (deps/SNOWTRAN3D.F90).
"""
function accum!(
        fsm::FSM{Tf, Ti}, w::SnowTransport{Tf, Ti},
        snowdepth0::Matrix{Tf}, Sice0::Matrix{Tf},
        dSWE_salt::Matrix{Tf}, dSWE_susp::Matrix{Tf}, dSWE_subl::Matrix{Tf},
        bs_flag::Tf, delta_WE::Tf, delta_SN::Tf, Tm::Tf
    ) where {Tf <: Real, Ti <: Integer}

    (; Nx, Ny) = fsm.grid
    (; dt) = fsm.params
    (; Ds, fsnow, Sice, Sliq) = fsm.state
    (; vegsnowd_xy, forestfrac, tiled_trans_run) = w
    (; dSWE_tot_subl, dSWE_tot_salt, dSWE_tot_susp) = w
    (; snowthickness, Ds_soft, Qsalt_u, Qsalt_v, Qsusp_u, Qsusp_v, Qsubl) = w

    if bs_flag == Tf(1)

        # SALTATION
        getnewdepth!(
            fsm, w, Qsalt_u, Qsalt_v, dSWE_salt,
            snowdepth0, Sice0, delta_WE, delta_SN, Tm
        )

        # SUSPENSION
        getnewdepth!(
            fsm, w, Qsusp_u, Qsusp_v, dSWE_susp,
            snowdepth0, Sice0, delta_WE, delta_SN, Tm
        )

        # SUBLIMATION
        @inbounds for j in 1:Ny, i in 1:Nx

            # Make adjustments for the case where there is no snow available
            # on the ground (or captured within the vegetation) to be eroded.
            Ds_hard = snowthickness[i, j] - Ds_soft[i, j]
            snowdmin = max(vegsnowd_xy[i, j], Ds_hard)
            swe_loc = sum(Sice[k, i, j] + Sliq[k, i, j] for k in axes(Sice, 1))

            # Convert Qsubl to sublimated snow depth dh_subl_loss
            dSWE_subl[i, j] = Qsubl[i, j] * dt * fsnow[i, j]
            dSWE_subl_loss = -dSWE_subl[i, j]
            # The following min should actually never be necessary, because handled
            # by the snowdmin check later. Kept for clarity and safety.
            dSWE_subl_loss = min(dSWE_subl_loss, swe_loc)
            dSWE_subl[i, j] = -dSWE_subl_loss
            if dSWE_subl_loss > eps(Tf)
                dh_subl_loss = hs_from_swe(fsm, w, dSWE_subl_loss, i, j)
            else
                dSWE_subl_loss = zero(Tf)
                dh_subl_loss = zero(Tf)
                dSWE_subl[i, j] = zero(Tf)
            end

            if snowthickness[i, j] > snowdmin && fsnow[i, j] > eps(Tf)
                if snowthickness[i, j] - dh_subl_loss / fsnow[i, j] <= snowdmin
                    dh_subl_loss = (snowthickness[i, j] - snowdmin) * fsnow[i, j]
                    if dh_subl_loss > eps(Tf)
                        dSWE_subl_loss = swe_from_hs(fsm, w, dh_subl_loss, i, j)
                        dSWE_subl[i, j] = -dSWE_subl_loss
                    else
                        dh_subl_loss = zero(Tf)
                        dSWE_subl_loss = zero(Tf)
                        dSWE_subl[i, j] = zero(Tf)
                    end
                end
            else
                dSWE_subl[i, j] = zero(Tf)
                dSWE_subl_loss = zero(Tf)
                dh_subl_loss = zero(Tf)
            end

            # The SWE loss due to wind only occurs on the open part of the pixel.
            # If it is a tiled run or open run, the transferred SWE is then weighted differently.
            # If tiled run, in the open tile, SWE loss stays the same since it will be
            # weighted later when combining the tiles.
            if !tiled_trans_run
                dSWE_subl[i, j] = dSWE_subl[i, j] * (one(Tf) - forestfrac[i, j])
                dSWE_subl_loss = dSWE_subl_loss * (one(Tf) - forestfrac[i, j])
                # Now, we may need to adjust dh_subl_loss.
                if dSWE_subl_loss > eps(Tf)
                    dh_subl_loss = hs_from_swe(fsm, w, dSWE_subl_loss, i, j)
                else
                    dSWE_subl_loss = zero(Tf)
                    dSWE_subl[i, j] = zero(Tf)
                    dh_subl_loss = zero(Tf)
                end
            end

            if dSWE_subl_loss > eps(Tf) && dh_subl_loss > eps(Tf)
                snow_ablation!(fsm, dh_subl_loss, dSWE_subl_loss, i, j, Tm)
            end

            # Update the snow layer thicknesses
            snowthickness[i, j] = layer_sum(Ds, i, j)

        end

    end

    # Update cumulated sublimation and transport
    @inbounds for j in 1:Ny, i in 1:Nx
        dSWE_tot_subl[i, j] += dSWE_subl[i, j]
        dSWE_tot_salt[i, j] += dSWE_salt[i, j]
        dSWE_tot_susp[i, j] += dSWE_susp[i, j]
    end

    return nothing
end

"""
    snowtran3d_julia!(fsm, met, snowdepth0, Sice0, dSWE_salt, dSWE_susp, dSWE_subl)

Snow transport by wind using Liston's SnowTran3D model.

Pure Julia implementation of the Liston and Sturm (1998) and Liston et
al. (2007) SnowTran3D model, translated from the Fortran routine `SNOWTRAN3D`
(deps/SNOWTRAN3D.F90) and producing the same results as [`snowtran3d!`](@ref)
up to Float32 round-off (the two implementations use different math libraries
for `exp`/`log`/`^` etc., so results are not guaranteed bit-identical).
Reference: Quéno et al. (2024)

Note on orientation (as in the Fortran code): y is the W->E axis, while x is
the S->N axis, i.e. South of (i,j) is (i-1,j) and West of (i,j) is (i,j-1).

# Arguments
- `fsm::FSM`: Model state structure
- `met::MET`: Meteo variable structure
- `snowdepth0::Matrix`: Snow depth of deposited snow (m) - modified in-place
- `Sice0::Matrix`: Ice content of deposited snow (kg/m²) - modified in-place
- `dSWE_salt::Matrix`: SWE change due to saltation (kg/m²) - output
- `dSWE_susp::Matrix`: SWE change due to suspension (kg/m²) - output
- `dSWE_subl::Matrix`: SWE change due to sublimation (kg/m²) - output
"""
function snowtran3d_julia!(
        fsm::FSM{Tf, Ti}, met::MET{Tf, Ti}, w::SnowTransport{Tf, Ti}, snowdepth0::Matrix{Tf}, Sice0::Matrix{Tf},
        dSWE_salt::Matrix{Tf}, dSWE_susp::Matrix{Tf},
        dSWE_subl::Matrix{Tf}
    ) where {Tf <: Real, Ti <: Integer}

    (; Nx, Ny) = fsm.grid
    (; Ds) = fsm.state
    (; Ld) = fsm.surface
    Ua_eff = fsm.diag.Uaeff
    (; vegsnowd_xy) = w
    (; Udir) = met
    (; uwind, vwind, snowthickness, veg_z0, Ds_soft, Utau_t) = w
    (; Qsalt, Qsalt_u, Qsalt_v, conc_salt, Qsusp, Qsusp_u, Qsusp_v, Qsubl) = w
    @unpack_constants(Tf)

    # Constants from the Fortran PARAM_SNOWTRAN3D and CONSTANTS_SNOWTRAN3D modules
    flag_variable_Utau_t = true  # Flag for variable threshold friction velocity (true) or constant (false)
    Utau_t_const = Tf(0.25)      # Constant threshold friction velocity (m/s) used if flag_variable_Utau_t = false
    wind_min = Tf(4)             # Minimum wind speed to compute snow transport (m/s)

    # Initialization of delta_WE and delta_SN
    delta_WE = Ld[1, 1]
    delta_SN = Ld[1, 1]

    # Initialize Ds_soft and snowdepth0
    fill!(Ds_soft, zero(Tf))
    fill!(snowdepth0, zero(Tf))

    # Initialize all fluxes at 0
    fill!(Qsalt, zero(Tf))
    fill!(Qsalt_u, zero(Tf))
    fill!(Qsalt_v, zero(Tf))
    fill!(dSWE_salt, zero(Tf))
    fill!(conc_salt, zero(Tf))
    fill!(Qsusp, zero(Tf))
    fill!(Qsusp_u, zero(Tf))
    fill!(Qsusp_v, zero(Tf))
    fill!(dSWE_susp, zero(Tf))
    fill!(Qsubl, zero(Tf))
    fill!(dSWE_subl, zero(Tf))

    # Initialization of uwind and vwind
    @. uwind = Ua_eff * cos(-Udir * pi / Tf(180) - pi / Tf(2))
    @. vwind = Ua_eff * sin(-Udir * pi / Tf(180) - pi / Tf(2))

    # Initialization of maximum wind speed on the domain (m/s)
    # and of snowthickness
    windspd_flag = zero(Tf)
    @inbounds for j in 1:Ny, i in 1:Nx
        windspd_flag = max(windspd_flag, Ua_eff[i, j])
        snowthickness[i, j] = layer_sum(Ds, i, j)
    end

    # Define the roughness lengths for each of the vegetation types
    @. veg_z0 = Tf(0.25) * vegsnowd_xy

    # Update the thicknesses of the hard and soft layers.
    compute_soft_snow!(fsm, w)

    # Update the threshold friction velocity
    if !flag_variable_Utau_t
        fill!(Utau_t, Utau_t_const)
    else
        surface_snow!(fsm, w)
    end

    # Set the blowing snow flag to zero until it is clear that we will
    # have blowing snow.
    bs_flag = zero(Tf)

    # If the wind speed is lower than some threshold, then don't
    # need to do any of the snow transport computations.
    if windspd_flag >= wind_min

        # Get the wind direction indexing arrays for this particular
        # wind event (time step).
        getdirection!(w)

        # Solve for Utau and z_0 if snow is saltating, else solve assuming
        # z_0 is known from snow depth and/or veg type, and solve for Utau.
        bs_flag = solve_utau!(fsm, w, met)

        # Update Ds_soft by removing too dense layers compared to U_tau
        update_soft_snow!(fsm, w)

        # If the blowing snow flag indicates wind transported snow
        # somewhere within the domain (bs_flag = 1.0), run the saltation
        # and suspension models.
        if bs_flag == Tf(1)
            saltation!(fsm, w, delta_WE, delta_SN)
            suspension!(fsm, w, met)
        end

    end

    # Compute the new snow depth due to accumulation from saltation and
    # suspension, and the mass loss due to sublimation.
    accum!(
        fsm, w, snowdepth0, Sice0, dSWE_salt, dSWE_susp, dSWE_subl,
        bs_flag, delta_WE, delta_SN, Tm
    )

    return nothing
end
