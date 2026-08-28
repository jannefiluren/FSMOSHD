"""
    thermal!(fsm)

Thermal property calculations for snow and soil layers.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
"""
function thermal!(fsm::FSM{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    @unpack_constants(Tf)

    @unpack Dzsoil, Nsmax, Nsoil, Nx, Ny = fsm

    @unpack bthr, gsat, kfix, rhof = fsm

    @unpack b, hcap_soil, hcon_soil, sathh, Vcrit, Vsat = fsm

    @unpack Ds, Nsnow, fsnow, Sice, Sliq, theta, Tsnow, Tsoil, Tveg = fsm

    @unpack tilefrac, glacierfrac, tthresh = fsm

    @unpack CONDCT, DENSTY, TILE = fsm

    @unpack ksnow, csoil, ksoil, gs1, Ds1, Ts1, ks1, Tveg0 = fsm

    # Thermal conductivity of snow
    for j in 1:Ny
        for i in 1:Nx

            if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

                # Fixed
                ksnow[:, i, j] .= kfix
                if (CONDCT == 1)
                    # Density function
                    for k in 1:Nsnow[i, j]
                        rhos = rhof
                        if ((DENSTY != 0) && (Ds[k, i, j] > eps(Tf)) && fsnow[i, j] > eps(Tf))
                            rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
                        end
                        ksnow[k, i, j] = hcon_ice * (rhos / rho_ice)^bthr
                    end
                end

            end

        end
    end


    # Heat capacity and thermal conductivity of soil
    dPsidT = -rho_ice * Lf / (rho_wat * grav * Tm)
    for j in 1:Ny
        for i in 1:Nx

            if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

                for k in 1:Nsoil

                    if (TILE == "glacier") # Glacier soil properties

                        # Note that hcap_ice is specific heat capacity and has to be converted to volumetric heat capacity
                        csoil[k, i, j] = hcap_ice * rho_ice * Dzsoil[k]
                        # Use pure ice thermal conductivity
                        ksoil[k, i, j] = hcon_ice
                        # Consider that ice surface behaves like saturated soil for surface moisture conductance
                        gs1[i, j] = gsat

                    else # Normal soil properties
                        csoil[k, i, j] = hcap_soil[i, j] * Dzsoil[k]
                        ksoil[k, i, j] = hcon_soil[i, j]
                        if (theta[k, i, j] > eps(Tf))
                            dthudT = Tf(0.0)
                            sthu = theta[k, i, j]
                            sthf = Tf(0.0)
                            Tc = Tsoil[k, i, j] - Tm
                            Tmax = Tm + (sathh[i, j] / dPsidT) * (Vsat[i, j] / theta[k, i, j])^b[i, j]
                            if (Tsoil[k, i, j] < Tmax)
                                dthudT = (-dPsidT * Vsat[i, j] / (b[i, j] * sathh[i, j])) * (dPsidT * Tc / sathh[i, j])^(Tf(-1) / b[i, j] - Tf(1))
                                sthu = Vsat[i, j] * (dPsidT * Tc / sathh[i, j])^(Tf(-1) / b[i, j])
                                sthu = min(sthu, theta[k, i, j])
                                sthf = (theta[k, i, j] - sthu) * rho_wat / rho_ice
                            end
                            Mf = rho_ice * Dzsoil[k] * sthf
                            Mu = rho_wat * Dzsoil[k] * sthu
                            csoil[k, i, j] = hcap_soil[i, j] * Dzsoil[k] + hcap_ice * Mf + hcap_wat * Mu + rho_wat * Dzsoil[k] * ((hcap_wat - hcap_ice) * Tc + Lf) * dthudT
                            Smf = rho_ice * sthf / (rho_wat * Vsat[i, j])
                            Smu = sthu / Vsat[i, j]
                            thice = Tf(0.0)
                            if (Smf > eps(Tf))
                                thice = Vsat[i, j] * Smf / (Smu + Smf)
                            end
                            thwat = Tf(0.0)
                            if (Smu > eps(Tf))
                                thwat = Vsat[i, j] * Smu / (Smu + Smf)
                            end
                            hcon_sat = hcon_soil[i, j] * (hcon_wat^thwat) * (hcon_ice^thice) / (hcon_air^Vsat[i, j])
                            ksoil[k, i, j] = (hcon_sat - hcon_soil[i, j]) * (Smf + Smu) + hcon_soil[i, j]
                            if (k == 1)
                                gs1[i, j] = gsat * max((Smu * Vsat[i, j] / Vcrit[i, j])^Tf(2), Tf(1.0))
                            end

                        end

                    end

                end

            end

        end

    end


    # Surface layer
    for j in 1:Ny
        for i in 1:Nx

            if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

                # GM/LQ: the following lines define the properties of the layer that interacts with the surface in EBALSRF.
                # to maintain numerical stability, this layer is always at least as thick at the top soil layer (10cm in default),
                # and layer properties incorporate soil properties for thin snowpacks.
                # IMPORTANT consequence of this trick: when adapting the snow layering, we need to ensure that the thickness
                # of the top soil layer does not exceed the max thickness of the first snow layer, otherwise we create artefacts in
                # the surface energy balance (thermal properties of first layer affected by soil even when they shouldnt be)
                # Note that this 'trick' has not yet been tested for top layers < 10cm!
                Ds1[i, j] = max(Dzsoil[1], Ds[1, i, j])
                Ts1[i, j] = Tsoil[1, i, j] + (Tsnow[1, i, j] - Tsoil[1, i, j]) * Ds[1, i, j] / Dzsoil[1]
                ks1[i, j] = Dzsoil[1] / (Tf(2) * Ds[1, i, j] / ksnow[1, i, j] + (Dzsoil[1] - Tf(2) * Ds[1, i, j]) / ksoil[1, i, j])
                if (Ds[1, i, j] > Tf(0.5) * Dzsoil[1])
                    ks1[i, j] = ksnow[1, i, j]
                end
                if (Ds[1, i, j] > Dzsoil[1])
                    Ts1[i, j] = Tsnow[1, i, j]
                end
                Tveg0[i, j] = Tveg[i, j]

            end

        end

    end

    return nothing
end
