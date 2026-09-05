@kwdef struct FixedConductivity{Tf} <: AbstractConductivity{Tf}
    kfix::Tf = 0.24        # Fixed thermal conductivity of snow (W/m/K)
end

@kwdef struct DensityConductivity{Tf} <: AbstractConductivity{Tf}
    bthr::Tf = 2           # Snow thermal conductivity exponent (-)
end

FixedConductivity{Tf}(Nx, Ny; kwargs...) where {Tf} = FixedConductivity{Tf}(; kwargs...)
DensityConductivity{Tf}(Nx, Ny; kwargs...) where {Tf} = DensityConductivity{Tf}(; kwargs...)

"""
    snow_conductivity!(scheme, i, j, state, diag, params)

Fill the snow thermal conductivity `ksnow[1:Nsnow, i, j]` for cell `(i, j)`,
implemented for every `AbstractConductivity`.
"""
function snow_conductivity! end

@inline function snow_conductivity!(c::FixedConductivity, i, j, state, diag, params)
    (; Nsnow) = state
    (; ksnow) = diag
    for k in 1:Nsnow[i, j]
        ksnow[k, i, j] = c.kfix
    end
    return nothing
end

@inline function snow_conductivity!(c::DensityConductivity{Tf}, i, j, state, diag, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, fsnow, Nsnow) = state
    (; ksnow) = diag
    (; rhof) = params
    for k in 1:Nsnow[i, j]
        rhos = rhof
        if ((Ds[k, i, j] > eps(Tf)) && fsnow[i, j] > eps(Tf))
            rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
        end
        ksnow[k, i, j] = hcon_ice * (rhos / rho_ice)^c.bthr
    end
    return nothing
end

"""
    thermal!(fsm)

Thermal property calculations for snow and soil layers.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
"""
function thermal!(fsm::FSM{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    (; CONDCT, SUBSTR) = fsm.physics

    backend = get_backend(fsm.diag.gs1)
    kernel! = thermal_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.landuse, fsm.grid, fsm.params,
        CONDCT, SUBSTR;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function thermal_kernel!(
        state, diag, landuse, grid, params::Parameters{Tf},
        CONDCT::AbstractConductivity{Tf}, SUBSTR::AbstractSubstrate{Tf},
    ) where {Tf}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    (; Dzsoil, Nsoil) = grid
    (; tthresh, gsat) = params
    (; b, hcap_soil, hcon_soil, sathh, Vcrit, Vsat, tilefrac) = landuse
    (; Ds, theta, Tsnow, Tsoil, Tveg) = state
    (; ksnow, csoil, ksoil, gs1, Ds1, Ts1, ks1, Tveg0) = diag

    if (tilefrac[i, j] >= tthresh)

        # Thermal conductivity of snow
        snow_conductivity!(CONDCT, i, j, state, diag, params)

        # Heat capacity and thermal conductivity of soil
        dPsidT = -rho_ice * Lf / (rho_wat * grav * Tm)

        for k in 1:Nsoil

            if SUBSTR isa IceSubstrate # Ice properties

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

        # Surface layer (always at least as thick as the top soil layer) that
        # requires Dzsnow[1] >= Dzsoil[1] since otherwise Ts1 is blended with
        # Tsoil even under a deep snowpack
        Ds1[i, j] = max(Dzsoil[1], Ds[1, i, j])
        Ts1[i, j] = Tsoil[1, i, j] + (Tsnow[1, i, j] - Tsoil[1, i, j]) * Ds[1, i, j] / Dzsoil[1]

        # Series resistance of the composite surface layer with (a) guard to avoid zero
        # division for cells that never held snow and (b) a soil resistance that turns negative
        # once snow fills over half the layer and therefore ks1 is overridden by below
        snow_R = Ds[1, i, j] > zero(Tf) ? Tf(2) * Ds[1, i, j] / ksnow[1, i, j] : zero(Tf)
        soil_R = (Dzsoil[1] - Tf(2) * Ds[1, i, j]) / ksoil[1, i, j]

        ks1[i, j] = Dzsoil[1] / (snow_R + soil_R)
        if (Ds[1, i, j] > Tf(0.5) * Dzsoil[1])
            ks1[i, j] = ksnow[1, i, j]
        end
        if (Ds[1, i, j] > Dzsoil[1])
            Ts1[i, j] = Tsnow[1, i, j]
        end
        Tveg0[i, j] = Tveg[i, j]

    end
end
