@kwdef struct FixedConductivity{Tf} <: AbstractConductivity{Tf}
    kfix::Tf = 0.24                                          # Fixed thermal conductivity of snow (W/m/K)
end

@kwdef struct DensityConductivity{Tf} <: AbstractConductivity{Tf}
    bthr::Tf = 2                                             # Snow thermal conductivity exponent (-)
end

FixedConductivity{Tf}(Nx, Ny; kwargs...) where {Tf} = FixedConductivity{Tf}(; kwargs...)
DensityConductivity{Tf}(Nx, Ny; kwargs...) where {Tf} = DensityConductivity{Tf}(; kwargs...)

"""
    snow_conductivity!(scheme, ksnow, Ds, Sice, Sliq, fsnow, Nsnow, params, i, j)

Fill `ksnow[1:Nsnow[i,j], i, j]` for one column. Every `AbstractConductivity` implements this.
Only the active layers are written; `thermal!` guards the one read that can reach below the pack.

Both schemes hold scalars only, so they are `isbits` and cross into a kernel by value.
"""
function snow_conductivity! end

@inline function snow_conductivity!(c::FixedConductivity, ksnow, Ds, Sice, Sliq, fsnow, Nsnow, params, i, j)
    for k in 1:Nsnow[i, j]
        ksnow[k, i, j] = c.kfix
    end
    return nothing
end

@inline function snow_conductivity!(c::DensityConductivity, ksnow, Ds, Sice, Sliq, fsnow, Nsnow, params, i, j)
    (; rhof, hcon_ice, rho_ice, DENSTY) = params
    Tf = eltype(Ds)
    for k in 1:Nsnow[i, j]
        rhos = rhof
        # TODO the DENSTY test goes away with DENSTY == 0 (roadmap Stage 8)
        if ((DENSTY != 0) && (Ds[k, i, j] > eps(Tf)) && fsnow[i, j] > eps(Tf))
            rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
        end
        ksnow[k, i, j] = hcon_ice * (rhos / rho_ice)^c.bthr
    end
    return nothing
end

"""
    thermal!(fsm)

Thermal property calculations for snow and soil layers.

The per-cell physics lives in `thermal_kernel!`, a KernelAbstractions kernel
launched over the whole grid (see `ebalsrf!` for the pattern). The three
former grid loops (snow conductivity, soil properties, surface layer) are
fused into one kernel; they only communicate through values of the same grid
cell, so the fusion is exact.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
"""
function thermal!(fsm::FSM{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    @unpack Dzsnow, Dzsoil, Nsmax, Nsoil, Nx, Ny = fsm

    @unpack gsat, rhof = fsm

    @unpack b, hcap_soil, hcon_soil, sathh, Vcrit, Vsat = fsm

    @unpack Ds, Nsnow, fsnow, Sice, Sliq, theta, Tsnow, Tsoil, Tveg = fsm

    @unpack tilefrac, tthresh = fsm

    @unpack CONDCT, DENSTY, SUBSTR = fsm

    @unpack ksnow, csoil, ksoil, gs1, Ds1, Ts1, ks1, Tveg0 = fsm

    # Strings cannot cross into kernels: resolve the tile test here

    backend = get_backend(gs1)
    kernel! = thermal_kernel!(backend)
    kernel!(
        ksnow, csoil, ksoil, gs1, Ds1, Ts1, ks1, Tveg0,
        Dzsoil, b, hcap_soil, hcon_soil, sathh, Vcrit, Vsat,
        Ds, Nsnow, fsnow, Sice, Sliq, theta, Tsnow, Tsoil, Tveg,
        tilefrac,
        tthresh, gsat, rhof,
        Nsoil, CONDCT, DENSTY, SUBSTR;
        ndrange = (Int(Nx), Int(Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function thermal_kernel!(
        ksnow, csoil, ksoil, gs1, Ds1, Ts1, ks1, Tveg0,
        Dzsoil, b, hcap_soil, hcon_soil,
        sathh, Vcrit, Vsat,
        Ds, Nsnow, fsnow, Sice, Sliq,
        theta, Tsnow, Tsoil, Tveg,
        tilefrac,
        tthresh::Tf, gsat::Tf, rhof::Tf,
        Nsoil::Ti, CONDCT::AbstractConductivity{Tf}, DENSTY::Ti, SUBSTR::AbstractSubstrate{Tf},
    ) where {Tf, Ti}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        # Thermal conductivity of snow

        cond_params = (; rhof, hcon_ice, rho_ice, DENSTY)
        snow_conductivity!(CONDCT, ksnow, Ds, Sice, Sliq, fsnow, Nsnow, cond_params, i, j)

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

        # Surface layer

        # GM/LQ: the following lines define the properties of the layer that interacts with the surface in EBALSRF.
        # to maintain numerical stability, this layer is always at least as thick at the top soil layer (10cm in default),
        # and layer properties incorporate soil properties for thin snowpacks.
        # IMPORTANT consequence of this trick: when adapting the snow layering, we need to ensure that the thickness
        # of the top soil layer does not exceed the max thickness of the first snow layer, otherwise we create artefacts in
        # the surface energy balance (thermal properties of first layer affected by soil even when they shouldnt be)
        # Note that this 'trick' has not yet been tested for top layers < 10cm!
        Ds1[i, j] = max(Dzsoil[1], Ds[1, i, j])
        Ts1[i, j] = Tsoil[1, i, j] + (Tsnow[1, i, j] - Tsoil[1, i, j]) * Ds[1, i, j] / Dzsoil[1]
        # Snow thermal resistance is zero when there is no snow in the first layer.
        # Required because the conductivity schemes fill only the active layers.
        snow_R = Ds[1, i, j] > zero(Tf) ? Tf(2) * Ds[1, i, j] / ksnow[1, i, j] : zero(Tf)
        ks1[i, j] = Dzsoil[1] / (snow_R + (Dzsoil[1] - Tf(2) * Ds[1, i, j]) / ksoil[1, i, j])
        if (Ds[1, i, j] > Tf(0.5) * Dzsoil[1])
            ks1[i, j] = ksnow[1, i, j]
        end
        if (Ds[1, i, j] > Dzsoil[1])
            Ts1[i, j] = Tsnow[1, i, j]
        end
        Tveg0[i, j] = Tveg[i, j]

    end
end
