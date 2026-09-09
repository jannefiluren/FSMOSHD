struct SoilSubstrate{Tf} <: AbstractSubstrate{Tf} end
struct IceSubstrate{Tf} <: AbstractSubstrate{Tf} end

SoilSubstrate{Tf}(grid::Grid; kwargs...) where {Tf} = SoilSubstrate{Tf}()
IceSubstrate{Tf}(grid::Grid; kwargs...) where {Tf} = IceSubstrate{Tf}()

"""
    soil_properties!(substrate, i, j, state, diag, surface, grid, params)

Fill the soil heat capacity `diag.csoil[1:Nsoil, i, j]`, thermal conductivity
`diag.ksoil[1:Nsoil, i, j]` and surface moisture conductance `diag.gs1[i, j]` for cell
`(i, j)`, implemented for every `AbstractSubstrate`. Called from the `thermal!` kernel.
"""
function soil_properties! end

@inline function soil_properties!(::IceSubstrate{Tf}, i, j, state, diag, surface, grid, params) where {Tf}
    @unpack_constants(Tf)
    (; Dzsoil, Nsoil) = grid
    (; gsat) = params
    (; csoil, ksoil, gs1) = diag

    for k in 1:Nsoil
        # hcap_ice is a specific heat capacity and needs converting to a volumetric one
        csoil[k, i, j] = hcap_ice * rho_ice * Dzsoil[k]
        ksoil[k, i, j] = hcon_ice
        # An ice surface behaves like saturated soil for surface moisture conductance
        gs1[i, j] = gsat
    end
    return nothing
end

@inline function soil_properties!(::SoilSubstrate{Tf}, i, j, state, diag, surface, grid, params) where {Tf}
    @unpack_constants(Tf)
    (; Dzsoil, Nsoil) = grid
    (; gsat) = params
    (; b, hcap_soil, hcon_soil, sathh, Vcrit, Vsat) = surface
    (; theta, Tsoil) = state
    (; csoil, ksoil, gs1) = diag

    dPsidT = -rho_ice * Lf / (rho_wat * grav * Tm)

    for k in 1:Nsoil
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
    return nothing
end

"""
    soil_temperature!(i, j, state, diag, grid, params, ::Val{Nsoil})

Advance the soil column temperature `state.Tsoil[1:Nsoil, i, j]` at cell `(i, j)` by one
step, solving the tridiagonal heat conduction system driven by `diag.Gsoil`.
"""
@inline function soil_temperature!(i, j, state, diag, grid, params, ::Val{Nsoil}) where {Nsoil}

    (; dt) = params
    (; Dzsoil) = grid
    (; Tsoil) = state
    (; csoil, ksoil, Gsoil) = diag
    Tf = eltype(Tsoil)

    # Kernel-local scratch
    a = zero(MVector{Nsoil, Tf})
    b = zero(MVector{Nsoil, Tf})
    c = zero(MVector{Nsoil, Tf})
    dTs = zero(MVector{Nsoil, Tf})
    Gs = zero(MVector{Nsoil, Tf})
    rhs = zero(MVector{Nsoil, Tf})
    gamma = zero(MVector{Nsoil, Tf})

    for k in 1:(Nsoil - 1)
        Gs[k] = Tf(2) / (Dzsoil[k] / ksoil[k, i, j] + Dzsoil[k + 1] / ksoil[k + 1, i, j])
    end
    a[1] = Tf(0)
    b[1] = csoil[1, i, j] + Gs[1] * dt
    c[1] = -Gs[1] * dt
    rhs[1] = (Gsoil[i, j] - Gs[1] * (Tsoil[1, i, j] - Tsoil[2, i, j])) * dt
    for k in 2:(Nsoil - 1)
        a[k] = c[k - 1]
        b[k] = csoil[k, i, j] + (Gs[k - 1] + Gs[k]) * dt
        c[k] = -Gs[k] * dt
        rhs[k] = Gs[k - 1] * (Tsoil[k - 1, i, j] - Tsoil[k, i, j]) * dt + Gs[k] * (Tsoil[k + 1, i, j] - Tsoil[k, i, j]) * dt
    end
    k = Nsoil
    Gs[k] = ksoil[k, i, j] / Dzsoil[k]
    a[k] = c[k - 1]
    b[k] = csoil[k, i, j] + (Gs[k - 1] + Gs[k]) * dt
    c[k] = Tf(0)
    rhs[k] = Gs[k - 1] * (Tsoil[k - 1, i, j] - Tsoil[k, i, j]) * dt
    tridiag!(dTs, Nsoil, gamma, Nsoil, a, b, c, rhs)
    for k in 1:Nsoil
        Tsoil[k, i, j] = Tsoil[k, i, j] + dTs[k]
    end
    return nothing
end

"""
    cap_soil_temperature!(substrate, i, j, state, grid)

Cap the substrate temperature at cell `(i, j)`, implemented for every
`AbstractSubstrate`. Glacier ice cannot exceed the melting point, so `IceSubstrate`
clamps it there and discards the excess energy; `SoilSubstrate` is a no-op.
"""
function cap_soil_temperature! end

@inline cap_soil_temperature!(::SoilSubstrate, i, j, state, grid) = nothing

@inline function cap_soil_temperature!(::IceSubstrate{Tf}, i, j, state, grid) where {Tf}
    @unpack_constants(Tf)
    (; Nsoil) = grid
    (; Tsoil) = state

    for k in 1:Nsoil
        Tsoil[k, i, j] = min(Tsoil[k, i, j], Tm)
    end
    return nothing
end

"""
    soil!(fsm)

Soil thermal processes: the temperature of the soil or glacier ice column.

# Arguments
- `fsm::FSM`: Model state structure
"""
function soil!(fsm::FSM{Tf}) where {Tf <: Real}

    (; substrate) = fsm.physics
    (; Nsoil) = fsm.grid

    backend = get_backend(fsm.state.Tsoil)
    kernel! = soil_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.grid, fsm.params,
        substrate, Val(Int(Nsoil));
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel inbounds = true function soil_kernel!(
        state, diag, surface, grid, params::Parameters{Tf},
        substrate::AbstractSubstrate{Tf}, ::Val{Nsoil},
    ) where {Tf, Nsoil}

    i, j = @index(Global, NTuple)

    (; tthresh) = params
    (; tilefrac) = surface

    if (tilefrac[i, j] >= tthresh)

        soil_temperature!(i, j, state, diag, grid, params, Val(Nsoil))
        cap_soil_temperature!(substrate, i, j, state, grid)

    end
end
