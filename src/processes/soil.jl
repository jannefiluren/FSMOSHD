# Soil process: soil / glacier-ice column temperature

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
    tridiag!(dTs, Nsoil, a, b, c, rhs)
    for k in 1:Nsoil
        Tsoil[k, i, j] = Tsoil[k, i, j] + dTs[k]
    end
    return nothing
end
