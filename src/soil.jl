struct SoilSubstrate{Tf} <: AbstractSubstrate{Tf} end
struct IceSubstrate{Tf} <: AbstractSubstrate{Tf} end

SoilSubstrate{Tf}(Nx, Ny; kwargs...) where {Tf} = SoilSubstrate{Tf}()
IceSubstrate{Tf}(Nx, Ny; kwargs...) where {Tf} = IceSubstrate{Tf}()

function soil!(fsm::FSM{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    (; SUBSTR) = fsm.physics
    (; Nsoil) = fsm.grid

    backend = get_backend(fsm.state.Tsoil)
    kernel! = soil_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.landuse, fsm.grid, fsm.params,
        SUBSTR, Val(Int(Nsoil));
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel inbounds = true function soil_kernel!(
        state, diag, landuse, grid, params::Parameters{Tf},
        SUBSTR::AbstractSubstrate{Tf}, ::Val{Nsoil},
    ) where {Tf, Nsoil}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)
    (; dt, tthresh) = params
    (; Dzsoil) = grid
    (; Tsoil) = state
    (; csoil, ksoil, Gsoil) = diag
    (; tilefrac) = landuse

    if (tilefrac[i, j] >= tthresh)

        a = zero(MVector{Nsoil, Tf})
        b = zero(MVector{Nsoil, Tf})
        c = zero(MVector{Nsoil, Tf})
        dTs = zero(MVector{Nsoil, Tf})
        Gs = zero(MVector{Nsoil, Tf})
        rhs = zero(MVector{Nsoil, Tf})
        gamma = zero(MVector{Nsoil, Tf})

        # Soil temperature update
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

        # In case of ice substracte cap temperatures to melting point (not energy conserving)
        if SUBSTR isa IceSubstrate
            for k in 1:Nsoil
                Tsoil[k, i, j] = min(Tsoil[k, i, j], Tm)
            end
        end

    end
end
