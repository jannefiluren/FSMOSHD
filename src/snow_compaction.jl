# Snow compaction parameterizations. Fieldless dispatch schemes: the compaction
# parameters are shared with the fresh-snow-density routine, so they stay on
# Parameters and are passed in.

struct AgeCompaction{Tf} <: AbstractCompaction{Tf} end
struct OverburdenCompaction{Tf} <: AbstractCompaction{Tf} end
struct CrocusCompaction{Tf} <: AbstractCompaction{Tf} end

AgeCompaction{Tf}(grid::Grid; kwargs...) where {Tf} = AgeCompaction{Tf}()
OverburdenCompaction{Tf}(grid::Grid; kwargs...) where {Tf} = OverburdenCompaction{Tf}()
CrocusCompaction{Tf}(grid::Grid; kwargs...) where {Tf} = CrocusCompaction{Tf}()

"""
    compact_snow!(scheme, i, j, state, params)

Compact the snow column at cell `(i, j)`: rescale the layer thicknesses `Ds` in
place to the compacted density, for every layer. A kernel point function (see
`.claude/rules/kernel-point-functions.md`); every `AbstractCompaction`
implements it. The compaction parameters are shared with fresh snow density, so
they stay on `Parameters` and are passed in rather than held on the scheme.
"""
function compact_snow! end

# Snow compaction with age
@inline function compact_snow!(::AgeCompaction{Tf}, i, j, state, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, fsnow, Nsnow, Tsnow) = state
    (; dt, rmlt, rcld, trho) = params
    for k in 1:Nsnow[i, j]
        if (Ds[k, i, j] > eps(Tf))
            rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
            if (Tsnow[k, i, j] >= Tm)
                if (rhos < rmlt)
                    rhos = rmlt + (rhos - rmlt) * exp(-dt / trho)
                end
            else
                if (rhos < rcld)
                    rhos = rcld + (rhos - rcld) * exp(-dt / trho)
                end
            end
            Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rhos / fsnow[i, j]
        end
    end
    return nothing
end

# Snow compaction by overburden
@inline function compact_snow!(::OverburdenCompaction{Tf}, i, j, state, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, fsnow, Nsnow, Tsnow) = state
    (; dt, eta0, snda, rhos_max) = params
    mass = Tf(0.0)
    for k in 1:Nsnow[i, j]
        mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
        if (Ds[k, i, j] > eps(Tf))
            rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
            rhos = rhos + (rhos * grav * mass * dt / (eta0 * exp(-(Tsnow[k, i, j] - Tm) / Tf(12.4) + rhos / Tf(55.6))) + dt * rhos * snda * exp((Tsnow[k, i, j] - Tm) / Tf(23.8) - max(rhos - Tf(150), Tf(0.0)) / Tf(21.7)))
            rhos = min(rhos, rhos_max)
            Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rhos / fsnow[i, j]
        end
        mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
    end
    return nothing
end

# Snow compaction by overburden, dependent on liquid water content (Crocus B92)
@inline function compact_snow!(::CrocusCompaction{Tf}, i, j, state, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, fsnow, Nsnow, Tsnow) = state
    (; dt, eta1, a_eta, b_eta, c_eta, rhos_max) = params
    mass = Tf(0.0)
    for k in 1:Nsnow[i, j]
        mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
        if (Ds[k, i, j] > eps(Tf))
            rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
            f1 = Tf(1) / (Tf(1) + Tf(600) * Sliq[k, i, j] / (rho_wat * Ds[k, i, j] * fsnow[i, j]))
            f2 = Tf(1.0)
            eta = f1 * f2 * eta1 * (rhos / c_eta) * exp(a_eta * (Tm - Tsnow[k, i, j]) + b_eta * rhos)
            rhos = rhos + rhos * grav * mass * dt / eta
            rhos = min(rhos, rhos_max)
            Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rhos / fsnow[i, j]
        end
        mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
    end
    return nothing
end
