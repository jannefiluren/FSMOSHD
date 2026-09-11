# Snow compaction parameterizations.

@kwdef struct AgeCompaction{Tf} <: AbstractCompaction{Tf}
    rmlt::Tf = 500              # Maximum density for melting snow (kg/m^3)
    rcld::Tf = 300              # Maximum density for cold snow (kg/m^3)
    trho::Tf = 3600 * 200       # Snow compaction time scale (s)
end

@kwdef struct OverburdenCompaction{Tf} <: AbstractCompaction{Tf}
    eta0::Tf = 3.7e7            # Reference snow viscosity (Pa s)
    snda::Tf = 2.8e-6           # Thermal metamorphism parameter (1/s)
    rhos_max::Tf = 750          # Maximum snow density (kg/m^3)
end

@kwdef struct CrocusCompaction{Tf} <: AbstractCompaction{Tf}
    eta1::Tf = 7.62237e6        # Reference snow viscosity for Crocus B92 compaction (Pa s)
    a_eta::Tf = 0.1             # Temperature factor for Crocus B92 compaction (K^-1)
    b_eta::Tf = 0.023           # First density factor for Crocus B92 compaction (m^3/kg)
    c_eta::Tf = 250             # Second density factor for Crocus B92 compaction (kg/m^3)
    rhos_max::Tf = 750          # Maximum snow density (kg/m^3)
end

AgeCompaction{Tf}(grid::Grid; kwargs...) where {Tf} = AgeCompaction{Tf}(; kwargs...)
OverburdenCompaction{Tf}(grid::Grid; kwargs...) where {Tf} = OverburdenCompaction{Tf}(; kwargs...)
CrocusCompaction{Tf}(grid::Grid; kwargs...) where {Tf} = CrocusCompaction{Tf}(; kwargs...)

"""
    compact_snow!(scheme, i, j, state, params)

Compact the snow column at cell `(i, j)`: rescale the layer thicknesses `Ds` in
place to the compacted density, for every layer.
"""
function compact_snow! end

# Snow compaction with age
@inline function compact_snow!(c::AgeCompaction{Tf}, i, j, state, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, fsnow, Nsnow, Tsnow) = state
    (; dt) = params
    (; rmlt, rcld, trho) = c
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
@inline function compact_snow!(c::OverburdenCompaction{Tf}, i, j, state, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, fsnow, Nsnow, Tsnow) = state
    (; dt) = params
    (; eta0, snda, rhos_max) = c
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
@inline function compact_snow!(c::CrocusCompaction{Tf}, i, j, state, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, fsnow, Nsnow, Tsnow) = state
    (; dt) = params
    (; eta1, a_eta, b_eta, c_eta, rhos_max) = c
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
