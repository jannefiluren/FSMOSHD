@kwdef struct FixedConductivity{Tf} <: AbstractConductivity{Tf}
    kfix::Tf = 0.24        # Fixed thermal conductivity of snow (W/m/K)
end

@kwdef struct DensityConductivity{Tf} <: AbstractConductivity{Tf}
    bthr::Tf = 2           # Snow thermal conductivity exponent (-)
end

FixedConductivity{Tf}(grid::Grid; kwargs...) where {Tf} = FixedConductivity{Tf}(; kwargs...)
DensityConductivity{Tf}(grid::Grid; kwargs...) where {Tf} = DensityConductivity{Tf}(; kwargs...)

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
    surface_layer_properties!(i, j, state, diag, grid)

Thickness, temperature and thermal conductivity of the layer the surface energy balance
sees (`diag.Ds1`, `diag.Ts1`, `diag.ks1`) for cell `(i, j)`. The layer is always at least
as thick as the top soil layer and mixes in soil properties for thin snowpacks, so it
requires `Dzsnow[1] >= Dzsoil[1]` (checked in types.jl) - a thinner first snow layer
leaves `Ts1` blended with `Tsoil` even under a deep snowpack.
"""
@inline function surface_layer_properties!(i, j, state, diag, grid)
    (; Dzsoil) = grid
    (; Ds, Tsnow, Tsoil) = state
    (; ksnow, ksoil, Ds1, Ts1, ks1) = diag
    Tf = eltype(Ds1)

    Ds1[i, j] = max(Dzsoil[1], Ds[1, i, j])
    Ts1[i, j] = Tsoil[1, i, j] + (Tsnow[1, i, j] - Tsoil[1, i, j]) * Ds[1, i, j] / Dzsoil[1]

    # Series resistance of the composite layer with (a) a guard against zero division for
    # cells that never held snow and (b) a soil resistance that turns negative once snow
    # fills over half the layer, where ks1 is overridden below
    snow_R = Ds[1, i, j] > zero(Tf) ? Tf(2) * Ds[1, i, j] / ksnow[1, i, j] : zero(Tf)
    soil_R = (Dzsoil[1] - Tf(2) * Ds[1, i, j]) / ksoil[1, i, j]

    ks1[i, j] = Dzsoil[1] / (snow_R + soil_R)
    if (Ds[1, i, j] > Tf(0.5) * Dzsoil[1])
        ks1[i, j] = ksnow[1, i, j]
    end
    if (Ds[1, i, j] > Dzsoil[1])
        Ts1[i, j] = Tsnow[1, i, j]
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
        fsm.state, fsm.diag, fsm.surface, fsm.grid, fsm.params,
        CONDCT, SUBSTR;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function thermal_kernel!(
        state, diag, surface, grid, params::Parameters{Tf},
        CONDCT::AbstractConductivity{Tf}, SUBSTR::AbstractSubstrate{Tf},
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; tthresh) = params
    (; tilefrac) = surface
    (; Tveg) = state
    (; Tveg0) = diag

    if (tilefrac[i, j] >= tthresh)

        snow_conductivity!(CONDCT, i, j, state, diag, params)
        soil_properties!(SUBSTR, i, j, state, diag, surface, grid, params)
        surface_layer_properties!(i, j, state, diag, grid)

        Tveg0[i, j] = Tveg[i, j]

    end
end
