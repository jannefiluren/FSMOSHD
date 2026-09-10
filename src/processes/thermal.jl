# Thermal property calculations for snow and soil layers

"""
    thermal!(fsm)

Thermal property calculations for snow and soil layers.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
"""
function thermal!(fsm::FSM{Tf}) where {Tf <: Real}

    (; conductivity, substrate) = fsm.physics

    backend = get_backend(fsm.diag.gs1)
    kernel! = thermal_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.grid, fsm.params,
        conductivity, substrate;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function thermal_kernel!(
        state, diag, surface, grid, params::Parameters{Tf},
        conductivity::AbstractConductivity{Tf}, substrate::AbstractSubstrate{Tf},
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; tthresh) = params
    (; tilefrac) = surface
    (; Tveg) = state
    (; Tveg0) = diag

    if (tilefrac[i, j] >= tthresh)

        snow_conductivity!(conductivity, i, j, state, diag, params)
        soil_properties!(substrate, i, j, state, diag, surface, grid, params)
        surface_layer_properties!(i, j, state, diag, grid)

        Tveg0[i, j] = Tveg[i, j]

    end
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
