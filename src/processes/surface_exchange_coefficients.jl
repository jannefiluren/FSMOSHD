# Computation of turbulent eddy diffusivities for heat and moisture

"""
    surface_exchange_coefficients!(fsm, meteo)

Eddy diffusivities for turbulent transfer of heat and moisture between the ground,
the canopy and the atmosphere.

# Arguments
- `fsm::FSM`: Model state structure
- `meteo::MET`: Current meteorological conditions
"""
function surface_exchange_coefficients!(fsm::FSM{Tf}, meteo::MET{Tf}) where {Tf <: Real}

    (; reference_height, surface_layer, snow_fraction) = fsm.physics

    backend = get_backend(fsm.state.Tsrf)
    kernel! = surface_exchange_coefficients_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.params, meteo,
        reference_height, surface_layer, snow_fraction;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function surface_exchange_coefficients_kernel!(
        state, diag, surface, params::Parameters{Tf}, meteo,
        reference_height::AbstractReferenceHeight{Tf},
        surface_layer::AbstractSurfaceLayer{Tf},
        snow_fraction::AbstractSnowFraction{Tf},
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; tthresh, zT, zU) = params
    (; tilefrac, hcan) = surface

    if (tilefrac[i, j] >= tthresh)

        zU1, zT1 = reference_heights(reference_height, zU, zT, hcan[i, j])
        z0g = ground_roughness(snow_fraction, i, j, state, surface)

        exchange_coefficients!(surface_layer, i, j, state, diag, surface, params, meteo, zU1, zT1, z0g)

    end
end
