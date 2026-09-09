# Canopy process: snow interception, sublimation and unloading

"""
    canopy!(fsm)

Snow interception, sublimation, and unloading from the vegetation canopy.

# Arguments
- `fsm::FSM`: Model state structure
"""
function canopy!(fsm::FSM{Tf}) where {Tf <: Real}

    (; canopy) = fsm.physics

    backend = get_backend(fsm.state.Sveg)
    kernel! = canopy_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.params,
        canopy;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function canopy_kernel!(
        state, diag, surface, params::Parameters{Tf},
        canopy::AbstractCanopy{Tf},
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; tthresh) = params
    (; tilefrac) = surface

    if (tilefrac[i, j] >= tthresh)

        canopy_snow!(canopy, i, j, state, diag, surface, params)

    end
end
