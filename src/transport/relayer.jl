# Grid-level relayering pass used only by the snow-transport step (transport!): after transport
# deposits redistributed snow into diag, this layers it into the pack by reusing snow_layering!.

"""
    relayer!(fsm, met, t; update_hist = false)

Grid-level relayering pass: accumulate the deposit currently in `diag.snowdepth0` /
`diag.Sice0` into the snowpack, update the snow cover fraction and relayer, at every cell above
the tile threshold. Reuses the [`snow_layering!`](@ref) point function that `snow_kernel!` runs
for new snow; [`transport!`](@ref) calls this to layer in redistributed snow. `update_hist`
should be `false` here so the 14-day history is rolled only once per step (by `snow!`).
"""
function relayer!(fsm::FSM{Tf}, met::MET{Tf}, t; update_hist::Bool = false) where {Tf}

    (; Nsmax) = fsm.grid

    backend = get_backend(fsm.state.Tsnow)
    kernel! = relayer_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.grid, fsm.params, met,
        fsm.physics.layering, fsm.physics.snow_fraction, update_hist, Val(Int(Nsmax));
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel inbounds = true function relayer_kernel!(
        state, diag, surface, grid, params::Parameters{Tf}, meteo,
        layering::AbstractLayering{Tf}, snow_fraction::AbstractSnowFraction{Tf},
        update_hist::Bool, ::Val{Nsmax},
    ) where {Tf, Nsmax}

    i, j = @index(Global, NTuple)

    (; tthresh) = params
    (; tilefrac) = surface

    if (tilefrac[i, j] >= tthresh)
        snow_layering!(
            layering, snow_fraction, i, j, state, diag, surface, grid, params, meteo,
            update_hist, Val(Nsmax)
        )
    end
end
