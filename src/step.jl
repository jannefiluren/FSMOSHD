"""
    step!(fsm, met, t)

Execute one complete physics time step of the snow model.

This function encapsulates the standard model execution sequence:
1. Meteorological data processing
2. Radiation calculations
3. Thermal property updates
4. Iterative energy balance
5. Canopy processes
6. Snow processes
7. Horizontal snow transport
8. Soil thermal processes

# Arguments
- `fsm::FSM`: Model state structure
- `met::MET`: Current meteorological conditions
- `t::DateTime`: Current simulation time
- `transport::Union{SnowTransport, Nothing}` (keyword): when a workspace is passed, run the
  horizontal snow-transport step ([`transport!`](@ref)) after `snow!`; when `nothing` (default),
  the step is identical to a run without transport. Transport is CPU-only.

# Example
```julia
fsm = setup(Float32, surface, Nx, Ny, Dict("tile" => "forest"))
met = MET{Float32}(Nx = Nx, Ny = Ny)
step!(fsm, met, DateTime(2023, 12, 1, 12))
```
"""
function step!(fsm::FSM{Tf}, met::MET{Tf}, t; transport = nothing) where {Tf}

    # 1. Meteorological data processing
    drive!(fsm, met)

    # 2. Radiation calculations
    radiation!(fsm, met, t)

    # 3. Thermal property updates
    thermal!(fsm)

    # 4. Iterative energy balance solution
    for _ in 1:fsm.params.Nitr
        surface_exchange_coefficients!(fsm, met)
        surface_energy_balance!(fsm, met)
    end

    # 5. Canopy interception / unloading
    canopy!(fsm)

    # 6. Snow processes
    snow!(fsm, met, t)

    # 7 Horizontal snow transport
    transport === nothing || transport!(fsm, met, transport, t)

    # 8. Soil thermal processes
    soil!(fsm)

    return nothing
end
