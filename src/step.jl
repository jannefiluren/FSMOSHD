"""
    step!(fsm, met, t)

Execute one complete physics time step of the snow model.

This function encapsulates the standard model execution sequence:
1. Meteorological data processing (drive!)
2. Radiation calculations  
3. Thermal property updates
4. Iterative energy balance (tile-specific)
5. Canopy processes (forest tiles only)
6. Snow processes
7. Soil thermal processes

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `met::MET`: Current meteorological conditions (read-only, never modified)
- `t::DateTime`: Current simulation time

# Example
```julia
fsm = setup(Float32, Int32, landuse, Nx, Ny, TILE="forest")  
met = MET{Float32, Int32}(Nx=Nx, Ny=Ny)
step!(fsm, met, DateTime(2023, 12, 1, 12))
```
"""
function step!(fsm::FSM{Tf, Ti}, met::MET{Tf, Ti}, t) where {Tf, Ti}

    # 1. Meteorological data processing
    drive!(fsm, met)

    # 2. Radiation calculations
    radiation!(fsm, met, t)

    # 3. Thermal property updates
    thermal!(fsm)

    # 4. Iterative energy balance solution
    tile_type = fsm.TILE
    for _ in 1:fsm.Nitr
        sfexch!(fsm, met)
        if tile_type == "forest"
            ebalfor!(fsm, met)
        else
            ebalsrf!(fsm, met)
        end
    end

    # 5. Forest-specific canopy processing
    if tile_type == "forest"
        canopy!(fsm, met)
    end

    # 6. Snow processes
    snow!(fsm, met, t)

    # 7. Soil thermal processes
    soil!(fsm)

    return nothing
end
