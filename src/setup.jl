"""
    setup([arch], Tf, Ti, landuse, Nx, Ny, settings)

Initialize the FSM snow model with specified configuration and domain properties.

Creates and configures the FSM model state structure with landuse data, model parameters,
and domain dimensions. Sets up soil properties, surface characteristics, and applies
configuration-specific settings for different surface types and model behaviors.

# Arguments
- `arch::AbstractArchitecture`: Architecture the model arrays live on
  (optional, default `CPU()`); pass e.g. `GPU(CUDABackend())` for GPU runs
- `Tf`: Floating-point precision type (typically Float32 or Float64)
- `Ti`: Integer type for array indices (typically Int32 or Int64)
- `landuse::Dict`: Landuse data dictionary with topographic and surface properties
- `Nx::Int, Ny::Int`: Model domain dimensions
- `settings::Dict`: Configuration dictionary containing:
  - `"tile"`: Surface tile type ("open", "forest", "glacier")
  - `"config"` (optional): Model configuration flags (SNFRAC, CANMOD, EXCHNG, ZOFFST, etc.)
  - `"params"` (optional): Parameter overrides (hfsn etc.)

# Returns
- `FSM`: Initialized model state structure ready for simulation
"""
function setup(Tf, Ti, landuse::Dict, Nx::Int, Ny::Int, settings::Dict)
    return setup(CPU(), Tf, Ti, landuse, Nx, Ny, settings)
end

function setup(arch::AbstractArchitecture, Tf, Ti, landuse::Dict, Nx::Int, Ny::Int, settings::Dict)

    @unpack_constants(Tf)

    # Create fsm object. Parameterizations are type parameters, so they must be
    # chosen at construction: setfield! cannot change a field's type afterwards.
    config = get(settings, "config", Dict())
    params = copy(get(settings, "params", Dict()))
    schemes = (
        ALBEDO = build_scheme(Tf, get(config, "ALBEDO", PrognosticAlbedo), Nx, Ny, params),
        CONDCT = build_scheme(Tf, get(config, "CONDCT", DensityConductivity), Nx, Ny, params),
    )
    fsm = FSM{Tf, Ti}(; Nx = Nx, Ny = Ny, schemes...)

    for scheme in schemes
        check_grid(scheme, Nx, Ny)
    end

    # Set tile
    fsm.TILE = settings["tile"]

    # Apply model configuration
    for (key, value) in config
        haskey(schemes, Symbol(key)) && continue   # already applied at construction
        if value isa Int
            value = Ti(value)
        end
        setfield!(fsm, Symbol(key), value)
    end

    # Apply parameter overrides
    for (key, value) in params
        field = Symbol(key)
        existing = getfield(fsm, field)
        target_type = eltype(existing)
        if existing isa AbstractArray && !(value isa AbstractArray)
            # scalar overriding an array field: fill the entire array
            fill!(existing, target_type(value))
        else
            # scalar to scalar, or array to array: assign directly
            setfield!(fsm, field, target_type.(value))
        end
    end

    # Settings specific for FSNRHO=0 (fixed fresh snow density)
    if (fsm.FSNRHO == 0)
        fsm.rhof = fsm.rho0
    end

    # Derived soil parameters
    mask = fsm.fcly .+ fsm.fsnd .> Tf(1)
    fsm.fcly[mask] .= Tf(1) .- fsm.fsnd[mask]

    fsm.b .= Tf(3.1) .+ Tf(15.7) .* fsm.fcly .- Tf(0.3) .* fsm.fsnd
    fsm.hcap_soil .= (Tf(2.128) .* fsm.fcly .+ Tf(2.385) .* fsm.fsnd) .* Tf(1.0e6) ./ (fsm.fcly .+ fsm.fsnd)
    fsm.sathh .= Tf(10) .^ (Tf(0.17) .- Tf(0.63) .* fsm.fcly .- Tf(1.58) .* fsm.fsnd)
    fsm.Vsat .= Tf(0.505) .- Tf(0.037) .* fsm.fcly .- Tf(0.142) .* fsm.fsnd
    fsm.Vcrit .= fsm.Vsat .* (fsm.sathh ./ Tf(3.364)) .^ (Tf(1) ./ fsm.b)
    hcon_min = (hcon_clay .^ fsm.fcly) .* (hcon_sand .^ (Tf(1) .- fsm.fcly))
    fsm.hcon_soil .= (hcon_air .^ fsm.Vsat) .* (hcon_min .^ (Tf(1) .- fsm.Vsat))

    # Initial soil profiles
    for k in 1:fsm.Nsoil
        fsm.theta[k, :, :] .= fsm.fsat * fsm.Vsat[:, :]
        fsm.Tsoil[k, :, :] .= fsm.Tprof
    end

    # Cap surface and soil temperatures for glacier
    if (fsm.TILE == "glacier")
        fsm.Tsrf .= min.(fsm.Tsrf, Tm)
        fsm.Tsoil .= min.(fsm.Tsoil, Tm)
    end

    # Load terrain properties from landuse data
    fsm.fsky_terr .= Tf.(landuse["skyvf"]["data"])
    fsm.dem .= Tf.(landuse["elevation"]["data"])
    fsm.prec_multi .= landuse["prec_multi"]["data"]   # TODO hack float64

    # Set tile fractions non open tiles
    if (fsm.TILE != "open")
        fsm.tilefrac .= Tf.(landuse[lowercase(fsm.TILE)]["data"])
    end

    # Initialize snow cover fraction specific variables
    fsm.slopemu .= Tf.(landuse["slopemu"]["data"])
    fsm.xi .= Tf.(landuse["xi"]["data"])
    fsm.Ld .= Tf.(landuse["Ld"]["data"])

    # Canopy properties
    if (fsm.TILE == "forest")

        fsm.fveg .= Tf.(landuse["fveg"]["data"])
        fsm.hcan .= Tf.(landuse["hcan"]["data"])
        fsm.lai .= Tf.(landuse["lai"]["data"])
        fsm.vfhp .= Tf.(landuse["vfhp"]["data"])
        fsm.fves .= Tf.(landuse["fves"]["data"])

        fsm.pmultf .= Tf.((1 .- (1 .- landuse["prec_multi"]["data"]) .* (1 .- landuse["forest"]["data"] * fsm.pmultf_for)) ./ landuse["prec_multi"]["data"])   # TODO if this works, integrate with prec_multi instead...

        fsm.VAI[:, :] = fsm.lai[:, :]
        fsm.trcn[:, :] = Tf(1) .- Tf(0.9) .* fsm.fveg[:, :]
        fsm.fsky .= fsm.vfhp ./ fsm.trcn
        # Handle values where fsky > 1
        mask = fsm.fsky .> Tf(1)
        fsm.trcn[mask] .= fsm.vfhp[mask]
        fsm.fsky[mask] .= Tf(1)
    end

    fsm.canh[:, :] = Tf(12500) * fsm.VAI[:, :]
    fsm.scap[:, :] = fsm.cvai * fsm.VAI[:, :]

    # The whole setup above runs on the CPU (it uses scalar indexing); the
    # finished structure is moved to the target architecture in one step.
    if !(arch isa CPU)
        fsm = on_architecture(arch, fsm)
    end

    return fsm

end


