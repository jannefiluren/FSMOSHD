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
  - `"config"` (optional): Model configuration flags (CANMOD, SNFRAC, EXCHNG, ZOFFST, etc.)
  - `"params"` (optional): Parameter overrides (hfsn etc.)

# Returns
- `FSM`: Initialized model state structure ready for simulation
"""
function setup end

function setup(Tf, Ti, landuse::Dict, Nx::Int, Ny::Int, settings::Dict)
    return setup(CPU(), Tf, Ti, landuse, Nx, Ny, settings)
end

"""
    build_scheme_from_flag(process, requested, Tf, Nx, Ny, params)

Build the parameterization a configuration entry asks for. `requested` is either an integer
flag - translated to a scheme type through the table below, which is the FSM2oshd naming - or
a scheme type/instance, which `build_scheme` takes unchanged.

The integer-flag layer exists to keep `config` aligned with OSHDinternal while the two are
tested against each other; it is meant to go away once `config` names scheme types directly.
"""
function build_scheme_from_flag(process, requested, Tf, Nx, Ny, params)

    requested isa Number || return build_scheme(Tf, requested, Nx, Ny, params)

    lookup = Dict(
        "ALBEDO" => Dict(
            0 => DiagnosticAlbedo,
            1 => DecayAlbedo,
            2 => PrognosticAlbedo,
            ),
        "CANMOD" => Dict(
            0 => NoCanopy,
            1 => OneLayerCanopy,
            ),
        "CONDCT" => Dict(
            0 => FixedConductivity,
            1 => DensityConductivity,
            ),
        "DENSTY" => Dict(
            1 => AgeCompaction,
            2 => OverburdenCompaction,
            3 => CrocusCompaction,
            ),
        "EXCHNG" => Dict(
            0 => NoStabilityCorrection,
            1 => LouisStabilityCorrection,
            ),
        "HYDROL" => Dict(
            0 => FreeDrainingHydrology,
            1 => BucketHydrology,
            2 => DensityBucketHydrology,
            ),
        "SNFRAC" => Dict(
            0 => SeasonalSnowFraction,
            1 => HelbigSnowFraction,
            2 => HelbigMaxSnowFraction,
            3 => PointSnowFraction,
            4 => TanhSnowFraction,
            ),
        "ZOFFST" => Dict(
            0 => AboveGround,
            1 => AboveCanopy,
            ),
        "FSNRHO" => Dict(
            0 => FixedFreshSnowDensity,
            1 => ClimateFreshSnowDensity,
            2 => ElevationFreshSnowDensity,
            ),
        "SNOLAY" => Dict(
            0 => OriginalLayering,
            1 => DensityLayering,
            ),
        )

    flags = lookup[process]
    flag = Int(requested)
    haskey(flags, flag) ||
        error("$process=$flag is not supported (use $(join(sort!(collect(keys(flags))), ", ")))")

    return build_scheme(Tf, flags[flag], Nx, Ny, params)

end

function setup(arch::AbstractArchitecture, Tf, Ti, landuse::Dict, Nx::Int, Ny::Int, settings::Dict)

    @unpack_constants(Tf)

    config = copy(get(settings, "config", Dict()))
    params = copy(get(settings, "params", Dict()))

    tile = settings["tile"]
    tile in ("open", "forest", "glacier") || error("tile requires open, forest or glacier (got tile = $tile)")

    # Default configuration
    ALBEDO = pop!(config, "ALBEDO", 2)
    CANMOD = pop!(config, "CANMOD", 0)
    CONDCT = pop!(config, "CONDCT", 1)
    DENSTY = pop!(config, "DENSTY", 3)
    EXCHNG = pop!(config, "EXCHNG", 1)
    HYDROL = pop!(config, "HYDROL", 2)
    SNFRAC = pop!(config, "SNFRAC", 3)
    ZOFFST = pop!(config, "ZOFFST", 0)
    FSNRHO = pop!(config, "FSNRHO", 2)
    SNOLAY = pop!(config, "SNOLAY", 0)

    canopy = build_scheme_from_flag("CANMOD", CANMOD, Tf, Nx, Ny, params)

    # Validate combinations of configurations
    if tile == "forest"
        canopy isa OneLayerCanopy || error("forest tile requires CANMOD == 1 (got CANMOD = $CANMOD)")
        EXCHNG == 2 || error("forest tile requires EXCHNG == 2 (got EXCHNG = $EXCHNG)")
    else
        EXCHNG in (0, 1) || error("open/glacier tile requires EXCHNG 0 or 1 (got EXCHNG = $EXCHNG)")
    end

    # Define surface and substrate layer given tile class
    surface_layer, substrate_layer = if tile == "forest"
        (build_scheme(Tf, ForestSurfaceLayer, Nx, Ny, params),
         build_scheme(Tf, SoilSubstrate, Nx, Ny, params))
    else
        stability = build_scheme_from_flag("EXCHNG", EXCHNG, Tf, Nx, Ny, params)
        (OpenSurfaceLayer{Tf}(; stability = stability),
         build_scheme(Tf, tile == "open" ? SoilSubstrate : IceSubstrate, Nx, Ny, params))
    end

    schemes = (
        surface_layer    = surface_layer,
        SUBSTR           = substrate_layer,
        ALBEDO           = build_scheme_from_flag("ALBEDO", ALBEDO, Tf, Nx, Ny, params),
        CANOPY           = canopy,
        CONDCT           = build_scheme_from_flag("CONDCT", CONDCT, Tf, Nx, Ny, params),
        COMPACT          = build_scheme_from_flag("DENSTY", DENSTY, Tf, Nx, Ny, params),
        HYDROL           = build_scheme_from_flag("HYDROL", HYDROL, Tf, Nx, Ny, params),
        SNFRAC           = build_scheme_from_flag("SNFRAC", SNFRAC, Tf, Nx, Ny, params),
        reference_height = build_scheme_from_flag("ZOFFST", ZOFFST, Tf, Nx, Ny, params),
        FSNRHO           = build_scheme_from_flag("FSNRHO", FSNRHO, Tf, Nx, Ny, params),
        LAYERING         = build_scheme_from_flag("SNOLAY", SNOLAY, Tf, Nx, Ny, params),
    )

    fsm = FSM{Tf, Ti}(; Nx = Nx, Ny = Ny, schemes...)

    for scheme in schemes
        check_grid(scheme, Nx, Ny)
    end

    # Apply config flags and parameter overrides to the right sub-struct.
    apply_config!(fsm, config)
    apply_params!(fsm, params)

    lu = fsm.landuse
    st = fsm.state

    # Settings specific for fixed fresh snow density
    if fsm.physics.FSNRHO isa FixedFreshSnowDensity
        fsm.params = reconstruct(fsm.params; rhof = fsm.params.rho0)
    end

    # Derived soil parameters
    mask = lu.fcly .+ lu.fsnd .> Tf(1)
    lu.fcly[mask] .= Tf(1) .- lu.fsnd[mask]

    lu.b .= Tf(3.1) .+ Tf(15.7) .* lu.fcly .- Tf(0.3) .* lu.fsnd
    lu.hcap_soil .= (Tf(2.128) .* lu.fcly .+ Tf(2.385) .* lu.fsnd) .* Tf(1.0e6) ./ (lu.fcly .+ lu.fsnd)
    lu.sathh .= Tf(10) .^ (Tf(0.17) .- Tf(0.63) .* lu.fcly .- Tf(1.58) .* lu.fsnd)
    lu.Vsat .= Tf(0.505) .- Tf(0.037) .* lu.fcly .- Tf(0.142) .* lu.fsnd
    lu.Vcrit .= lu.Vsat .* (lu.sathh ./ Tf(3.364)) .^ (Tf(1) ./ lu.b)
    hcon_min = (hcon_clay .^ lu.fcly) .* (hcon_sand .^ (Tf(1) .- lu.fcly))
    lu.hcon_soil .= (hcon_air .^ lu.Vsat) .* (hcon_min .^ (Tf(1) .- lu.Vsat))

    # Initial soil profiles
    for k in 1:fsm.grid.Nsoil
        st.theta[k, :, :] .= fsm.params.fsat * lu.Vsat[:, :]
        st.Tsoil[k, :, :] .= fsm.params.Tprof
    end

    # Cap surface and soil temperatures for glacier
    if fsm.physics.SUBSTR isa IceSubstrate
        st.Tsrf .= min.(st.Tsrf, Tm)
        st.Tsoil .= min.(st.Tsoil, Tm)
    end

    # Load terrain properties from landuse data
    lu.fsky_terr .= Tf.(landuse["skyvf"]["data"])
    lu.dem .= Tf.(landuse["elevation"]["data"])
    lu.prec_multi .= landuse["prec_multi"]["data"]   # TODO hack float64

    # Set tile fractions non open tiles
    if tile != "open"
        lu.tilefrac .= Tf.(landuse[lowercase(tile)]["data"])
    end

    # Initialize snow cover fraction specific variables
    lu.slopemu .= Tf.(landuse["slopemu"]["data"])
    lu.xi .= Tf.(landuse["xi"]["data"])
    lu.Ld .= Tf.(landuse["Ld"]["data"])

    # Canopy properties
    if tile == "forest"

        lu.fveg .= Tf.(landuse["fveg"]["data"])
        lu.hcan .= Tf.(landuse["hcan"]["data"])
        lu.lai .= Tf.(landuse["lai"]["data"])
        lu.vfhp .= Tf.(landuse["vfhp"]["data"])
        lu.fves .= Tf.(landuse["fves"]["data"])

        lu.pmultf .= Tf.((1 .- (1 .- landuse["prec_multi"]["data"]) .* (1 .- landuse["forest"]["data"] * fsm.params.pmultf_for)) ./ landuse["prec_multi"]["data"])   # TODO if this works, integrate with prec_multi instead...

        lu.VAI[:, :] = lu.lai[:, :]
        lu.trcn[:, :] = Tf(1) .- Tf(0.9) .* lu.fveg[:, :]
        lu.fsky .= lu.vfhp ./ lu.trcn
        # Handle values where fsky > 1
        mask = lu.fsky .> Tf(1)
        lu.trcn[mask] .= lu.vfhp[mask]
        lu.fsky[mask] .= Tf(1)
    end

    # Narrow the tile mask by the configuration's data requirement (canopy tile: fveg > 0).
    if tile == "forest"
        canopy_free = (lu.tilefrac .>= fsm.params.tthresh) .& (lu.fveg .<= 0)
        dropped = count(canopy_free)
        if dropped > 0
            @warn "forest tile: $dropped active cell(s) have fveg == 0 and are excluded from the tile"
            lu.tilefrac[canopy_free] .= Tf(0)
        end
    end

    lu.canh[:, :] = Tf(12500) * lu.VAI[:, :]
    lu.scap[:, :] = fsm.params.cvai * lu.VAI[:, :]

    # The whole setup above runs on the CPU (it uses scalar indexing); the
    # finished structure is moved to the target architecture in one step.
    if !(arch isa CPU)
        fsm = on_architecture(arch, fsm)
    end

    return fsm

end


