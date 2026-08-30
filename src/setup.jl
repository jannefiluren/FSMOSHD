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
  - `"config"` (optional): Model configuration flags (SNFRAC, EXCHNG, ZOFFST, etc.)
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
    config = copy(get(settings, "config", Dict()))
    params = copy(get(settings, "params", Dict()))
    # EXCHNG/ZOFFST are integer flags that select surface-exchange schemes; consume
    # them here so they are not re-applied as Parameters below.
    EXCHNG = Int(pop!(config, "EXCHNG", 1))
    ZOFFST = Int(pop!(config, "ZOFFST", 0))
    FSNRHO = Int(pop!(config, "FSNRHO", 2))
    DENSTY = Int(pop!(config, "DENSTY", 3))
    HYDROL = Int(pop!(config, "HYDROL", 2))
    SNOLAY = Int(pop!(config, "SNOLAY", 0))
    schemes = (
        ALBEDO = build_scheme(Tf, get(config, "ALBEDO", PrognosticAlbedo), Nx, Ny, params),
        CANOPY = build_scheme(Tf, get(config, "CANOPY",
            settings["tile"] == "forest" ? OneLayerCanopy : NoCanopy), Nx, Ny, params),
        SUBSTR = build_scheme(Tf, get(config, "SUBSTR",
            settings["tile"] == "glacier" ? IceSubstrate : SoilSubstrate), Nx, Ny, params),
        CONDCT = build_scheme(Tf, get(config, "CONDCT", DensityConductivity), Nx, Ny, params),
        reference_height = build_scheme(Tf, ZOFFST == 0 ? AboveGround : AboveCanopy, Nx, Ny, params),
        surface_layer = build_scheme(Tf, EXCHNG == 2 ? ForestSurfaceLayer : OpenSurfaceLayer, Nx, Ny, params),
        stability = build_scheme(Tf, EXCHNG == 1 ? LouisStabilityCorrection : NoStabilityCorrection, Nx, Ny, params),
        FSNRHO = build_scheme(Tf, FSNRHO == 0 ? FixedFreshSnowDensity :
            FSNRHO == 1 ? ClimateFreshSnowDensity : ElevationFreshSnowDensity, Nx, Ny, params),
        COMPACT = build_scheme(Tf, DENSTY == 1 ? AgeCompaction :
            DENSTY == 2 ? OverburdenCompaction :
            DENSTY == 3 ? CrocusCompaction :
            error("DENSTY=$DENSTY is not supported (constant density was removed; use 1, 2, or 3)"),
            Nx, Ny, params),
        HYDROL = build_scheme(Tf, HYDROL == 0 ? FreeDrainingHydrology :
            HYDROL == 1 ? BucketHydrology : DensityBucketHydrology, Nx, Ny, params),
        LAYERING = build_scheme(Tf, SNOLAY == 0 ? OriginalLayering : DensityLayering, Nx, Ny, params),
    )
    fsm = FSM{Tf, Ti}(; Nx = Nx, Ny = Ny, schemes...)

    for scheme in schemes
        check_grid(scheme, Nx, Ny)
    end

    # Tile type is a setup-local input, not stored on the model (Stage 7).
    tile = settings["tile"]

    # Apply config flags and parameter overrides to the right sub-struct.
    apply_config!(fsm, config, schemes)
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


