# Model construction: the setup routine that assembles an FSM from a grid, a landuse domain and
# a physics/params selection, followed by the construction helpers it uses.

"""
    setup([arch], grid, landuse; tile, physics = Dict(), params = Dict())

Initialize the FSM snow model on `grid` for a landuse domain.

Builds and configures the model state, selecting physics parameterizations from `physics`
(falling back to defaults) and applying `params` overrides. Float precision and domain size are
taken from `grid` (`eltype(grid)`, `grid.Nx`, `grid.Ny`).

# Arguments
- `arch::AbstractArchitecture`: architecture the model arrays live on (optional, default `CPU()`);
  pass e.g. `GPU(CUDABackend())` for GPU runs.
- `grid::Grid`: model grid, carrying float precision and `Nx`/`Ny`.
- `landuse::Dict`: landuse data with topographic and surface properties.

# Keyword arguments
- `tile`: surface tile type (`"open"`, `"forest"`, `"glacier"`).
- `physics::Dict` (optional): scheme selection, mapping a physics name to a scheme type or
  instance, e.g. `Dict("snow_fraction" => TanhSnowFraction, "canopy" => OneLayerCanopy)`. Keys:
  `snow_albedo`, `canopy`, `substrate`, `conductivity`, `compaction`, `hydrology`,
  `new_snow_density`, `layering`, `snow_fraction`, `reference_height`, `surface_layer`.
  Unspecified schemes use defaults; `canopy`/`surface_layer`/`substrate` default from `tile`.
- `params::Dict` (optional): parameter overrides routed by field name to `Parameters` (scalars)
  or `Surface` (per-cell fields), e.g. `dt`, `z0_snow`. Scheme parameters are set at construction
  via the `physics` schemes, not here.

# Returns
- `FSM`: initialized model state ready for simulation.
"""
function setup end

function setup(grid::Grid, landuse::Dict; kwargs...)
    return setup(CPU(), grid, landuse; kwargs...)
end

# Convenience: a settings dict with "tile" (required) and optional "physics"/"params" keys.
function setup(arch::AbstractArchitecture, grid::Grid, landuse::Dict, settings::AbstractDict)
    return setup(
        arch, grid, landuse;
        tile = settings["tile"],
        physics = get(settings, "physics", Dict()),
        params = get(settings, "params", Dict()),
    )
end

setup(grid::Grid, landuse::Dict, settings::AbstractDict) = setup(CPU(), grid, landuse, settings)

const PHYSICS_KEYS = (
    "snow_albedo", "canopy", "substrate", "conductivity", "compaction", "hydrology",
    "new_snow_density", "layering", "snow_fraction", "reference_height", "surface_layer",
)

function setup(
        arch::AbstractArchitecture,
        grid::Grid,
        landuse::Dict;
        tile::String,
        physics::AbstractDict = Dict(),
        params::AbstractDict = Dict(),
    )

    Tf = eltype(grid)
    Nx, Ny = grid.Nx, grid.Ny
    @unpack_constants(Tf)

    tile in ("open", "forest", "glacier") || error("tile requires open, forest or glacier (got tile = $tile)")

    for key in keys(physics)
        key in PHYSICS_KEYS || throw(ArgumentError("unknown physics key \"$key\" (known: $(join(PHYSICS_KEYS, ", ")))"))
    end

    # canopy / surface_layer / substrate default from the tile; the user may override via physics.
    canopy = instantiate(get(physics, "canopy", tile == "forest" ? OneLayerCanopy : NoCanopy), grid)
    if tile == "forest"
        canopy isa OneLayerCanopy || error("forest tile requires a OneLayerCanopy canopy")
        default_surface_layer = ForestSurfaceLayer
        default_substrate = SoilSubstrate
    else
        default_surface_layer = OpenSurfaceLayer{Tf}(; stability = LouisStabilityCorrection{Tf}())
        default_substrate = tile == "open" ? SoilSubstrate : IceSubstrate
    end

    schemes = (
        snow_albedo        = instantiate(get(physics, "snow_albedo", PrognosticAlbedo), grid),
        canopy             = canopy,
        substrate          = instantiate(get(physics, "substrate", default_substrate), grid),
        conductivity       = instantiate(get(physics, "conductivity", DensityConductivity), grid),
        compaction         = instantiate(get(physics, "compaction", CrocusCompaction), grid),
        hydrology          = instantiate(get(physics, "hydrology", DensityBucketHydrology), grid),
        new_snow_density   = instantiate(get(physics, "new_snow_density", ElevationFreshSnowDensity), grid),
        layering           = instantiate(get(physics, "layering", OriginalLayering), grid),
        snow_fraction      = instantiate(get(physics, "snow_fraction", PointSnowFraction), grid),
        reference_height   = instantiate(get(physics, "reference_height", AboveGround), grid),
        surface_layer      = instantiate(get(physics, "surface_layer", default_surface_layer), grid),
    )

    fsm = FSM(grid; schemes...)

    for scheme in schemes
        check_grid(scheme, Nx, Ny)
    end

    # Apply parameter overrides to the right sub-struct.
    apply_params!(fsm, params)

    sf = fsm.surface
    st = fsm.state

    # Settings specific for fixed fresh snow density
    if fsm.physics.new_snow_density isa FixedFreshSnowDensity
        fsm.params = reconstruct(fsm.params; rhof = fsm.params.rho0)
    end

    # Derived soil parameters
    mask = sf.fcly .+ sf.fsnd .> Tf(1)
    sf.fcly[mask] .= Tf(1) .- sf.fsnd[mask]

    sf.b .= Tf(3.1) .+ Tf(15.7) .* sf.fcly .- Tf(0.3) .* sf.fsnd
    sf.hcap_soil .= (Tf(2.128) .* sf.fcly .+ Tf(2.385) .* sf.fsnd) .* Tf(1.0e6) ./ (sf.fcly .+ sf.fsnd)
    sf.sathh .= Tf(10) .^ (Tf(0.17) .- Tf(0.63) .* sf.fcly .- Tf(1.58) .* sf.fsnd)
    sf.Vsat .= Tf(0.505) .- Tf(0.037) .* sf.fcly .- Tf(0.142) .* sf.fsnd
    sf.Vcrit .= sf.Vsat .* (sf.sathh ./ Tf(3.364)) .^ (Tf(1) ./ sf.b)
    hcon_min = (hcon_clay .^ sf.fcly) .* (hcon_sand .^ (Tf(1) .- sf.fcly))
    sf.hcon_soil .= (hcon_air .^ sf.Vsat) .* (hcon_min .^ (Tf(1) .- sf.Vsat))

    # Initial soil profiles
    for k in 1:fsm.grid.Nsoil
        st.theta[k, :, :] .= fsm.params.fsat * sf.Vsat[:, :]
        st.Tsoil[k, :, :] .= fsm.params.Tprof
    end

    # Cap surface and soil temperatures for glacier
    if fsm.physics.substrate isa IceSubstrate
        st.Tsrf .= min.(st.Tsrf, Tm)
        st.Tsoil .= min.(st.Tsoil, Tm)
    end

    # Load terrain properties from landuse data
    sf.fsky_terr .= Tf.(landuse["skyvf"]["data"])
    sf.dem .= Tf.(landuse["elevation"]["data"])
    sf.prec_multi .= landuse["prec_multi"]["data"]   # TODO hack float64

    # Set tile fractions non open tiles
    if tile != "open"
        sf.tilefrac .= Tf.(landuse[lowercase(tile)]["data"])
    end

    # Initialize snow cover fraction specific variables
    sf.slopemu .= Tf.(landuse["slopemu"]["data"])
    sf.xi .= Tf.(landuse["xi"]["data"])
    sf.Ld .= Tf.(landuse["Ld"]["data"])

    # Canopy properties
    if tile == "forest"

        sf.fveg .= Tf.(landuse["fveg"]["data"])
        sf.hcan .= Tf.(landuse["hcan"]["data"])
        sf.lai .= Tf.(landuse["lai"]["data"])
        sf.vfhp .= Tf.(landuse["vfhp"]["data"])
        sf.fves .= Tf.(landuse["fves"]["data"])

        sf.pmultf .= Tf.((1 .- (1 .- landuse["prec_multi"]["data"]) .* (1 .- landuse["forest"]["data"] * fsm.params.pmultf_for)) ./ landuse["prec_multi"]["data"])   # TODO if this works, integrate with prec_multi instead...

        sf.VAI[:, :] = sf.lai[:, :]
        sf.trcn[:, :] = Tf(1) .- Tf(0.9) .* sf.fveg[:, :]
        sf.fsky .= sf.vfhp ./ sf.trcn
        # Handle values where fsky > 1
        mask = sf.fsky .> Tf(1)
        sf.trcn[mask] .= sf.vfhp[mask]
        sf.fsky[mask] .= Tf(1)
    end

    # Narrow the tile mask by the configuration's data requirement (canopy tile: fveg > 0).
    if tile == "forest"
        canopy_free = (sf.tilefrac .>= fsm.params.tthresh) .& (sf.fveg .<= 0)
        dropped = count(canopy_free)
        if dropped > 0
            @warn "forest tile: $dropped active cell(s) have fveg == 0 and are excluded from the tile"
            sf.tilefrac[canopy_free] .= Tf(0)
        end
    end

    sf.canh[:, :] = Tf(12500) * sf.VAI[:, :]
    sf.scap[:, :] = fsm.params.cvai * sf.VAI[:, :]

    # The whole setup above runs on the CPU (it uses scalar indexing); the
    # finished structure is moved to the target architecture in one step.
    if !(arch isa CPU)
        fsm = on_architecture(arch, fsm)
    end

    return fsm

end

# --------------------------------------------------------------------------------------------
# Construction helpers used by setup: scheme instantiation, grid checks, parameter overrides.

"""
    grid_array(Tf, x, Nx, Ny)

Materialize a parameterization parameter as an `Nx` by `Ny` array of element type `Tf`.
A scalar is broadcast over the whole grid; an array is converted element-wise.
"""
grid_array(::Type{Tf}, x::Number, Nx, Ny) where {Tf} = fill(Tf(x), Nx, Ny)
grid_array(::Type{Tf}, x::AbstractArray, Nx, Ny) where {Tf} = convert(Array{Tf, 2}, x)

"""
    check_grid(scheme, Nx, Ny)

Assert that any grid-shaped parameter held by `scheme` matches the `Nx` by `Ny` model grid.
The fallback accepts anything, so a parameterization built only from scalars needs no method.

Without this a mismatch is not caught at setup; it surfaces later as a `BoundsError` from
inside a kernel, which says nothing about the actual cause.
"""
check_grid(scheme, Nx, Ny) = nothing

function check_grid(scheme::AbstractParameterization, Nx, Ny)
    for name in fieldnames(typeof(scheme))
        value = getfield(scheme, name)
        value isa AbstractArray || continue
        size(value) == (Nx, Ny) || throw(
            DimensionMismatch(
                "$(nameof(typeof(scheme))) field `$name` is $(size(value)) but the model grid is ($Nx, $Ny)"
            )
        )
    end
    return nothing
end

"""
    instantiate(scheme, grid)

Return a physics parameterization ready for the model: a scheme **type** is default-constructed at
the grid's precision (`Scheme{eltype(grid)}(grid)`), while a ready-made **instance** is returned
unchanged. Non-default schemes are constructed by the caller, e.g.
`PrognosticAlbedo{Float32}(grid; adm = 200, adc = my_array)`.
"""
instantiate(scheme::Type, grid) = scheme{eltype(grid)}(grid)
instantiate(scheme, grid) = scheme

"""
    reconstruct(x; kwargs...)

Copy the immutable struct `x` with the named fields replaced. `Parameters` is rebuilt
rather than mutated so that it stays isbits and can cross into a kernel by value.
"""
function reconstruct(x::T; kwargs...) where {T}
    names = fieldnames(T)
    fields = NamedTuple{names}(map(f -> getfield(x, f), names))
    return T(; merge(fields, NamedTuple(kwargs))...)
end

# Route a scalar into the immutable Parameters via a functional update.
@inline function set_param!(fsm, sym::Symbol, value)
    v = convert(fieldtype(typeof(fsm.params), sym), value)
    fsm.params = reconstruct(fsm.params; NamedTuple{(sym,)}((v,))...)
    return fsm
end

"""
    apply_params!(fsm, params)

Apply parameter overrides: a `Parameters` scalar is reconstructed onto `fsm.params`;
a `Surface` per-cell array is filled (scalar) or copied (array) in place.
"""
function apply_params!(fsm, params)
    for (key, value) in params
        sym = Symbol(key)
        if hasfield(typeof(fsm.params), sym)
            set_param!(fsm, sym, value)
        elseif hasfield(typeof(fsm.surface), sym)
            arr = getfield(fsm.surface, sym)
            value isa AbstractArray ? (arr .= eltype(arr).(value)) : fill!(arr, eltype(arr)(value))
        else
            throw(ArgumentError("unknown parameter override \"$key\""))
        end
    end
    return fsm
end
