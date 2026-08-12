"""
    setup(Tf, Ti, landuse, Nx, Ny, settings)

Initialize the FSM snow model with specified configuration and domain properties.

Creates and configures the FSM model state structure with landuse data, model parameters,
and domain dimensions. Sets up soil properties, surface characteristics, and applies 
configuration-specific settings for different surface types and model behaviors.

# Arguments
- `Tf`: Floating-point precision type (typically Float32 or Float64)
- `Ti`: Integer type for array indices (typically Int32 or Int64)
- `landuse::Dict`: Landuse data dictionary with topographic and surface properties
- `Nx::Int, Ny::Int`: Model domain dimensions
- `settings::Dict`: Configuration dictionary containing:
  - `"tile"`: Surface tile type ("open", "forest", "glacier")
  - `"config"` (optional): Model configuration flags (SNFRAC, CANMOD, EXCHNG, ZOFFST, etc.)
  - `"params"` (optional): Parameter overrides (hfsn, z0sn, etc.)

# Returns
- `FSM`: Initialized model state structure ready for simulation
"""
function setup(Tf, Ti, landuse::Dict, Nx::Int, Ny::Int, settings::Dict)

    @unpack_constants(Tf)

    # Create fsm object
    fsm = FSM{Tf, Ti}(Nx = Nx, Ny = Ny)

    # Set tile
    fsm.TILE = settings["tile"]

    # Apply model configuration
    if haskey(settings, "config")
        for (key, value) in settings["config"]
            if value isa Int
                value = Ti(value)
            end
            setfield!(fsm, Symbol(key), value)
        end
    end

    # Apply parameter overrides
    if haskey(settings, "params")
        for (key, value) in settings["params"]
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
    end

    # this should be set in the settings parameters like settings["params"]["rhof"] = 300.0
    # and not happen automatically if FSNRHO=0
    # Settings specific for FSNRHO=0 (fixed fresh snow density)
    # if (fsm.FSNRHO == 0)
    #     fsm.rhof = fsm.rho0
    # end

    # Derived soil parameters
    # should we move all these soil params into types.jl or is it fine if they are hard coded?

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
    # if (fsm.TILE != "forest")
    #     # these moved to types.jl as default parameters dictionary
    #     fsm.VAI[:, :] .= Tf(0)
    #     fsm.hcan[:, :] .= Tf(0)
    #     fsm.fsky[:, :] .= Tf(1)
    #     fsm.trcn[:, :] .= exp.(-fsm.kdif .* fsm.VAI[:, :])
    #     fsm.fveg[:, :] .= Tf(1) .- exp.(-fsm.kveg .* fsm.VAI[:, :])
    #     fsm.fves[:, :] .= Tf(1) .- exp.(-fsm.kveg .* fsm.VAI[:, :])
    # else

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

    # Tuned snow surface properties (same for glacier, but there given in settings["params"]["z0_snow"])
    if (fsm.SNTRAN == 1 || fsm.SNSLID == 1)
        fsm.z0_snow[fsm.glacierfrac .> eps(Tf)] .= fsm.z0gl
    end

    # Initialize SnowSlide arrays if enabled    TODO this is computed using Float64 in matlab originally...
    if fsm.SNSLID == 1
        # Allocate SnowSlide arrays
        fsm.slope = zeros(Tf, Nx, Ny)
        fsm.Shd = zeros(Tf, Nx, Ny)
        fsm.dSWE_tot_slide = zeros(Tf, Nx, Ny)
        fsm.index_sorted_dem = zeros(Ti, Nx * Ny, 2)

        # Sort DEM indices for processing order (highest to lowest elevation)
        sort_dem_indices!(fsm.index_sorted_dem, fsm.dem)

        # Load slope data
        if haskey(landuse, "slope")
            fsm.slope = Tf.(landuse["slope"]["data"])
        else
            error("SNSLID=1 but no 'slope' data found in landuse dictionary.")     # TODO maybe remove, should always be there...
        end

        # Calculate snow holding depth from slope
        slope_thres = copy(fsm.slope)
        slope_thres[slope_thres .< fsm.snow_slide_slope_floor] .= fsm.snow_slide_slope_floor  # Limit to >10 degrees to avoid inf

        # Snow holding depth, normal to the slope
        shd_norm = fsm.snow_slide_shd_a .* slope_thres .^ (fsm.snow_slide_shd_b)

        # Convert to vertical snow holding depth
        cos_slope_thres = cosd.(slope_thres)
        cos_slope_thres[cos_slope_thres .< fsm.snow_slide_cos_floor] .= fsm.snow_slide_cos_floor  # Avoid division by zero
        fsm.Shd = shd_norm .* cos_slope_thres

    end

    return fsm

end


"""
    sort_dem_indices!(index_sorted_dem, dem)

Sort DEM indices from highest to lowest elevation for SnowSlide processing order.

This function sorts a DEM of size N x Ny by decreasing elevation.
It returns a (Nx*Ny,2) array: dimension 1 is sorted by decreasing
elevation; dimension 2 contains row and column of the given pixels.

# Arguments
- `index_sorted_dem::Matrix{Ti}`: Output array for sorted (i,j) indices
- `dem::Matrix{Tf}`: Digital elevation model
"""
function sort_dem_indices!(index_sorted_dem::Matrix{Ti}, dem::Matrix{Tf}) where {Tf <: Real, Ti <: Integer}
    # Get sorted indices by decreasing elevation
    ind = sortperm(vec(dem), rev = true)

    # Get dimensions
    rows, cols = size(dem)

    # Create 2D index arrays
    rows2D = repeat((1:rows)', outer = (1, cols))
    cols2D = repeat((1:cols)', outer = (rows, 1))

    # Extract row and column indices for sorted positions
    ind_rows = rows2D[ind]
    ind_cols = cols2D[ind]

    # Fill output array
    index_sorted_dem[:, 1] = Ti.(ind_rows)
    index_sorted_dem[:, 2] = Ti.(ind_cols)

    return nothing
end
