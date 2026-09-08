# One-time setup for the snow-transport operators: build a SnowTransport workspace from the
# model grid + a landuse dictionary. Mirrors the transport setup block that used to live in
# setup.jl on origin/jan-cpu_gpu_refactor, but writes into the standalone workspace rather than
# onto FSM. Both the Julia and Fortran transport paths consume this same workspace, so any
# array computed here is shared and cannot make the two paths disagree.

"""
    setup_transport(fsm, landuse; wind = false, slide = false, use_fortran = false, tiled = false) -> SnowTransport

Allocate and initialise a [`SnowTransport`](@ref) workspace for `fsm`:

- loads `slope` (deg) and `forest` fraction from the `landuse` dictionary,
- derives the snow holding depth `Shd` from slope,
- sorts the DEM indices (highest → lowest) for the SnowSlide processing order,
- sets the vegetation snow-holding capacity and the tuning constants (defaults match
  `deps/MODULES.F90`).

`wind`/`slide` select which processes [`transport!`](@ref) runs; `use_fortran` chooses the
Fortran ccall path over the Julia port; `tiled` sets the tiled-transport flag. Transport is
CPU-only, so `fsm` must hold host arrays.
"""
function setup_transport(
        fsm::FSM{Tf, Ti}, landuse::Dict;
        wind::Bool = false, slide::Bool = false, use_fortran::Bool = false, tiled::Bool = false,
    ) where {Tf, Ti}

    Nx = Int(fsm.grid.Nx)
    Ny = Int(fsm.grid.Ny)

    w = SnowTransport{Tf, Ti}(Nx = Ti(Nx), Ny = Ti(Ny))
    w.wind = wind
    w.slide = slide
    w.use_fortran = use_fortran
    w.tiled_trans_run = tiled

    # Vegetation snow-holding capacity (constant on the operational setup)
    w.vegsnowd_xy .= Tf(0.1)

    # Forest fraction (optional; defaults to zero)
    if haskey(landuse, "forest")
        w.forestfrac .= Tf.(landuse["forest"]["data"])
    end

    # Slope (deg) — required
    haskey(landuse, "slope") ||
        throw(ArgumentError("setup_transport requires 'slope' data in the landuse dictionary"))
    w.slope .= Tf.(landuse["slope"]["data"])

    # Snow holding depth from slope: normal-to-slope value projected to the vertical.
    # Slope floored at 10 deg to avoid the divergence of the power law near 0.
    slope_thres = max.(w.slope, Tf(10))
    shd_norm = Tf(3178.4) .* slope_thres .^ Tf(-1.998)
    w.Shd .= shd_norm .* max.(cosd.(slope_thres), Tf(0.001))

    # Processing order for SnowSlide: cells from highest to lowest elevation
    sort_dem_indices!(w.index_sorted_dem, fsm.surface.dem)

    return w
end

"""
    sort_dem_indices!(index_sorted_dem, dem)

Sort DEM indices from highest to lowest elevation for the SnowSlide processing order. Fills
`index_sorted_dem`, an `(Nx*Ny, 2)` array whose rows are `(i, j)` pairs ordered by decreasing
`dem`.
"""
function sort_dem_indices!(index_sorted_dem::Matrix{Ti}, dem::AbstractMatrix{Tf}) where {Tf <: Real, Ti <: Integer}
    ind = sortperm(vec(dem), rev = true)

    rows, cols = size(dem)
    rows2D = repeat((1:rows), outer = (1, cols))
    cols2D = repeat((1:cols)', outer = (rows, 1))

    index_sorted_dem[:, 1] = Ti.(rows2D[ind])
    index_sorted_dem[:, 2] = Ti.(cols2D[ind])

    return nothing
end
