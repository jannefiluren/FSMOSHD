# Snow-transport workspace.
#
# Transport (SnowSlide + SnowTran3D) is a neighbour-coupled, grid-global operator, run as a
# standalone step rather than a per-cell `step!` stage. ALL transport-specific state lives here
# instead of on the model, keeping the decomposed, GPU-adaptable FSM untouched: the static
# per-cell setup arrays, the cumulative-change accumulators, the tuning constants (defaults match
# deps/MODULES.F90 so the Julia and Fortran paths agree), and the SnowTran3D working arrays.
# CPU-only (concrete `Matrix`): transport does neighbour-coupled scalar indexing, not ported to
# the GPU. Build one with `setup_transport(fsm, landuse)`.

@kwdef mutable struct SnowTransport{Tf}

    # Domain size
    Nx::Int = 1                                               # Size of first array dimension (rows)
    Ny::Int = 1                                               # Size of second array dimension (columns)

    # Which processes run, and which implementation, when transport! is called from step!
    wind::Bool = false                                      # Run wind transport (SnowTran3D)
    slide::Bool = false                                     # Run snow slides (SnowSlide)
    use_fortran::Bool = false                               # Use the Fortran ccall path (else the Julia port)

    # Tuning constants (defaults match deps/MODULES.F90)
    rhos_min::Tf = 50                                        # Minimum snow density (kg/m^3)
    rhos_max::Tf = 750                                       # Maximum snow density (kg/m^3)
    rho_snow::Tf = 300                                       # Constant snow density for transport (kg/m^3)
    dyn_ratio::Tf = 0.7                                      # Dynamic snow holding depth ratio (-)
    trig_ratio::Tf = 1.5                                     # Hysteretic triggering snow holding depth ratio (-)
    rho_deposit::Tf = 300                                    # Constant snow avalanche deposit density (kg/m^3)
    slope_min::Tf = 25                                       # Minimum slope for snow slide occurrence (deg)
    Shd_min::Tf = 0.05                                       # Minimum snow holding depth (m)
    tiled_trans_run::Bool = false                            # Tiled transport run flag

    # Static per-cell setup arrays (built in setup_transport)
    slope::Matrix{Tf} = fill(Tf(NaN), Nx, Ny)               # Slope angles (deg)
    Shd::Matrix{Tf} = fill(Tf(NaN), Nx, Ny)                 # Snow holding depth for gravitational transport (m)
    vegsnowd_xy::Matrix{Tf} = fill(Tf(0.1), Nx, Ny)         # Vegetation snow holding capacity (m)
    forestfrac::Matrix{Tf} = zeros(Tf, Nx, Ny)              # Forest fraction (-)
    index_sorted_dem::Matrix{Int} = zeros(Int, Nx * Ny, 2)    # DEM indices sorted by decreasing elevation

    # Cumulative SWE-change accumulators (kg/m^2)
    dSWE_tot_slide::Matrix{Tf} = zeros(Tf, Nx, Ny)          # ... due to slides
    dSWE_tot_salt::Matrix{Tf} = zeros(Tf, Nx, Ny)           # ... due to saltation
    dSWE_tot_susp::Matrix{Tf} = zeros(Tf, Nx, Ny)           # ... due to suspension
    dSWE_tot_subl::Matrix{Tf} = zeros(Tf, Nx, Ny)           # ... due to sublimation

    # Per-step SWE-change outputs (scratch for transport!, one process each)
    dSWE_salt::Matrix{Tf} = zeros(Tf, Nx, Ny)               # This step's saltation SWE change
    dSWE_susp::Matrix{Tf} = zeros(Tf, Nx, Ny)               # This step's suspension SWE change
    dSWE_subl::Matrix{Tf} = zeros(Tf, Nx, Ny)               # This step's sublimation SWE change
    dSWE_slide::Matrix{Tf} = zeros(Tf, Nx, Ny)              # This step's slide SWE change

    # SnowSlide work arrays
    snow_depo::Matrix{Bool} = zeros(Bool, Nx, Ny)           # Pixels receiving slide deposits
    Shd_corr::Matrix{Tf} = zeros(Tf, Nx, Ny)                # Thresholded snow holding depth (m)

    # SnowTran3D work arrays for snowtran3d_julia!
    uwind::Matrix{Tf} = zeros(Tf, Nx, Ny)                   # x component of wind speed (m/s)
    vwind::Matrix{Tf} = zeros(Tf, Nx, Ny)                   # y component of wind speed (m/s)
    snowthickness::Matrix{Tf} = zeros(Tf, Nx, Ny)           # Snowpack depth, not scaled by fsnow (m)
    veg_z0::Matrix{Tf} = zeros(Tf, Nx, Ny)                  # Vegetation roughness length (m)
    Ds_soft::Matrix{Tf} = zeros(Tf, Nx, Ny)                 # Soft snow thickness (m)
    Utau::Matrix{Tf} = zeros(Tf, Nx, Ny)                    # Friction velocity (m/s)
    Utau_t::Matrix{Tf} = zeros(Tf, Nx, Ny)                  # Threshold friction velocity (m/s)
    h_star::Matrix{Tf} = zeros(Tf, Nx, Ny)                  # Height of the saltation layer (m)
    z_0::Matrix{Tf} = zeros(Tf, Nx, Ny)                     # Surface roughness length (m)
    Qsalt::Matrix{Tf} = zeros(Tf, Nx, Ny)                   # Saltation flux (kg/m/s)
    Qsalt_u::Matrix{Tf} = zeros(Tf, Nx, Ny)                 # x component of saltation flux (kg/m/s)
    Qsalt_v::Matrix{Tf} = zeros(Tf, Nx, Ny)                 # y component of saltation flux (kg/m/s)
    Qsalt_max::Matrix{Tf} = zeros(Tf, Nx, Ny)              # Maximum possible saltation flux (kg/m/s)
    Qsalt_maxu::Matrix{Tf} = zeros(Tf, Nx, Ny)             # x component of Qsalt_max (kg/m/s)
    Qsalt_maxv::Matrix{Tf} = zeros(Tf, Nx, Ny)             # y component of Qsalt_max (kg/m/s)
    conc_salt::Matrix{Tf} = zeros(Tf, Nx, Ny)              # Saltation reference-level mass concentration (kg/m^3)
    Qsusp::Matrix{Tf} = zeros(Tf, Nx, Ny)                  # Suspension flux (kg/m/s)
    Qsusp_u::Matrix{Tf} = zeros(Tf, Nx, Ny)               # x component of suspension flux (kg/m/s)
    Qsusp_v::Matrix{Tf} = zeros(Tf, Nx, Ny)               # y component of suspension flux (kg/m/s)
    Qsubl::Matrix{Tf} = zeros(Tf, Nx, Ny)                 # Sublimation flux (kg/m^2/s)

    # SnowTran3D work arrays for getnewdepth!
    dh_s_u::Matrix{Tf} = zeros(Tf, Nx, Ny)                 # x component of snow depth change (m)
    dh_s_v::Matrix{Tf} = zeros(Tf, Nx, Ny)                 # y component of snow depth change (m)
    dSWE_s_u::Matrix{Tf} = zeros(Tf, Nx, Ny)              # x component of SWE change (kg/m^2)
    dSWE_s_v::Matrix{Tf} = zeros(Tf, Nx, Ny)              # y component of SWE change (kg/m^2)
    dSWE_s_u_loss::Matrix{Tf} = zeros(Tf, Nx, Ny)         # x component of SWE loss (kg/m^2)
    dSWE_s_v_loss::Matrix{Tf} = zeros(Tf, Nx, Ny)         # y component of SWE loss (kg/m^2)
    dh_s_u_loss::Matrix{Tf} = zeros(Tf, Nx, Ny)           # x component of snow depth loss (m)
    dh_s_v_loss::Matrix{Tf} = zeros(Tf, Nx, Ny)           # y component of snow depth loss (m)
    dSWE_s_u_gain::Matrix{Tf} = zeros(Tf, Nx, Ny)         # x component of SWE gain (kg/m^2)
    dSWE_s_v_gain::Matrix{Tf} = zeros(Tf, Nx, Ny)         # y component of SWE gain (kg/m^2)

    # Wind direction index arrays
    index_ue::Matrix{Int} = zeros(Int, Nx, 2 * Ny + 1)        # Wind index array E
    index_uw::Matrix{Int} = zeros(Int, Nx, 2 * Ny + 1)        # Wind index array W
    index_vn::Matrix{Int} = zeros(Int, Ny, 2 * Nx + 1)        # Wind index array N
    index_vs::Matrix{Int} = zeros(Int, Ny, 2 * Nx + 1)        # Wind index array S

end
