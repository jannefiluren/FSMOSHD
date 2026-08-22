@with_kw mutable struct FSM{Tf, Ti}

    # Layer configuration

    Dzsnow::Vector{Tf} = [0.1, 0.2, 0.4]                     # Maximum snow layer thicknesses (m)
    Dzsoil::Vector{Tf} = [0.1, 0.2, 0.4, 0.8]                # Maximum soil layer thicknesses (m)
    Nsmax::Ti = length(Dzsnow)                               # Number of snow layers
    Nsoil::Ti = length(Dzsoil)                               # Number of soil layers

    # Domain size

    Nx::Ti = 1                                               # Size of first array dimension (rows)
    Ny::Ti = 1                                               # Size of second array dimension (columns)

    # Driving data

    dt::Tf = 3600                                            # Time step (s)
    zT::Tf = 10                                              # Temperature measurement height (m)
    zU::Tf = 10                                              # Wind speed measurement height (m)
    zRH::Tf = 10                                             # Relative humidity measurement height (m)

    # Model configuration

    ALBEDO::Ti = 2                                           # Snow albedo (0, 1, 2)
    CANMOD::Ti = 0                                           # Forest canopy (0, 1)
    CONDCT::Ti = 1                                           # Snow thermal conductivity (0, 1)
    DENSTY::Ti = 3                                           # Snow density (0, 1, 2, 3)
    EXCHNG::Ti = 1                                           # Turbulent exchange (0, 1)
    HYDROL::Ti = 2                                           # Snow hydraulics (0, 1, 2)
    SNFRAC::Ti = 3                                           # Snow cover fraction (0, 1, 2, 3, 4)
    ZOFFST::Ti = 0                                           # Measurement height offset (0, 1)
    FSNRHO::Ti = 2                                           # Fresh snow density (0, 1, 2)
    ALRADT::Ti = 1                                           # Albedo decay as function of incoming direct shortwave radiation (0, 1)
    SNTRAN::Ti = 0                                           # Wind-driven snow transport (0, 1)
    SNSLID::Ti = 0                                           # Snow slides (0, 1)
    SNOLAY::Ti = 0                                           # Density-dependent layering (0, 1)
    HN_ON::Bool = false                                      # Activate new snow model
    Z0PERT::Bool = false                                     # Activate snow roughness length perturbations
    WCPERT::Bool = false                                     # Activate liquid water capacity perturbations
    FSPERT::Bool = false                                     # Activate fresh snow density perturbations
    ALPERT::Bool = false                                     # Activate albedo perturbations
    SLPERT::Bool = false                                     # Activate settling perturbations

    # Tile options

    TILE::String = "open"                                    # Tile type
    tthresh::Tf = 0.1                                        # Tile threshold

    # Numerical solution parameters

    Nitr = 4                                                 # Number of iterations for surface energy balance

    # Canopy parameters

    avg0::Tf = 0.1                                           # Snow-free vegetation albedo (-)
    avgs::Tf = 0.4                                           # Snow-covered vegetation albedo (-)
    cden::Tf = 0.004                                         # Dense canopy turbulent transfer coefficient (-)
    cvai::Tf = 4.4                                           # Canopy snow capacity per unit vegetation area index (kg/m^2)
    cveg::Tf = 20                                            # Vegetation turbulent transfer coefficient ((s/m)^0.5)
    Gcn1::Tf = 0.5                                           # Leaf angle distribution parameter (-)
    Gcn2::Tf = 0                                             # Leaf angle distribution parameter (-)
    gsnf::Tf = 0                                             # Snow-free vegetation moisture conductance (m/s)
    kdif::Tf = 0.5                                           # Diffuse radiation extinction coefficient (-)
    kveg::Tf = 1                                             # Canopy cover coefficient (-)
    rchd::Tf = 0.67                                          # Ratio of displacement height to canopy height (-)
    rchz::Tf = 0.2                                           # Ratio of roughness length to canopy height (-)
    tcnc::Tf = 3600 * 240                                    # Canopy unloading time scale for cold snow (s)
    tcnm::Tf = 3600 * 48                                     # Canopy unloading time scale for melting snow (s)
    pmultf_for::Tf = 0.5                                     # Multiplier for snowfall in forest (-)

    # Snow parameters

    a_eta::Tf = 0.1                                          # Temperature factor for Crocus B92 compaction (K^-1)
    asmx::Tf = 0.86                                          # Maximum albedo for fresh snow (-)
    asmn::Tf = 0.6                                           # Minimum albedo for melting snow (-)
    b_eta::Tf = 0.023                                        # First density factor for Crocus B92 compaction (m^3/kg)
    bthr::Tf = 2                                             # Snow thermal conductivity exponent (-)
    c_eta::Tf = 250                                          # Second density factor for Crocus B92 compaction (kg/m^3)
    eta0::Tf = 3.7e7                                         # Reference snow viscosity (Pa s)
    eta1::Tf = 7.62237e6                                     # Reference snow viscosity for Crocus B92 compaction (Pa s)
    hfsn::Tf = 0.1                                           # Snowcover fraction depth scale (m)
    kfix::Tf = 0.24                                          # Fixed thermal conductivity of snow (W/m/K)
    rho0::Tf = 300                                           # Fixed snow density (kg/m^3)
    rhob::Tf = 6                                             # Temperature factor in fresh snow density (kg/m^3/K)
    rhoc::Tf = 26                                            # Wind factor in fresh snow density (kg s^0.5/m^3.5)
    rhof::Tf = 109                                           # Fresh snow density (kg/m^3)
    rhos_min::Tf = 50                                        # Minimum snow density (kg/m^3)
    rhos_max::Tf = 750                                       # Maximum snow density (kg/m^3)
    rcld::Tf = 300                                           # Maximum density for cold snow (kg/m^3)
    rgr0::Tf = 5.0e-5                                        # Fresh snow grain radius (m)
    rmlt::Tf = 500                                           # Maximum density for melting snow (kg/m^3)
    Salb::Tf = 10                                            # Albedo decay constant (kg/m^2)
    snda::Tf = 2.8e-6                                        # Thermal metamorphism parameter (1/s)
    Talb::Tf = -2                                            # Albedo decay temperature threshold (C)
    tcld::Tf = 3600 * 1000                                   # Cold snow albedo decay time scale (s)
    tmlt::Tf = 3600 * 100                                    # Melting snow albedo decay time scale (s)
    trho::Tf = 3600 * 200                                    # Snow compaction time scale (s)
    Wirr::Tf = 0.03                                          # Irreducible liquid water content of snow (-)
    z0gl::Tf = 0.0009                                        # Roughness length of snow on glaciers (m)
    Sfmin::Tf = 10                                           # Minimum snowfall over 24h needed to refresh albedo (kg/m^2)

    # Snow layering parameters

    Ds_min::Tf = 0.01                                        # Minimum possible snow layer thickness (m)
    Ds_surflay::Tf = 0.5                                     # Maximum thickness of surface fine snow layering (m)

    # Snow transport parameters

    dyn_ratio::Tf = 0.09                                     # Dynamic snow holding depth ratio (-)
    rho_deposit::Tf = 300.0                                  # Constant snow avalanche deposit density (kg/m³)
    slope_min::Tf = 30.0                                     # Minimum slope for snow slide occurrence (deg)
    Shd_min::Tf = 0.01                                       # Minimum snow holding depth (m)
    rho_snow::Tf = 300.0                                     # Constant snow density for transport (kg/m³)
    snow_slide_slope_floor::Tf = 10.0                        # Slope threshold for calculating snow holding depth normal to slope (deg)
    snow_slide_shd_a::Tf = 3178.4                            # Parameter for computing snow holding depth (-)
    snow_slide_shd_b::Tf = -1.998                            # Exponent for computing snow holding depth (-)
    snow_slide_cos_floor::Tf = 0.001                         # Cosine threshold for calculating vertical snow holding depth (-)

    # Ground surface parameters

    bstb::Tf = 5                                             # Atmospheric stability parameter (-)
    gsat::Tf = 0.01                                          # Surface conductance for saturated soil (m/s)

    # Additional forest snow process parameters

    adfs::Tf = 3                                             # Snow albedo adjustment dependent on shortwave radiation (-)
    adfl::Tf = 2                                             # Snow albedo adjustment dependent on longwave radiation (-)
    fsar::Tf = 0.1                                           # Snow albedo adjustment range dependent on vegetation fraction (-)
    psf::Tf = 1                                              # Solid precipitation multiplier in forest at minimum canopy cover (-)
    psr::Tf = 0.1                                            # Additional multiplier range across canopy cover (-)
    wcan::Tf = 2.5                                           # Parameter of exponential wind profile (-)
    zsub::Tf = 2                                             # Sub-canopy reference height (m)
    zgf::Tf = 1                                              # Roughness length adjustment factor depending on vegetation fraction (-)
    zgr::Tf = 0                                              # Roughness length adjustment range depending on vegetation fraction (-)
    khcf::Tf = 3                                             # Diffusivity adjustment for canopy effects (-)

    # Surface parameters

    adm::Tf = 100                                            # Melting snow albedo decay time (h)
    adc::Array{Tf, 2} = Tf(1000) * ones(Nx, Ny)              # Cold snow albedo decay time (h)
    afs::Array{Tf, 2} = asmx * ones(Nx, Ny)                  # Maximum albedo for fresh snow
    z0_snow::Array{Tf, 2} = 0.002 * ones(Nx, Ny)             # Roughness length of snow (m)

    # Surface properties

    alb0::Array{Tf, 2} = 0.2 * ones(Nx, Ny)                  # Snow-free ground albedo (-)
    z0sf::Array{Tf, 2} = 0.2 * ones(Nx, Ny)                  # Snow-free roughness length (m)
    fcly::Array{Tf, 2} = 0.3 * ones(Nx, Ny)                  # Soil clay fraction (-)
    fsnd::Array{Tf, 2} = 0.6 * ones(Nx, Ny)                  # Soil sand fraction (-)
    fsat::Tf = 0.5                                           # Initial moisture content of soil layers as fractions of saturation
    Tprof::Tf = 285                                          # Initial soil layer temperatures (K)

    # Canopy parameters (dummy values should be filled from landuse data)

    VAI::Array{Tf, 2} = zeros(Nx, Ny)                        # Vegetation area index (-)
    vfhp::Array{Tf, 2} = fill(NaN, Nx, Ny)                   # Hemispherical sky-view fraction including canopy (-)
    canh::Array{Tf, 2} = fill(NaN, Nx, Ny)                   # Canopy heat capacity (J/K/m^2)
    fsky::Array{Tf, 2} = ones(Nx, Ny)                        # Sky view fraction (-)
    fveg::Array{Tf, 2} = Tf(1) .- exp.(-kveg .* VAI[:, :])   # Canopy cover fraction (-)
    fves::Array{Tf, 2} = Tf(1) .- exp.(-kveg .* VAI[:, :])   # Stand-scale canopy cover fraction (-)
    hcan::Array{Tf, 2} = zeros(Nx, Ny)                       # Canopy height (m)
    lai::Array{Tf, 2} = fill(NaN, Nx, Ny)                    # Leaf area index (-)
    pmultf::Array{Tf, 2} = fill(NaN, Nx, Ny)                 # Precipitation multiplier to revert correction applied to open area (-)
    scap::Array{Tf, 2} = fill(NaN, Nx, Ny)                   # Canopy snow capacity (kg/m^2)
    trcn::Array{Tf, 2} = exp.(-kdif .* VAI[:, :])            # Canopy transmissivity (-)


    # Terrain properties

    slopemu::Array{Tf, 2} = fill(NaN, Nx, Ny)                # Slope parameter (-)
    xi::Array{Tf, 2} = fill(NaN, Nx, Ny)                     # Terrain correlation length (m)
    Ld::Array{Tf, 2} = fill(NaN, Nx, Ny)                     # Grid cell size (m)
    fsky_terr::Array{Tf, 2} = fill(NaN, Nx, Ny)              # Sky view fraction terrain (-)
    dem::Array{Tf, 2} = fill(NaN, Nx, Ny)                    # Grid elevation (m)
    tilefrac::Array{Tf, 2} = ones(Nx, Ny)                    # Tile fraction (-)
    glacierfrac::Array{Tf, 2} = fill(NaN, Nx, Ny)            # Glacier fraction (-)
    vegsnowd_xy::Array{Tf, 2} = Tf(0.1) * ones(Nx, Ny)       # Vegetation snow holding capacity (m)
    prec_multi::Array{Float64, 2} = fill(NaN, Nx, Ny)        # Precipitation multiplier (-)    TODO use float64 to match matlab/fortran version - change precision later
    slope::Array{Tf, 2} = fill(NaN, Nx, Ny)                  # Slope angles (deg)
    Shd::Array{Tf, 2} = fill(NaN, Nx, Ny)                    # Snow holding depth for gravitational transport (m)

    # Derived soil parameters

    b::Array{Tf, 2} = zeros(Nx, Ny)                          # Clapp-Hornberger exponent (-)
    hcap_soil::Array{Tf, 2} = zeros(Nx, Ny)                  # Volumetric heat capacity of dry soil (J/K/m^3)
    hcon_soil::Array{Tf, 2} = zeros(Nx, Ny)                  # Thermal conductivity of dry soil (W/m/K)
    sathh::Array{Tf, 2} = zeros(Nx, Ny)                      # Saturated soil water pressure (m)
    Vsat::Array{Tf, 2} = zeros(Nx, Ny)                       # Volumetric soil moisture at saturation (-)
    Vcrit::Array{Tf, 2} = zeros(Nx, Ny)                      # Volumetric soil moisture at critical point (-)

    # State variables

    albs::Array{Tf, 2} = Tf(0.85) * ones(Nx, Ny)             # Snow albedo (-)
    Ds::Array{Tf, 3} = zeros(Nsmax, Nx, Ny)                  # Snow layer thicknesses (m)
    Nsnow::Array{Ti, 2} = zeros(Ti, Nx, Ny)                  # Number of snow layers
    Qcan::Array{Tf, 2} = zeros(Nx, Ny)                       # Canopy air space humidity (kg/kg)
    rgrn::Array{Tf, 3} = zeros(Nsmax, Nx, Ny)                # Snow layer grain radius (m)
    Sice::Array{Tf, 3} = zeros(Nsmax, Nx, Ny)                # Ice content of snow layers (kg/m^2)
    Sliq::Array{Tf, 3} = zeros(Nsmax, Nx, Ny)                # Liquid content of snow layers (kg/m^2)
    Sveg::Array{Tf, 2} = zeros(Nx, Ny)                       # Snow mass on vegetation (kg/m^2)
    Tcan::Array{Tf, 2} = Tf(285) * ones(Nx, Ny)              # Canopy air space temperature (K)
    theta::Array{Tf, 3} = zeros(Nsoil, Nx, Ny)               # Volumetric moisture content of soil layers (-)
    Tsnow::Array{Tf, 3} = Tf(273.15) * ones(Nsmax, Nx, Ny)   # Snow layer temperatures (K)
    Tsoil::Array{Tf, 3} = Tf(285) * ones(Nsoil, Nx, Ny)      # Soil layer temperatures (K)
    Tsrf::Array{Tf, 2} = Tf(285) * ones(Nx, Ny)              # Surface skin temperature (K)
    fsnow::Array{Tf, 2} = zeros(Nx, Ny)                      # Snow cover fraction (-)
    Tveg::Array{Tf, 2} = Tf(285) * ones(Nx, Ny)              # Vegetation temperature (K)
    snowdepthmin::Array{Tf, 2} = zeros(Nx, Ny)               # Minimum snow depth at time step of swemin (m)
    snowdepthmax::Array{Tf, 2} = zeros(Nx, Ny)               # Maximum snow depth at time step of swemax (m)
    snowdepthhist::Array{Tf, 3} = zeros(14, Nx, Ny)          # History of snow depth during last 14 days with most recent entries first (m)
    swemin::Array{Tf, 2} = zeros(Nx, Ny)                     # Minimum SWE during the season (kg/m^2)
    swemax::Array{Tf, 2} = zeros(Nx, Ny)                     # Maximum SWE during the season (kg/m^2)
    swehist::Array{Tf, 3} = zeros(14, Nx, Ny)                # History of SWE during last 14 days with most recent entries first (kg/m^2)
    histowet::Array{Tf, 3} = zeros(Nsmax, Nx, Ny)            # Historical variable for past wetting of a layer (-)

    # Variables used in drive-function

    es::Array{Tf, 2} = zeros(Nx, Ny)                         # Saturation vapour pressure (Pa)
    Qa::Array{Tf, 2} = zeros(Nx, Ny)                         # Specific humidity (kg/kg)
    Uaeff::Array{Tf, 2} = zeros(Nx, Ny)                     # Wind speed with lower bound applied (m/s)
    Sfeff::Array{Tf, 2} = zeros(Nx, Ny)                     # Snowfall rate reaching the surface, adjusted by canopy processes (kg/m^2/s)

    # Variables used in radiation-function

    alb::Array{Tf, 2} = zeros(Nx, Ny)                        # Albedo (-)
    asrf_out::Array{Tf, 2} = zeros(Nx, Ny)                   # Surface albedo (-)
    SWveg::Array{Tf, 2} = zeros(Nx, Ny)                      # Net short wave radiation absorbed by vegetation (W/m^2)
    SWsrf::Array{Tf, 2} = zeros(Nx, Ny)                      # Net short wave radiation absorbed by the surface (W/m^2)
    SWsci::Array{Tf, 2} = zeros(Nx, Ny)                      # Subcanopy incoming shortwave radiation (W/m^2)
    LWeff::Array{Tf, 2} = zeros(Nx, Ny)                      # Incoming longwave radiation used in the energy balance, terrain-corrected where applicable (W/m^2)

    # Variables used in thermal-function

    ksnow::Array{Tf, 3} = zeros(Nsmax, Nx, Ny)               # Thermal conductivity of snow (W/m/K)
    csoil::Array{Tf, 3} = zeros(Nsoil, Nx, Ny)               # Areal heat capacity of soil (J/K/m^2)
    ksoil::Array{Tf, 3} = zeros(Nsoil, Nx, Ny)               # Thermal conductivity of soil (W/m/K)
    gs1::Array{Tf, 2} = zeros(Nx, Ny)                        # Surface moisture conductance (m/s)
    Ds1::Array{Tf, 2} = zeros(Nx, Ny)                        # Surface layer thickness (m)
    Ts1::Array{Tf, 2} = zeros(Nx, Ny)                        # Surface layer temperature (K)
    ks1::Array{Tf, 2} = zeros(Nx, Ny)                        # Surface thermal conductivity (W/m/K)
    Tveg0::Array{Tf, 2} = zeros(Nx, Ny)                      # Vegetation temperature at start of timestep (K)

    # Variables used in sfexch-function

    KH::Array{Tf, 2} = zeros(Nx, Ny)                         # Eddy diffusivity for heat to the atmosphere (m/s)
    KHa::Array{Tf, 2} = zeros(Nx, Ny)                        # Eddy diffusivity from the canopy air space (m/s)
    KHg::Array{Tf, 2} = zeros(Nx, Ny)                        # Eddy diffusivity for heat from the ground (m/s)
    KHv::Array{Tf, 2} = zeros(Nx, Ny)                        # Eddy diffusivity for heat from vegetation (m/s)
    KWg::Array{Tf, 2} = zeros(Nx, Ny)                        # Eddy diffusivity for water from the ground (m/s)
    KWv::Array{Tf, 2} = zeros(Nx, Ny)                        # Eddy diffusivity for water from vegetation (m/s)
    Usc::Array{Tf, 2} = zeros(Nx, Ny)                        # Wind speed in canopy layer (m/s)

    # Variables used in ebalsrf-function

    Esrf::Array{Tf, 2} = zeros(Nx, Ny)                       # Moisture flux from the surface (kg/m^2/s)
    Eveg::Array{Tf, 2} = zeros(Nx, Ny)                       # Moisture flux from vegetation (kg/m^2/s)
    G::Array{Tf, 2} = zeros(Nx, Ny)                          # Heat flux into the surface (W/m^2)
    H::Array{Tf, 2} = zeros(Nx, Ny)                          # Sensible heat flux to the atmosphere (W/m^2)
    Hsrf::Array{Tf, 2} = zeros(Nx, Ny)                       # Sensible heat flux from the surface (W/m^2)
    LE::Array{Tf, 2} = zeros(Nx, Ny)                         # Latent heat flux to the atmosphere (W/m^2)
    LEsrf::Array{Tf, 2} = zeros(Nx, Ny)                      # Latent heat flux from the surface (W/m^2)
    LWsci::Array{Tf, 2} = zeros(Nx, Ny)                      # Subcanopy incoming longwave radiation (W/m^2)
    LWveg::Array{Tf, 2} = zeros(Nx, Ny)                      # Net longwave radiation absorbed by vegetation (W/m^2)
    Melt::Array{Tf, 2} = zeros(Nx, Ny)                       # Surface melt rate (kg/m^2/s)
    Rnet::Array{Tf, 2} = zeros(Nx, Ny)                       # Net radiation (W/m^2)
    Rsrf::Array{Tf, 2} = zeros(Nx, Ny)                       # Net radiation at surface (W/m^2)

    # Variables used in ebalfor-function

    A_ebal::Array{Tf, 2} = zeros(4, 4)                       # Energy balance matrix for forest
    Acp_ebal::Array{Tf, 2} = zeros(4, 4)                     # Copy of energy balance matrix for LU decomposition
    b_ebal::Vector{Tf} = zeros(4)                            # Right-hand side vector for energy balance
    x_ebal::Vector{Tf} = zeros(4)                            # Solution vector for energy balance
    vv_ebal::Vector{Tf} = zeros(4)                           # Scaling vector for LU decomposition
    indx_ebal::Vector{Ti} = zeros(4)                         # Pivot indices for LU decomposition

    # Variables used in canopy-function

    intcpt::Array{Tf, 2} = zeros(Nx, Ny)                     # Canopy interception (kg/m^2)
    Sbveg::Array{Tf, 2} = zeros(Nx, Ny)                      # Sublimation from vegetation (kg/m^2)
    unload::Array{Tf, 2} = zeros(Nx, Ny)                     # Snow mass unloaded from canopy (kg/m^2)

    # Variables used in snow-function

    Gsoil::Array{Tf, 2} = zeros(Nx, Ny)                      # Heat flux into soil (W/m^2)
    Roff::Array{Tf, 2} = zeros(Nx, Ny)                       # Total runoff (kg/m^2)
    meltflux_out::Array{Tf, 2} = zeros(Nx, Ny)               # Runoff from snowmelt at base of snow (kg/m^2)
    Sbsrf::Array{Tf, 2} = zeros(Nx, Ny)                      # Sublimation from the snow surface (kg/m^2)
    Roff_bare::Array{Tf, 2} = zeros(Nx, Ny)                  # Bare soil runoff (kg/m^2)
    Roff_snow::Array{Tf, 2} = zeros(Nx, Ny)                  # Runoff at base of snow (kg/m^2)
    snowdepth0::Array{Tf, 2} = zeros(Nx, Ny)                 # Snow depth at start of timestep (m)
    Sice0::Array{Tf, 2} = zeros(Nx, Ny)                      # Ice content at start of timestep (kg/m^2)

    a::Vector{Tf} = zeros(Nsmax)                             # Tridiagonal matrix lower diagonal
    bsnow::Vector{Tf} = zeros(Nsmax)                         # Tridiagonal matrix main diagonal
    c::Vector{Tf} = zeros(Nsmax)                             # Tridiagonal matrix upper diagonal
    csnow::Vector{Tf} = zeros(Nsmax)                         # Areal heat capacity of snow layers (J/K/m^2)
    dTssnow::Vector{Tf} = zeros(Nsmax)                       # Snow layer temperature increments (K)
    D::Vector{Tf} = zeros(Nsmax)                             # Layer thickness (m)
    E::Vector{Tf} = zeros(Nsmax)                             # Energy flux (W/m^2)
    Gs::Vector{Tf} = zeros(Nsmax)                            # Inter-layer thermal conductance (W/m^2/K)
    rhs::Vector{Tf} = zeros(Nsmax)                           # Right-hand side for tridiagonal solver
    R::Vector{Tf} = zeros(Nsmax)                             # Liquid water flux between layers (kg/m^2/s)
    S::Vector{Tf} = zeros(Nsmax)                             # Layer source term (W/m^2)
    U::Vector{Tf} = zeros(Nsmax)                             # Layer internal energy (J/m^2)
    W::Vector{Tf} = zeros(Nsmax)                             # Layer liquid water content (kg/m^2)

    SWEbuffer::Vector{Tf} = zeros(15)                        # Buffer for SWE history (kg/m^2)
    snowdepthbuffer::Vector{Tf} = zeros(15)                  # Buffer for snow depth history (m)
    diffSWEbuffer::Vector{Tf} = zeros(14)                    # Buffer for SWE differences (kg/m^2)

    # Variables used in snow_layering-function

    Ds0::Array{Tf, 2} = zeros(Nx, Ny)                        # Snow layer thickness at start of timestep (m)
    hw::Vector{Tf} = zeros(Nsmax)                            # Liquid water equivalent height (m)
    rho::Vector{Tf} = zeros(Nsmax + 1)                       # Snow density (kg/m^3)
    diff_rho::Vector{Tf} = zeros(Nsmax)                      # Density difference between layers (kg/m^3)
    csnow_loc::Vector{Tf} = zeros(Nsmax + 1)                 # Local heat capacity of snow layers (J/K/m^2)
    Sice_loc::Vector{Tf} = zeros(Nsmax + 1)                  # Local ice content of snow layers (kg/m^2)
    Sliq_loc::Vector{Tf} = zeros(Nsmax + 1)                  # Local liquid content of snow layers (kg/m^2)
    Ds_loc::Vector{Tf} = zeros(Nsmax + 1)                    # Local snow layer thicknesses (m)
    histowet_loc::Vector{Tf} = zeros(Nsmax + 1)              # Local historical wetting variable (-)
    U_loc::Vector{Tf} = zeros(Nsmax + 1)                     # Local layer internal energy (J/m^2)
    Tsnow_loc::Vector{Tf} = zeros(Nsmax + 1)                 # Local snow layer temperatures (K)

    # Variables used in soil-function

    asoil::Vector{Tf} = zeros(Nsoil)                         # Tridiagonal matrix lower diagonal for soil
    bsoil::Vector{Tf} = zeros(Nsoil)                         # Tridiagonal matrix main diagonal for soil
    cssoil::Vector{Tf} = zeros(Nsoil)                        # Tridiagonal matrix upper diagonal for soil
    dTssoil::Vector{Tf} = zeros(Nsoil)                       # Soil layer temperature increments (K)
    Gssoil::Vector{Tf} = zeros(Nsoil)                        # Inter-layer thermal conductance for soil (W/m^2/K)
    rhssoil::Vector{Tf} = zeros(Nsoil)                       # Right-hand side for soil tridiagonal solver

    # Variables used in tridiag-function

    gammasnow::Vector{Tf} = zeros(Nsmax)                     # Tridiagonal solver work array for snow
    gammasoil::Vector{Tf} = zeros(Nsoil)                     # Tridiagonal solver work array for soil

    # Variables used in snowslide-function

    dSWE_tot_slide::Array{Tf, 2} = zeros(Nx, Ny)             # Cumulated SWE change due to slides (kg/m^2)
    index_sorted_dem::Array{Ti, 2} = zeros(Ti, Nx * Ny, 2)   # Sorted indices of digital elevation model

    # Variables used in snowtran3d-function

    dSWE_tot_subl::Array{Tf, 2} = zeros(Nx, Ny)              # Cumulated SWE change due to sublimation (kg/m^2)
    dSWE_tot_salt::Array{Tf, 2} = zeros(Nx, Ny)              # Cumulated SWE change due to saltation (kg/m^2)
    dSWE_tot_susp::Array{Tf, 2} = zeros(Nx, Ny)              # Cumulated SWE change due to suspension (kg/m^2)

end

@with_kw mutable struct MET{Tf, Ti}

    # Domain size

    Nx::Ti = 1                                               # Size of first array dimension (rows)
    Ny::Ti = 1                                               # Size of second array dimension (columns)

    # Meteorological variables

    Sdir::Array{Tf, 2} = fill(NaN, Nx, Ny)                   # Direct shortwave radiation per inclined surface area (W/m^2)
    Sdif::Array{Tf, 2} = fill(NaN, Nx, Ny)                   # Diffuse shortwave radiation (W/m^2)
    Sdird::Array{Tf, 2} = fill(NaN, Nx, Ny)                  # Direct shortwave radiation per horizontal surface area (W/m^2)
    LW::Array{Tf, 2} = fill(NaN, Nx, Ny)                     # Incoming longwave radiation (W/m^2)
    Sf::Array{Tf, 2} = fill(NaN, Nx, Ny)                     # Snowfall rate (kg/m^2/s)
    Rf::Array{Tf, 2} = fill(NaN, Nx, Ny)                     # Rainfall rate (kg/m^2/s)
    Sf24h::Array{Tf, 2} = fill(NaN, Nx, Ny)                  # Total snowfall over 24h (kg/m^2)
    Ta::Array{Tf, 2} = fill(NaN, Nx, Ny)                     # Air temperature (K)
    RH::Array{Tf, 2} = fill(NaN, Nx, Ny)                     # Relative humidity (%)
    Ua::Array{Tf, 2} = fill(NaN, Nx, Ny)                     # Wind speed (m/s)
    Udir::Array{Tf, 2} = fill(NaN, Nx, Ny)                   # Wind direction (degrees, clockwise from North)
    Ps::Array{Tf, 2} = fill(NaN, Nx, Ny)                     # Surface air pressure (Pa)
    Tv::Array{Tf, 2} = fill(NaN, Nx, Ny)                     # Time-varying transmissivity for direct shortwave radiation (-)

    # Snowfall tracking variables

    Sf24h_f64::Array{Float64, 2} = zeros(Nx, Ny)             # Total snowfall over 24h (kg/m^2)  TODO intermediate variable using Float64 to match matlab/fortran code - remove later
    Sf_history_f64::Array{Float64, 3} = zeros(Nx, Ny, 24)    # History of snowfall over the last 24h (kg/m^2)  TODO using Float64 to match matlab/fortran code - change precision later

end
