@with_kw mutable struct FSM{Tf, Ti}

  # Maximum snow and soil layer thicknesses (m)

  Dzsnow::Vector{Tf} = [0.1, 0.2, 0.4]
  Dzsoil::Vector{Tf} = [0.1, 0.2, 0.4, 0.8]

  # Number of snow and soil layers

  Nsmax::Ti = length(Dzsnow)
  Nsoil::Ti = length(Dzsoil)
  
  # Domain size

  Nx::Ti = 1
  Ny::Ti = 1

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
  FSNRHO::Ti = 2                                           # Fresh snow density (0=fixed, 1=climate-dependent, 2=climate+elevation)
  ALRADT::Ti = 1                                           # Activate tuning of albedo decay as function of incoming direct SWR (0, 1)
  SNOPRP::Ti = 1                                           # Snow surface properties (0=constant, 1=elevation-dependent)
  SNTRAN::Ti = 0                                           # Snow transport (0, 1)
  SNSLID::Ti = 0                                           # Snow slides (0, 1)
  SNOLAY::Ti = 0                                           # Density-dependent layering (0, 1)
  HN_ON::Bool = false                                      # TODO remove? Activate the new snow model
  Z0PERT::Bool = false                                     # TODO remove? Activate z0 perturbations
  WCPERT::Bool = false                                     # TODO remove? Activate liquid water capacity perturbations
  FSPERT::Bool = false                                     # TODO remove? Activate fresh snow density perturbations
  ALPERT::Bool = false                                     # TODO remove? Activate albedo perturbations
  SLPERT::Bool = false                                     # TODO remove? Activate settling perturbations

  # Tile options

  TILE::String = "open"                                    # Tile type
  tthresh::Tf = 0.1                                        # Tile threshold

  # Numerical solution parameters

  Nitr = 4                                                 # Number of iterations for surface energy balance

  # Defaults for canopy parameters

  avg0::Tf = 0.1                                           # Snow-free vegetation albedo
  avgs::Tf = 0.4                                           # Snow-covered vegetation albedo
  cden::Tf = 0.004                                         # Dense canopy turbulent transfer coefficient
  cvai::Tf = 4.4                                           # Canopy snow capacity per unit VAI (kg/m^2)
  cveg::Tf = 20                                            # Vegetation turbulent transfer coefficient
  Gcn1::Tf = 0.5                                           # Leaf angle distribution parameter
  Gcn2::Tf = 0                                             # Leaf angle distribution parameter
  gsnf::Tf = 0                                             # Snow-free vegetation moisture conductance (m/s)
  kdif::Tf = 0.5                                           # Diffuse radiation extinction coefficient
  kveg::Tf = 1                                             # Canopy cover coefficient
  rchd::Tf = 0.67                                          # Ratio of displacement height to canopy height
  rchz::Tf = 0.2                                           # Ratio of roughness length to canopy height
  tcnc::Tf = 3600*240                                      # Canopy unloading time scale for cold snow (s)
  tcnm::Tf = 3600*48                                       # Canopy unloading time scale for melting snow (s)
  pmultf_for::Tf = 0.5                                         # TODO add description what this is...

  # Defaults for snow parameters

  a_eta::Tf = 0.1                                          # Temperature factor for Crocus B92 compaction (K^-1)
  asmx::Tf = 0.86                                          # Maximum albedo for fresh snow
  asmn::Tf = 0.6                                           # Minimum albedo for melting snow
  b_eta::Tf = 0.023                                        # First density factor for Crocus B92 compaction (m^3/kg)
  bthr::Tf = 2                                             # Snow thermal conductivity exponent
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
  rgr0::Tf = 5e-5                                          # Fresh snow grain radius (m)
  rmlt::Tf = 500                                           # Maximum density for melting snow (kg/m^3)
  Salb::Tf = 10                                            # Albedo decay constant (kg/m^2)
  snda::Tf = 2.8e-6                                        # Thermal metamorphism parameter (1/s)
  Talb::Tf = -2                                            # Albedo decay temperature threshold (C)
  tcld::Tf = 3600*1000                                     # Cold snow albedo decay time scale (s)
  tmlt::Tf = 3600*100                                      # Melting snow albedo decay time scale (s)
  trho::Tf = 3600*200                                      # Snow compaction time scale (s)
  Wirr::Tf = 0.03                                          # Irreducible liquid water content of snow
  z0sn::Tf = 0.002                                         # Snow roughness length (m)
  Sfmin::Tf = 10                                           # Minimum 24h snowfall to refresh albedo (kg/m^2)

  # Snow layering parameters
  Ds_min::Tf = 0.01                                        # Minimum possible snow layer thickness (m)
  Ds_surflay::Tf = 0.5                                     # Maximum thickness of surface fine snow layering (m)

  # SnowSlide parameters
  dyn_ratio::Tf = 0.09                                     # Dynamic snow holding depth ratio
  rho_deposit::Tf = 300.0                                  # Constant snow avalanche deposit density (kg/m³)
  slope_min::Tf = 30.0                                     # Minimum slope for snow slide occurrence (deg)
  Shd_min::Tf = 0.01                                       # Minimum snow holding depth (m)
  rho_snow::Tf = 300.0                                     # Constant snow density for transport (kg/m³)  TODO comes from PARAM_SNOWTRAN3D, should it be here?

  # Defaults for ground surface parameters

  bstb::Tf = 5                                             # Atmospheric stability parameter
  gsat::Tf = 0.01                                          # Surface conductance for saturated soil (m/s)

  # Defaults for additional forest snow process parametrization

  adfs::Tf = 3                                             # Snow albedo adjustment dependent on SWR
  adfl::Tf = 2                                             # Snow albedo adjustment dependent on LWR
  fsar::Tf = 0.1                                           # Snow albedo adjustment range dependent on vegetation fraction
  psf::Tf  = 1                                             # Scaling factor for solid precipitation (within forest stand, at min CC)
  psr::Tf  = 0.1                                           # Range of solid precipitation (within forest stand, spread min-max CC)
  wcan::Tf = 2.5                                           # Parameter of exponential wind profile
  zsub::Tf = 2                                             # Sub-canopy reference height (m)
  zgf::Tf = 1                                              # Roughness length adjustment factor depending on vegetation fraction
  zgr::Tf = 0                                              # Roughness length adjustment range depending on vegetation fraction
  khcf::Tf = 3                                             # Diffusivity adjustment for canopy effects (Finnigan 2000)

  # Surface parameters

  adm::Tf = -999999                                        # TODO defaults? Melting snow albedo decay time (h)
  adc::Array{Tf,2} = -999999*ones(Nx,Ny)                   # TODO defaults? Cold snow albedo decay time (h)
  afs::Array{Tf,2} = -999999*ones(Nx,Ny)                   # TODO defaults? Maximum albedo for fresh snow
  z0_snow::Array{Tf,2} = -999999*ones(Nx,Ny)               # TODO defaults? Roughness length of snow (m)

  # Surface properties

  alb0::Array{Tf,2} = 0.2*ones(Nx,Ny)                      # Snow-free ground albedo
  z0sf::Array{Tf, 2} = 0.2*ones(Nx,Ny)                     # Snow-free roughness length (m)
  fcly::Array{Tf, 2} = 0.3*ones(Nx,Ny)                     # Soil clay fraction
  fsnd::Array{Tf, 2} = 0.6*ones(Nx,Ny)                     # Soil sand fraction

  # Canopy parameters

  canh::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Canopy heat capacity (J/K/m^2)
  fsky::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Sky view fraction
  fveg::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Canopy cover fraction
  fves::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Stand-scale canopy cover fraction
  hcan::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Canopy height (m)
  lai::Array{Tf, 2} = -999999*ones(Nx,Ny)                  # TODO defaults? Leaf area index
  pmultf::Array{Tf, 2} = -999999*ones(Nx,Ny)               # TODO defaults? Precip multiplier to revert precip correction applied to open area
  scap::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Canopy snow capacity (kg/m^2)
  trcn::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Canopy transmissivity
  VAI::Array{Tf, 2} = -999999*ones(Nx,Ny)                  # TODO defaults? Vegetation area index
  vfhp::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Hemispherical sky-view fraction including canopy

  # Terrain properties

  slopemu::Array{Tf, 2} = -999999*ones(Nx,Ny)              # TODO defaults? Slope parameter
  xi::Array{Tf, 2} = -999999*ones(Nx,Ny)                   # TODO defaults? Terrain correlation length
  Ld::Array{Tf, 2} = -999999*ones(Nx,Ny)                   # TODO defaults? Grid cell size or domain size (m)
  fsky_terr::Array{Tf, 2} = -999999*ones(Nx,Ny)            # TODO defaults? Sky view fraction terrain
  lat::Array{Tf, 2} = -999999*ones(Nx,Ny)                  # TODO defaults? Latitude of each grid cell (center?)
  lon::Array{Tf, 2} = -999999*ones(Nx,Ny)                  # TODO defaults? Longitude of each grid cell (center?)
  dem::Array{Tf, 2} = -999999*ones(Nx,Ny)                  # TODO defaults? Grid elevation (m)
  tilefrac::Array{Tf,2} = -999999*ones(Nx,Ny)              # TODO defaults? Tile fraction
  glacierfrac::Array{Tf,2} = -999999*ones(Nx,Ny)           # TODO defaults? Glacier flag
  vegsnowd_xy::Array{Tf,2} = -999999*ones(Nx,Ny)           # TODO defaults? Vegetation snow holding capacity (m)

  prec_multi::Array{Float64,2} = -999999*ones(Nx,Ny)       # TODO defaults? Precipitation multiplier (-)    TODO HACK FLOAT64

  # Derived soil parameters

  b::Array{Tf, 2} = zeros(Nx,Ny)                           # Clapp-Hornberger exponent
  hcap_soil::Array{Tf, 2} = zeros(Nx,Ny)                   # Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil::Array{Tf, 2} = zeros(Nx,Ny)                   # Thermal conductivity of dry soil (W/m/K)
  sathh::Array{Tf, 2} = zeros(Nx,Ny)                       # Saturated soil water pressure (m)
  Vsat::Array{Tf, 2} = zeros(Nx,Ny)                        # Volumetric soil moisture at saturation
  Vcrit::Array{Tf, 2} = zeros(Nx,Ny)                       # Volumetric soil moisture at critical point

  # State variables

  albs::Array{Tf, 2} = Tf(0.85)*ones(Nx,Ny)                # Snow albedo
  Ds::Array{Tf,3} = zeros(Nsmax,Nx,Ny)                     # Snow layer thicknesses (m)
  Nsnow::Array{Ti,2} = zeros(Ti,Nx,Ny)                     # Number of snow layers
  Qcan::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Canopy air space humidity
  rgrn::Array{Tf,3} = -999999*ones(Nsmax,Nx,Ny)            # TODO defaults? Snow layer grain radius (m)
  Sice::Array{Tf,3} = zeros(Nsmax,Nx,Ny)                   # Ice content of snow layers (kg/m^2)
  Sliq::Array{Tf,3} = zeros(Nsmax,Nx,Ny)                   # Liquid content of snow layers (kg/m^2)
  Sveg::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Snow mass on vegetation (kg/m^2)
  Tcan::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Canopy air space temperature (K)
  theta::Array{Tf,3} = zeros(Nsoil,Nx,Ny)                  # Volumetric moisture content of soil layers
  Tsnow::Array{Tf,3} = Tf(273.15)*ones(Nsmax,Nx,Ny)        # Snow layer temperatures (K)
  Tsoil::Array{Tf,3} = Tf(285)*ones(Nsoil,Nx,Ny)           # Soil layer temperatures (K)
  Tsrf::Array{Tf, 2} = Tf(285)*ones(Nx,Ny)                 # Surface skin temperature (K)
  fsnow::Array{Tf, 2} = zeros(Nx,Ny)                       # Snow cover fraction terrain
  Tveg::Array{Tf, 2} = -999999*ones(Nx,Ny)                 # TODO defaults? Vegetation temperature (K)
  snowdepthmin::Array{Tf, 2} = zeros(Nx,Ny)                # Minimum snow depth at time step of swemin (m)
  snowdepthmax::Array{Tf, 2} = zeros(Nx,Ny)                # Maximum snow depth at time stemp of swemax (m)
  snowdepthhist::Array{Tf,3} = zeros(14,Nx,Ny)             # History of snow depth during last 14 days (m). Most recent entries first
  swemin::Array{Tf, 2} = zeros(Nx,Ny)                      # Minimum SWE during the season (m)
  swemax::Array{Tf, 2} = zeros(Nx,Ny)                      # Maximum SWE during the season (m)
  swehist::Array{Tf,3} = zeros(14,Nx,Ny)                   # History of SWE during last 14 days (kg/m^2). Most recent entries first
  histowet::Array{Tf,3} = -999999*ones(Nsmax,Nx,Ny)        # TODO defaults? Historical variable for past wetting of a layer (0-1)
  dSWE_tot_subl::Array{Tf, 2} = zeros(Nx,Ny)               # Cumulated SWE change due to sublimation (kg/m^2)
  dSWE_tot_salt::Array{Tf, 2} = zeros(Nx,Ny)               # Cumulated SWE change due to saltation (kg/m^2)
  dSWE_tot_susp::Array{Tf, 2} = zeros(Nx,Ny)               # Cumulated SWE change due to suspension (kg/m^2)

  # Radiation - temporary arrays

  alb::Array{Tf,2} = zeros(Nx,Ny)
  asrf_out::Array{Tf,2} = zeros(Nx,Ny)
  SWveg::Array{Tf,2} = zeros(Nx,Ny)
  SWsrf::Array{Tf,2} = zeros(Nx,Ny)
  SWsci::Array{Tf,2} = zeros(Nx,Ny)
  LWt::Array{Tf,2} = zeros(Nx,Ny)

  # Thermal - temporary arrays

  ksnow::Array{Tf,3} = zeros(Nsmax, Nx, Ny)
  csoil::Array{Tf,3} = zeros(Nsoil, Nx, Ny)
  ksoil::Array{Tf,3} = zeros(Nsoil, Nx, Ny)
  gs1::Array{Tf,2} = zeros(Nx, Ny)
  Ds1::Array{Tf,2} = zeros(Nx, Ny)
  Ts1::Array{Tf,2} = zeros(Nx, Ny)
  ks1::Array{Tf,2} = zeros(Nx, Ny)
  Tveg0::Array{Tf,2} = zeros(Nx, Ny)

  # Sfexch - temporary arrays

  KH::Array{Tf,2} = zeros(Nx, Ny)
  KHa::Array{Tf,2} = zeros(Nx, Ny)
  KHg::Array{Tf,2} = zeros(Nx, Ny)
  KHv::Array{Tf,2} = zeros(Nx, Ny)
  KWg::Array{Tf,2} = zeros(Nx, Ny)
  KWv::Array{Tf,2} = zeros(Nx, Ny)
  Usc::Array{Tf,2} = zeros(Nx, Ny)

  # Ebalsrf - temporary arrays

  Esrf::Array{Tf,2} = zeros(Nx,Ny)
  Eveg::Array{Tf,2} = zeros(Nx,Ny)
  G::Array{Tf,2} = zeros(Nx,Ny)
  H::Array{Tf,2} = zeros(Nx,Ny)
  Hsrf::Array{Tf,2} = zeros(Nx,Ny)
  LE::Array{Tf,2} = zeros(Nx,Ny)
  LEsrf::Array{Tf,2} = zeros(Nx,Ny)
  LWsci::Array{Tf,2} = zeros(Nx,Ny)
  LWveg::Array{Tf,2} = zeros(Nx,Ny)
  Melt::Array{Tf,2} = zeros(Nx,Ny)
  Rnet::Array{Tf,2} = zeros(Nx,Ny)
  Rsrf::Array{Tf,2} = zeros(Nx,Ny)
  Ssub::Tf = zero(Tf)

  # Ebalfor - temporary arrays

  A_ebal::Array{Tf,2} = zeros(4, 4)
  Acp_ebal::Array{Tf,2} = zeros(4, 4)
  b_ebal::Vector{Tf} = zeros(4)
  x_ebal::Vector{Tf} = zeros(4)
  vv_ebal::Vector{Tf} = zeros(4)
  indx_ebal::Vector{Ti} = zeros(4)

  # Canopy - temporary arrays

  intcpt::Array{Tf,2} = zeros(Nx, Ny)
  Sbveg::Array{Tf,2} = zeros(Nx, Ny)
  unload::Array{Tf, 2} = zeros(Nx, Ny)

  # Snow - temporary arrays

  Gsoil::Array{Tf,2} = zeros(Nx, Ny)
  Roff::Array{Tf,2} = zeros(Nx, Ny)
  meltflux_out::Array{Tf,2} = zeros(Nx, Ny)
  Sbsrf::Array{Tf,2} = zeros(Nx, Ny)
  Roff_bare::Array{Tf, 2} = zeros(Nx, Ny)
  Roff_snow::Array{Tf, 2} = zeros(Nx, Ny)
  snowdepth0::Array{Tf, 2} = zeros(Nx, Ny)
  Sice0::Array{Tf, 2} = zeros(Nx, Ny)

  a::Vector{Tf} = zeros(Nsmax)
  bsnow::Vector{Tf} = zeros(Nsmax)
  c::Vector{Tf} = zeros(Nsmax)
  csnow::Vector{Tf} = zeros(Nsmax)
  dTssnow::Vector{Tf} = zeros(Nsmax)
  D::Vector{Tf} = zeros(Nsmax)
  E::Vector{Tf} = zeros(Nsmax)
  Gs::Vector{Tf} = zeros(Nsmax)
  rhs::Vector{Tf} = zeros(Nsmax)
  R::Vector{Tf} = zeros(Nsmax)
  S::Vector{Tf} = zeros(Nsmax)
  U::Vector{Tf} = zeros(Nsmax)
  W::Vector{Tf} = zeros(Nsmax)

  SWEbuffer::Vector{Tf} = zeros(15)
  snowdepthbuffer::Vector{Tf} = zeros(15)
  diffSWEbuffer::Vector{Tf} = zeros(14)

  # Snow layering - temporary arrays

  Ds0::Array{Tf, 2} = zeros(Nx, Ny)
  hw::Vector{Tf} = zeros(Nsmax)
  rho::Vector{Tf} = zeros(Nsmax + 1)
  diff_rho::Vector{Tf} = zeros(Nsmax)
  csnow_loc::Vector{Tf} = zeros(Nsmax + 1)
  Sice_loc::Vector{Tf} = zeros(Nsmax + 1)
  Sliq_loc::Vector{Tf} = zeros(Nsmax + 1)
  Ds_loc::Vector{Tf} = zeros(Nsmax + 1)
  histowet_loc::Vector{Tf} = zeros(Nsmax + 1)
  U_loc::Vector{Tf} = zeros(Nsmax + 1)
  Tsnow_loc::Vector{Tf} = zeros(Nsmax + 1)

  # Soil - temporary vectors

  asoil::Vector{Tf} = zeros(Nsoil)
  bsoil::Vector{Tf} = zeros(Nsoil)
  cssoil::Vector{Tf} = zeros(Nsoil)
  dTssoil::Vector{Tf} = zeros(Nsoil)
  Gssoil::Vector{Tf} = zeros(Nsoil)
  rhssoil::Vector{Tf} = zeros(Nsoil)

  # Tridiag - temporary vectors

  gammasnow::Vector{Tf} = zeros(Nsmax)
  gammasoil::Vector{Tf} = zeros(Nsoil)

  # SnowSlide variables     TODO where to place these?
  slope::Matrix{Tf} = Matrix{Tf}(undef, 0, 0)             # Slope angles (deg)    TODO move to terrain params
  Shd::Matrix{Tf} = Matrix{Tf}(undef, 0, 0)               # Snow holding depth (m)    TODO move to terrain params
  dSWE_tot_slide::Matrix{Tf} = Matrix{Tf}(undef, 0, 0)    # Cumulated SWE change due to slides    TODO this seems to be some output and not a state
  index_sorted_dem::Matrix{Ti} = Matrix{Ti}(undef, 0, 0)  # Sorted DEM indices    TODO this seems to be some static input to the algorithm that should be computed in setup

end

@with_kw mutable struct MET{Tf, Ti}
  
  # Domain size

  Nx::Ti = 1
  Ny::Ti = 1
  
  # Time variables

  year::Array{Ti, 2} = -999999*ones(1, 1)                  # TODO defaults? 
  month::Array{Ti, 2} = -999999*ones(1, 1)                 # TODO defaults?
  day::Array{Ti, 2} = -999999*ones(1, 1)                   # TODO defaults?
  hour::Array{Tf, 2} = -999999*ones(1, 1)                  # TODO defaults?

  # Meteorological variables

  Sdir::Array{Tf, 2} = zeros(Nx, Ny)                       # Direct-beam shortwave radiation (W/m^2)
  Sdif::Array{Tf, 2} = zeros(Nx, Ny)                       # Diffuse shortwave radiation (W/m^2)
  Sdird::Array{Tf, 2} = zeros(Nx, Ny)                      # Direct-beam shortwave radiation, per horizontal surface area (W/m^2)
  LW::Array{Tf, 2} = zeros(Nx, Ny)                         # Incoming longwave radiation (W/m^2)
  Sf::Array{Tf, 2} = zeros(Nx, Ny)                         # Snowfall rate (kg/m^2/s)
  Rf::Array{Tf, 2} = zeros(Nx, Ny)                         # Rainfall rate (kg/m^2/s)
  Ta::Array{Tf, 2} = zeros(Nx, Ny)                         # Air temperature (K)
  RH::Array{Tf, 2} = zeros(Nx, Ny)                         # Relative humidity (%)
  Ua::Array{Tf, 2} = zeros(Nx, Ny)                         # Wind speed (m/s)
  Udir::Array{Tf, 2} = zeros(Nx, Ny)                       # Wind direction (degrees, clockwise from N)
  Ps::Array{Tf, 2} = zeros(Nx, Ny)                         # Surface pressure (Pa)
  Sf24h::Array{Tf, 2} = zeros(Nx, Ny)                      # Snowfall 24hr (kg/m^2)
  es::Array{Tf, 2} = zeros(Nx, Ny)                         # Saturation vapour pressure (Pa)
  Qa::Array{Tf, 2} = zeros(Nx, Ny)                         # Specific humidity (kg/kg)
  Tv::Array{Tf, 2} = ones(Nx, Ny)                          # Time-varying transmissivity for direct SWR (-)

  # Snowfall tracking variable

  Sf24h_f64::Array{Float64, 2} = zeros(Nx, Ny)             # TODO use Float64 to better match results from matlab/fortran code base
  Sf_history_f64::Array{Float64, 3} = zeros(Nx, Ny, 24)    # TODO use Float64 to better match results from matlab/fortran code base

end
