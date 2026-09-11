!-----------------------------------------------------------------------
! Physical constants
!-----------------------------------------------------------------------
module CONSTANTS
real, parameter :: &
  cp = 1005,         &! Specific heat capacity of air (J/K/kg)
  eps = 0.622,       &! Ratio of molecular weights of water and dry air
  e0 = 610.78,       &! Saturation vapour pressure at Tm (Pa)
  grav = 9.81,       &! Acceleration due to gravity (m/s^2)
  hcap_ice = 2100,   &! Specific heat capacity of ice (J/K/kg)
  hcap_wat = 4180,   &! Specific heat capacity of water (J/K/kg)
  hcon_air = 0.025,  &! Thermal conductivity of air (W/m/K)
  hcon_clay = 1.16,  &! Thermal conductivity of clay (W/m/K)
  hcon_ice = 2.24,   &! Thermal conducivity of ice (W/m/K)
  hcon_sand = 1.57,  &! Thermal conductivity of sand (W/m/K)
  hcon_wat = 0.56,   &! Thermal conductivity of water (W/m/K)
  I0 = 1367,         &! Solar constant (W/m^2)
  Lf = 0.334e6,      &! Latent heat of fusion (J/kg)
  Lv = 2.501e6,      &! Latent heat of vapourisation (J/kg)
  Ls = Lf + Lv,      &! Latent heat of sublimation (J/kg)
  pi = 3.14159,      &! pi
  Rair = 287,        &! Gas constant for air (J/K/kg)
  Rwat = 462,        &! Gas constant for water vapour (J/K/kg)
  Runi = 8313,       &! Universal gas constant (J/kmol/K)
  rho_ice = 917,     &! Density of ice (kg/m^3)
  rho_wat = 1000,    &! Density of water (kg/m^3)
  rho_air = 1.275,   &! Density of air (kg/m^3)
  sb = 5.67e-8,      &! Stefan-Boltzmann constant (W/m^2/K^4)
  em_snow = 0.99,    &! Emissivity snow for Stefan-Boltzmann
  em_soil = 0.90,    &! Emissivity soil for Stefan-Boltzmann
  Tm = 273.15,       &! Melting point (K)
  vkman = 0.4,       &! Von Karman constant
  xM = 18.01,        &! Molecular weight of water (kg/kmol)
  visc_air = 13.e-6, &! Kinematic viscosity of air (m^2/s)
  undef = 1.e+6       ! Initialization value for allocatables
integer, parameter :: &
  iundef = 1.e+6      ! Initialization value for integer allocatables
end module CONSTANTS

!-----------------------------------------------------------------------
! Model configuration
!-----------------------------------------------------------------------
module MODCONF
! Process options                                         : Possible values
integer :: &
  ALBEDO,     &! snow albedo                              : 0, 1, 2
  CANMOD,     &! forest canopy                            : 0, 1
  CONDCT,     &! snow thermal conductivity                : 0, 1
  DENSTY,     &! snow density                             : 0, 1, 2, 3
  EXCHNG,     &! turbulent exchange                       : 0, 1
  HYDROL,     &! snow hydraulics                          : 0, 1, 2
  SNFRAC,     &! snow cover fraction                      : 0, 1, 2, 3, 4
! Driving data options
  ZOFFST,     &! measurement height offset                : 0, 1
! OSHD-specific options
  OSHDTN,     &! oshd-specific tuning options             : 0, 1
               ! of fresh snow albedo, snow roughness lengths and fresh snow density
  ALRADT,     &! activate tuning of the albedo decay      : 0, 1
               ! as a function of incoming direct shortwave radiation
! FSM2trans specific options
  SNTRAN,     &! snow transport                           : 0, 1
  SNSLID,     &! snow slides                              : 0, 1
  SNOLAY,     &! density-dependent layering               : 0, 1
  WNDENS,     &! wnd compaction of snow surface           : 0, 1, 2
  CHECKS       ! check state variables at every time step : 0, 1, 2
logical :: &
  HN_ON,      &!activate the new snow model
  FOR_HN       !write 18h states for the hn model.
end module MODCONF

module MODPERT
logical:: &
  Z0PERT,     &! activate z0 perturbations
  WCPERT,     &! activate liquid water capacity perturbations
  FSPERT,     &! activate fresh snow density perturbations
  ALPERT,     &! activate perturbations of (hourly) albedo changes
  AFPERT,     &! activate fresh snow albedo perturbations
  SLPERT       ! activate settling perturbations
end module MODPERT

!-----------------------------------------------------------------------
! Model tile
!-----------------------------------------------------------------------
module MODTILE
! Model tiles                                   : Possible values
logical:: &
  tiled_trans_run      ! Tiled trans run flag   : .true., .false.
character(len=20) :: &
  TILE           ! model tile                   : 'open', 'forest', 'glacier' 
real :: & 
  tthresh        ! Tile fraction of grid cell required for tile to be considered 
end module MODTILE 

!-----------------------------------------------------------------------
! Model output configuration
!-----------------------------------------------------------------------
module MODOUTPUT
! list of diagnostic variables that can be written to output bin files.
! at the moment, only 2d real variables are handled.
character(len=3), dimension(39) :: &
  WRITE_DIAG_VARS = (/&
                      'asf', &  ! asrf_out     Surface albedo
                      'ava', &  ! dSWE_slide   SWE change due to snow slides (kg/m^2)
                      'eml', &  ! Melt         Surface melt rate (kg/m^2/s)
                      'esf', &  ! Esrf         Moisture flux from the surface (kg/m^2/s)
                      'evg', &  ! Eveg         Moisture flux from vegetation (kg/m^2/s)
                      'gsl', &  ! Gsoil        Heat flux into soil (W/m^2)
                      'hsf', &  ! Hsrf         Sensible heat flux from the surface (W/m^2)
                      'hsn', &  ! snowdepth    Total snowdepth (m)
                      'int', &  ! intcpt       Canopy interception (kg/m^2)
                      'kha', &  ! KH           Eddy diffusivity for heat to the atmosphere (m/s)
                      'khc', &  ! KHa          Eddy diffusivity for heat from the canopy air space (m/s)
                      'khg', &  ! KHg          Eddy diffusivity for heat from the ground (m/s)
                      'khv', &  ! KHv          Eddy diffusivity for heat from vegetation (m/s)
                      'kwg', &  ! KWg          Eddy diffusivity for water from the ground (m/s)
                      'kwv', &  ! KWv          Eddy diffusivity for water from vegetation (m/s)
                      'lha', &  ! LE           Latent heat flux to the atmosphere (W/m^2)
                      'lhs', &  ! LEsrf        Latent heat flux from the surface (W/m^2)
                      'lwc', &  ! LWsci        Subcanopy incoming longwave radiation (W/m^2)
                      'lwr', &  ! Lwt          Incoming longwave radiation corrected for subgrid topography (W/m^2)
                      'lwv', &  ! LWveg        Net longwave radiation absorbed by vegetation (W/m^2)
                      'rnt', &  ! Rnet         Net radiation (W/m^2)
                      'rom', &  ! meltflux_out Runoff from snowmelt at the base of snow (kg/m^2)
                      'rot', &  ! Roff         Total runoff, snow and bare soil (kg/m^2)
                      'rsf', &  ! Rsrf         Net radiation at surface (W/m^2)
                      'sal', &  ! dSWE_salt    SWE change due to saltation (kg/m^2)
                      'sbl', &  ! dSWE_subl    SWE change due to snowdrift sublimation (kg/m^2)
                      'sbs', &  ! Sbsrf        Sublimation from the snow surface (kg/m^2)
                      'sbv', &  ! Sbveg        Sublimation from vegetation (kg/m^2)
                      'sdf', &  ! Sdift        Incoming diffuse beam radiation corrected for subgrid topography (W/m^2)
                      'sdi', &  ! Sdirt        Incoming direct beam radiation corrected for subgrid topography (W/m^2)
                      'sha', &  ! H            Sensible heat flux to the atmosphere (W/m^2)
                      'slq', &  ! Sliq_out     Total LWC (kg/m^2)
                      'ssp', &  ! dSWE_susp    SWE change due to suspension (kg/m^2)
                      'swc', &  ! SWsci        Subcanopy incoming shortwave radiation (W/m^2)
                      'swe', &  ! SWE          Total SWE (kg/m^2)
                      'sws', &  ! SWsrf        Net SW radiation absorbed by the surface (W/m^2)
                      'swv', &  ! SWveg        Net SW radiation absorbed by vegetation (W/m^2)
                      'uac', &  ! Usc          Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)
                      'uld'  &  ! unload       Snow mass unloaded from canopy (kg/m^2)
                      /)
   ! GM Note: Energy-balance relevant diagnostics omitted in this list can be derived from existing variables:
                    !Hveg = H-Hsrf
                    !LEveg = LE-LEsrf
                    !LWsrf = Tss^4*sb
                    !Rveg  = SWveg + LWveg      
! list of state variables that can be written to output bin files.
! at the moment, only 2d real variables are handled.
character(len=3), dimension(15):: &
  WRITE_STATE_VARS = (/&
                       'als', &  ! albs
                       'avt', &  ! dSWE_tot_slide
                     ! 'dsn', &  ! Ds(3d var)
                     ! 'hsh', &  ! snowdepthhist (3d var, 1st dimension: time.)
                       'hmn', &  ! snowdepthmin
                       'hmx', &  ! snowdepthmax
                     ! 'nsn', &  ! Nsnow
                       'qca', &  ! Qcan
                       'sbt', &  ! dSWE_tot_subl
                       'scf', &  ! fsnow
                     ! 'sic', &  ! Sice (3d var)
                     ! 'slq', &  ! Sliq (3d var)
                       'sat', &  ! dSWE_tot_salt
                       'smn', &  ! swemin
                       'smx', &  ! swemax
                       'spt', &  ! dSWE_tot_susp
                       'sve', &  ! Sveg
                     ! 'swh', &  ! swehist (3d var, 1st dimension: time.)
                       'tca', &  ! Tcan
                       'tsf', &  ! Tsrf
                     ! 'tsl', &  ! Tsoil (3d var)
                     ! 'tsn', &  ! Tsnow (3d var)
                       'tve'  &  ! Tveg
                       /)
character(len=3), dimension(39) :: &
  LIST_DIAG_RESULTS  ! List of diagnostic variables that the user wants to write into output bin files (subset of WRITE_DIAG_VARS)
character(len=3), dimension(15) :: &
  LIST_STATE_RESULTS ! List of result variables that the user wants to write into output bin files (subset of WRITE_STATE_VARS)
logical, dimension(39) :: WRITE_DIAG_TABLE  ! table specifying whether the corresponding variables in WRITE_DIAG_VARS should be written or not.
logical, dimension(15) :: WRITE_STATE_TABLE ! table specifying whether the corresponding variables in WRITE_STATE_VARS should be written or not.
end module MODOUTPUT

!-----------------------------------------------------------------------
! Output diagnostics
!-----------------------------------------------------------------------
module DIAGNOSTICS
integer :: &
  Nave                ! Number of timesteps in average outputs
end module DIAGNOSTICS

!-----------------------------------------------------------------------
! Meteorological driving variables
!-----------------------------------------------------------------------
module DRIVING
integer :: &
  year,              &! Year
  month,             &! Month of year
  day                 ! Day of month
real :: &
  hour                ! Hour of day (UTC+1)
real :: &
  dt,                &! Timestep (s)
  zT,                &! Temperature measurement height (m)
  zU,                &! Wind speed measurement height (m)
  zRH                 ! Relative humidity measurement height (m)
real, allocatable :: &
  LW(:,:),           &! Incoming longwave radiation (W/m^2)
  Ps(:,:),           &! Surface pressure (Pa)
  Qa(:,:),           &! Specific humidity (kg/kg)
  RH(:,:),           &! Relative humidity (%)
  Rf(:,:),           &! Rainfall rate (kg/m^2/s)
  Sf(:,:),           &! Snowfall rate (kg/m^2/s)
  Sf24h(:,:),        &! Snowfall 24hr (kg/m^2)
  Sdif(:,:),         &! Diffuse shortwave radiation (W/m^2)
  Sdir(:,:),         &! Direct-beam shortwave radiation (W/m^2)
  Ta(:,:),           &! Air temperature (K)
  Tv(:,:),           &! Time-varying transmissivity for direct SWR (-)
  Ua(:,:),           &! Wind speed (m/s)
  Udir(:,:),         &! Wind direction (degrees, clockwise from N)
  z0P(:,:),          &! z0 perturbations
  wcP(:,:),          &! Liquid water capacity perturbations
  fsP(:,:),          &! Fresh snow density perturbations
  alP(:,:),          &! Albedo perturbations
  afP(:,:),          &! Fresh snow albedo perturbations
  slP(:,:),          &! Settling perturbations
  Sdird(:,:)          ! Direct-beam shortwave radiation, per horizontal surface area (W/m2)
end module DRIVING

!-----------------------------------------------------------------------
! Grid parameters
!-----------------------------------------------------------------------
module GRID
integer :: &
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions
real :: &
  Ds_min,            &! Minimum possible snow layer thickness (m)
  Ds_surflay          ! Maximum thickness of surface fine snow layering (m)
real, allocatable :: &
  Dzsnow(:),         &! Minimum snow layer thicknesses (m)
  Dzsoil(:)           ! Soil layer thicknesses (m)
end module GRID

!-----------------------------------------------------------------------
! Input / output unit numbers
!-----------------------------------------------------------------------
! module IOUNITS
! integer, parameter :: &
!   umta = 51           ! Metadata output file unit number
! end module IOUNITS

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
module PARAMETERS
! Numerical solution parameter
integer :: &
  Nitr                ! Number of iterations in energy balance calulation
! Vegetation parameters
real :: &
  avg0,              &! Snow-free vegetation albedo
  avgs,              &! Snow-covered vegetation albedo
  cden,              &! Dense canopy turbulent transfer coefficient
  cvai,              &! Canopy snow capacity per unit VAI (kg/m^2)
  Gcn1,              &! Leaf angle distribution parameter
  Gcn2,              &! Leaf angle distribution parameter
  gsnf,              &! Snow-free vegetation moisture conductance (m/s)
  kdif,              &! Diffuse radiation extinction coefficient
  kveg,              &! Canopy cover coefficient
  cveg,              &! Vegetation turbulent transfer coefficient
  rchd,              &! Ratio of displacement height to canopy height
  rchz,              &! Ratio of roughness length to canopy height
  tcnc,              &! Canopy unloading time scale for cold snow (s)
  tcnm                ! Canopy unloading time scale for melting snow (s)
! Snow parameters
real :: &
  a_eta,             &! Temperature factor for Crocus B92 compaction (K^-1)
  adm,               &! Melting snow albedo decay time (h)
  asmx,              &! Maximum albedo for fresh snow
  asmn,              &! Minimum albedo for melting snow
  b_eta,             &! First density factor for Crocus B92 compaction (m^3/kg)
  bstb,              &! Atmospheric stability parameter
  bthr,              &! Snow thermal conductivity exponent
  c_eta,             &! Second density factor for Crocus B92 compaction (kg/m^3)
  eta0,              &! Reference snow viscosity (Pa s)
  eta1,              &! Reference snow viscosity for Crocus B92 compaction (Pa s)
  hfsn,              &! Snowcover fraction depth scale (m)
  kfix,              &! Fixed thermal conductivity of snow (W/m/K)
  rho0,              &! Fixed snow density (kg/m^3)
  rhob,              &! Temperature factor in fresh snow density (kg/m^3/K)
  rhoc,              &! Wind factor in fresh snow density (kg s^0.5/m^3.5)
  rhof,              &! Fresh snow density (kg/m^3)
  rhos_min,          &! Minimum snow density (kg/m^3)
  rhos_max,          &! Maximum snow density (kg/m^3)
  rcld,              &! Maximum density for cold snow (kg/m^3)
  rmlt,              &! Maximum density for melting snow (kg/m^3)
  snda,              &! Thermal metamorphism parameter (1/s)
  Talb,              &! Albedo decay temperature threshold (C)
  tcld,              &! Cold snow albedo decay time scale (s)
  tmlt,              &! Melting snow albedo decay time scale (s)
  trho,              &! Snow compaction time scale (s)
  Wirr,              &! Irreducible liquid water content of snow
  z0sn,              &! Snow roughness length (m)
  Sfmin,             &! Minimum 24h snowfall to refresh albedo (kg/m^2)
  U_thresh,          &! Threshold wind speed for wind compaction of surface snow (m/s)
  A_enh0,            &! Wind compaction enhancement parameter (s^3/m^3)
  d_max               ! Maximal depth affected by wind compaction (m)

! Surface parameters
real :: &
  gsat                ! Surface conductance for saturated soil (m/s)
! Additional parameters used in forest snow process parametrizations
real :: &      
  adfs,              &! Snow albedo adjustment dependent on SWR
  adfl,              &! Snow albedo adjustment dependent on LWR
  fsar,              &! Snow albedo adjustment range dependent on vegetation fraction
  psf,               &! Scaling factor for solid precipitation (within forest stand, at min CC)
  psr,               &! Range of solid precipitation (within forest stand, spread min-max CC)
  wcan,              &! Parameter of exponential wind profile 
  zsub,              &! Sub-canopy reference height (m)
  zgf,               &! Roughness length adjustment factor depending on vegetation fraction  
  zgr,               &! Roughness length adjustment range depending on vegetation fraction  
  khcf                ! Diffusivity adjustment for canopy effects (Finnigan 2000)
end module PARAMETERS

!-----------------------------------------------------------------------
! Spatial surface characteristics
!-----------------------------------------------------------------------
module PARAMMAPS
real, allocatable :: &
  adc(:,:),          &! Cold snow albedo decay time (h)
  afs(:,:),          &! Maximum albedo for fresh snow
  alb0(:,:),         &! Snow-free ground albedo
  canh(:,:),         &! Canopy heat capacity (J/K/m^2)
  fcly(:,:),         &! Soil clay fraction
  fsnd(:,:),         &! Soil sand fraction
  fsky(:,:),         &! Sky view fraction
  scap(:,:),         &! Canopy snow capacity (kg/m^2)
  trcn(:,:),         &! Canopy transmissivity
  VAI(:,:),          &! Vegetation area index
  veg_shd(:,:),      &! Vegetation snow holding depth (m)
  z0sf(:,:),         &! Snow-free roughness length (m)
  z0_snow(:,:)        ! Roughness length of snow (m)
end module PARAMMAPS

!-----------------------------------------------------------------------
! Soil properties
!-----------------------------------------------------------------------
module SOILPARAMS
real, allocatable :: &
  b(:,:),            &! Clapp-Hornberger exponent
  hcap_soil(:,:),    &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil(:,:),    &! Thermal conductivity of dry soil (W/m/K)
  sathh(:,:),        &! Saturated soil water pressure (m)
  Vcrit(:,:),        &! Volumetric soil moisture at critical point
  Vsat(:,:)           ! Volumetric soil moisture at saturation
end module SOILPARAMS

!-----------------------------------------------------------------------
! Model state variables  
!-----------------------------------------------------------------------
module STATE_VARIABLES
! Canopy properties
real, allocatable :: &
  Qcan(:,:),         &! Canopy air space humidity
  Tcan(:,:),         &! Canopy air space temperature (K)
  Sveg(:,:),         &! Snow mass on vegetation (kg/m^2)
  Tveg(:,:)           ! Vegetation temperature (K)
  
! Surface state variables
real, allocatable :: &
  Tsrf(:,:),         &! Surface skin temperature (K)
  fsnow(:,:)          ! Snow cover fraction terrain

! Snow state variables
integer, allocatable :: &
  Nsnow(:,:)          ! Number of snow layers
real, allocatable ::     &
  albs(:,:),             &! Snow albedo
  Ds(:,:,:),             &! Snow layer thicknesses (m)
  histowet(:,:,:),       &! Historical variable for past wetting of a layer (0-1)
  Sice(:,:,:),           &! Ice content of snow layers (kg/m^2)
  Sliq(:,:,:),           &! Liquid content of snow layers (kg/m^2)
  Tsnow(:,:,:),          &! Snow layer temperatures (K)
  swehist(:,:,:),        &! History of SWE during last 14 days (kg/m^2). Most recent entries first.
  swemin(:,:),           &! Minimum SWE during the season (kg/m^2)
  swemax(:,:),           &! Maximum SWE during the season (kg/m^2)
  snowdepthhist(:,:,:),  &! History of snow depth during last 14 days (m). Most recent entries first.
  snowdepthmin(:,:),     &! Minimum snow depth at time step of swemin (m)
  snowdepthmax(:,:)       ! Maximum snow depth at time step of swemax (m)

! Soil state variables
real, allocatable :: &
  theta(:,:,:),          &! Volumetric moisture content of soil layers
  Tsoil(:,:,:)            ! Soil layer temperatures (K)
  
! SnowTran3D state variables
real, allocatable ::     &
  dSWE_tot_subl(:,:),    &! Cumulated SWE change due to snowdrift sublimation (kg/m^2)
  dSWE_tot_salt(:,:),    &! Cumulated SWE change due to saltation (kg/m^2)
  dSWE_tot_susp(:,:)      ! Cumulated SWE change due to suspension (kg/m^2)
  
! SnowSlide state variables
real, allocatable ::     &
  dSWE_tot_slide(:,:)     ! Cumulated SWE change due to snow slides (kg/m^2)
  
integer, allocatable ::  &
  index_sorted_dem(:,:)   ! Location (i,j) of sorted grid points
  
end module STATE_VARIABLES

!-----------------------------------------------------------------------
! Landuse information 
!-----------------------------------------------------------------------
module LANDUSE
! Canopy properties
real, allocatable :: &
  fveg(:,:),         &! Canopy cover fraction
  fves(:,:),         &! Stand-scale canopy cover fraction
  lai(:,:),          &! Leaf area index 
  vfhp(:,:),         &! Hemispherical sky-view fraction including canopy
  hcan(:,:)           ! Canopy height (m)

! Terrain properties
real, allocatable :: &
  fsky_terr(:,:),    &! Sky view fraction terrain
  slopemu(:,:),      &! slope parameter 
  xi(:,:),           &! terrain correlation length
  cellsize(:,:),     &! grid cell size (m)
  dem(:,:),          &! grid elevation (m)
  dem_sorted(:),     &! sorted grid elevation (m)
  slope(:,:),        &! slope (deg)
  Shd(:,:),          &! snow holding depth, for SnowSlide (m)
  pmultf(:,:)         ! precip multiplier to revert precip correction applied to open area

! Tile properties 
real, allocatable :: &
  tilefrac(:,:),     &! Tile fraction 
  glacierfrac(:,:),  &! Glacier flag
  forestfrac(:,:)     ! Forest fraction

end module LANDUSE

!-----------------------------------------------------------------------
! Parameters for SnowTran3D
!-----------------------------------------------------------------------
module PARAM_SNOWTRAN3D

logical, parameter :: &

  ! Define whether the threshold surface shear velocity will be
  ! constant during the simulation (flag_variable_Utau_t = .FALSE.), or whether
  ! it will evolve as a funtion of air temperature and wind speed
  ! (flag_variable_Utau_t = .TRUE.).
  !!! (recommended default value: flag_variable_Utau_t = .TRUE.)
  flag_variable_Utau_t = .TRUE.,   &! Flag for variable threshold friction velocitiy (.TRUE.) or constant (.FALSE.)

  
  ! Define whether you want the upwind boundaries to have zero
  ! incoming transport flux, or the equilibrium transport flux.
  ! flag_boundary_condition = .FALSE. corresponds to a zero flux, and flag_boundary_condition = .TRUE.
  ! corresponds to the equilibrium flux.
  !!! (recommended default value: flag_boundary_condition = .FALSE.)
  flag_boundary_condition = .FALSE. ! Boundary condition flag

real, parameter :: &

  ! For the case of flag_variable_Utau_t = 0.0, define what the threshold
  ! surface shear velocity (m/s) will be during the simulation.
  !!! (recommended default value: Utau_t_const = 0.25)
  Utau_t_const = 0.25,             &! Constant threshold friction velocity (m/s) used if flag_variable_Utau_t = 0.0

  ! Density of snow deposited by snowdrift
  !!! (recommended default value: rho_snow = 300.0)
  rho_snow = 300.0,                &! Constant snow density (kg/m^3)

  ! The blowby parameter is implemented to account for the erosion
  ! of the tops of deep snow accumulations.  It corrects a
  ! deficiency in the du*/dx* < 0 formulation.  It is a number that
  ! should range from 0 to 1.0, and represents the fraction of the
  ! upwind saltation flux that is transfered farther downwind into
  ! the next grid cell.  So, the bigger the number, the less
  ! peaked the drift accumulation profile is.  blowby = 0.0 is the
  ! original model.  The Tabler surfaces can be used to do the
  ! same kind of thing.
  !!! (recommended default value: blowby = 0.0)
  blowby = 0.75                     ! Fraction of the saltation flux transferred downwind
  
end module PARAM_SNOWTRAN3D

!-----------------------------------------------------------------------
! Constants for SnowTran3D
!-----------------------------------------------------------------------
module CONSTANTS_SNOWTRAN3D
real, parameter :: &

  ! Constants related to surface shear stress and saltation transport.
  fetch = 500.0,        &! Equilibrium fetch distance (m)
  xmu = 3.0,            &! Scaling constant for non-equilibrium saltation transport
  C_z = 0.12,           &! Coefficient 0.12 in Liston and Sturm (1998) eq. 5 p. 500
  h_const = 1.6,        &! Coefficient 1.6 in Liston and Sturm (1998) eq. 14 p. 501
  wind_min = 4.0,       &! Minimum wind speed to compute snow transport (m/s)
  
  ! Constants related to suspended snow profile.
  Up_const = 2.8,       &! Constant coefficient for calculation of U_p
  dz_susp = 0.20,       &! dz in suspension layer (m)
  ztop_susp = 2.0,      &! Height of the top of suspension layer (m)
  fall_vel = 0.3,       &! Particle-settling velocity, s in L&S 1998 eq. 12 (m/s)
  Ur_const = 0.5         ! Constant coefficient beta in L&S 1998, eq. 13
  
end module CONSTANTS_SNOWTRAN3D

!-----------------------------------------------------------------------
! Parameters for SnowSlide
!-----------------------------------------------------------------------
module PARAM_SNOWSLIDE
real, parameter :: &

  ! Dynamic ratio for the snow holding depth of pixels receiving an avalanche
  !!! (recommended default value: dyn_ratio = 0.7)
  dyn_ratio = 0.7,  &! Dynamic snow holding depth ratio

  ! Hysteretic triggering ratio for the snow holding depth to trigger an avalanche
  !!! (recommended default value: trig_ratio = 1.0)
  trig_ratio = 1.5,  &! Hysteretic triggering snow holding depth ratio

  ! Density of snow deposits by avalanches
  !!! (recommended default value: rho_snow = 300.0)
  rho_deposit = 300.0,  &! Constant snow avalanche deposit density (kg/m^3)  

  ! Snow slide can only be triggered above a certain slope
  !!! (recommended default value: slope_min = 25 deg)
  slope_min = 25.0,     &! Minimum slope for snow slide occurrence (deg)
  
  ! Snow holding depth has a min threshold 
  !!! (recommended default value: Shd_min = 0.05 m)
  Shd_min = 0.05         ! Minimum snow holding depth (m)
  
end module PARAM_SNOWSLIDE
