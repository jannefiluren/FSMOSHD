!-----------------------------------------------------------------------
! Physical constants
!-----------------------------------------------------------------------
module CONSTANTS

use, intrinsic :: iso_fortran_env, only: dp=>real64 

real*8, parameter :: &
  cp = 1005_dp,         &! Specific heat capacity of air (J/K/kg)
  eps = 0.622_dp,       &! Ratio of molecular weights of water and dry air
  e0 = 610.78_dp,       &! Saturation vapour pressure at Tm (Pa)
  grav = 9.81_dp,       &! Acceleration due to gravity (m/s^2)
  hcap_ice = 2100_dp,   &! Specific heat capacity of ice (J/K/kg)
  hcap_wat = 4180_dp,   &! Specific heat capacity of water (J/K/kg)
  hcon_air = 0.025_dp,  &! Thermal conductivity of air (W/m/K)
  hcon_clay = 1.16_dp,  &! Thermal conductivity of clay (W/m/K)
  hcon_ice = 2.24_dp,   &! Thermal conducivity of ice (W/m/K)
  hcon_sand = 1.57_dp,  &! Thermal conductivity of sand (W/m/K)
  hcon_wat = 0.56_dp,   &! Thermal conductivity of water (W/m/K)
  I0 = 1367_dp,         &! Solar constant (W/m^2)
  Lf = 0.334e6_dp,      &! Latent heat of fusion (J/kg)
  Lv = 2.501e6_dp,      &! Latent heat of vapourisation (J/kg)
  Ls = Lf + Lv,         &! Latent heat of sublimation (J/kg)
  pi = 3.14159_dp,      &! pi
  Rair = 287_dp,        &! Gas constant for air (J/K/kg)
  Rwat = 462_dp,        &! Gas constant for water vapour (J/K/kg)
  rho_ice = 917_dp,     &! Density of ice (kg/m^3)
  rho_wat = 1000_dp,    &! Density of water (kg/m^3)
  sb = 5.67e-8_dp,      &! Stefan-Boltzmann constant (W/m^2/K^4)
  em_snow = 0.99_dp,    &! Emissivity snow for Stefan-Boltzmann
  em_soil = 0.90_dp,    &! Emissivity soil for Stefan-Boltzmann
  Tm = 273.15_dp,       &! Melting point (K)
  vkman = 0.4_dp,       &! Von Karman constant
  undef =  1.e+6_dp      ! Initialization value for allocatables
integer, parameter :: &
  iundef = 1.e+6      ! Initialization value for integer allocatables
end module CONSTANTS

!-----------------------------------------------------------------------
! Model configuration
!-----------------------------------------------------------------------
module MODCONF
! Process options                            : Possible values
integer :: &
  ALBEDO,     &! snow albedo                 : 0, 1, 2
  CANMOD,     &! forest canopy               : 0, 1
  CONDCT,     &! snow thermal conductivity   : 0, 1
  DENSTY,     &! snow density                : 0, 1, 2, 3
  EXCHNG,     &! turbulent exchange          : 0, 1
  HYDROL,     &! snow hydraulics             : 0, 1, 2
  SNFRAC,     &! snow cover fraction         : 0, 1, 2, 3, 4
  RADSBG,     &! subgrid radiation param     : 0, 1
! Driving data options
  ZOFFST,     &! measurement height offset    : 0, 1
! OSHD-specific options
  OSHDTN       ! oshd-specific tuning options of fresh snow albedo, snow roughness lengths and fresh snow density: 0, 1
logical :: &
  HN_ON,      &!activate the new snow model
  FOR_HN       !write 18h states for the hn model.
end module MODCONF

module MODPERT
logical:: &
  Z0PERT,     &! activate z0 perturbations
  WCPERT,     &! activate liquid water capacity perturbations
  FSPERT,     &! activate fresh snow density perturbations
  ALPERT,     &! activate albedo perturbations
  SLPERT       ! activate settling perturbations
end module MODPERT

!-----------------------------------------------------------------------
! Model tile
!-----------------------------------------------------------------------
module MODTILE
! Model tiles                                   : Possible values
character(len=20) :: &
  TILE           ! model tile                   : 'open', 'forest', 'glacier' 
real*8 :: & 
  tthresh        ! Tile fraction of grid cell required for tile to be considered 
end module MODTILE 

!-----------------------------------------------------------------------
! Model output configuration
!-----------------------------------------------------------------------
module MODOUTPUT
! list of diagnostic variables that can be written to output bin files.
! at the moment, only 2d real variables are handled.
character(len=4), dimension(35) :: &
  WRITE_DIAG_VARS = (/'rotc', &  ! Roff      Total runoff, snow and bare soil (kg/m^2)
                      'hsnt', &  ! snowdepth Total snowdepth (m)
                      'swet', &  ! SWE       Total SWE (kg/m^2)
                      'slqt', &  ! Sliq_out  Total LWC (kg/m^2)
                      'swtb', &  ! Sdirt     Incoming direct beam radiation corrected for subgrid topography (W/m^2)
                      'swtd', &  ! Sdift     Incoming diffuse beam radiation corrected for subgrid topography (W/m^2)
                      'lwtr', &  ! Lwt       Incoming longwave radiation corrected for subgrid topography (W/m^2)
                      'romc', &  ! meltflux_out Runoff from snowmelt at the base of snow (kg/m^2)
                      'sbsc', &  ! Sbsrf     Snow sublimation rate (kg/m^2/s)
                      'asrf', &  ! asrf_out  Surface albedo
                      'emlt', &  ! Melt      Surface melt rate (kg/m^2/s)
                      'esrf', &  ! Esrf      Moisture flux from the surface (kg/m^2/s)
                      'eveg', &  ! Eveg      Moisture flux from vegetation (kg/m^2/s)
                      'ghsl', &  ! Gsoil     Heat flux into soil (W/m^2)
                      'hesr', &  ! Hsrf      Sensible heat flux from the surface (W/m^2)
                      'intc', &  ! intcpt    Canopy interception (kg/m^2)
                      'khag', &  ! KH        Eddy diffusivity for heat to the atmosphere (m/s)
                      'khac', &  ! KHa       Eddy diffusivity for heat from the canopy air space (m/s)
                      'khgr', &  ! KHg       Eddy diffusivity for heat from the ground (m/s)
                      'khve', &  ! KHv       Eddy diffusivity for heat from vegetation (m/s)
                      'kwgr', &  ! KWg       Eddy diffusivity for water from the ground (m/s)
                      'kwve', &  ! KWv       Eddy diffusivity for water from vegetation (m/s)
                      'lahe', &  ! LE        Latent heat flux to the atmosphere (W/m^2)
                      'lesr', &  ! LEsrf     Latent heat flux from the surface (W/m^2)
                      'lwsc', &  ! LWsci     Subcanopy incoming longwave radiation (W/m^2)
                      'lwve', &  ! LWveg     Net longwave radiation absorbed by vegetation (W/m^2)
                      'rnet', &  ! Rnet      Net radiation (W/m^2)
                      'rnsr', &  ! Rsrf      Net radiation at surface (W/m^2)
                      'sbve', &  ! Sbveg     Sublimation from vegetation (kg/m^2)
                      'sehe', &  ! H         Sensible heat flux to the atmosphere (W/m^2)
                      'swsc', &  ! SWsci     Subcanopy incoming shortwave radiation (W/m^2)
                      'swsr', &  ! SWsrf     Net SW radiation absorbed by the surface (W/m^2)
                      'swve', &  ! SWveg     Net SW radiation absorbed by vegetation (W/m^2)
                      'uasc', &  ! Usc       Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)
                      'unld'  &  ! unload    Snow mass unloaded from canopy (kg/m^2)
                      /)
   ! GM Note: Energy-balance relevant diagnostics omitted in this list can be derived from existing variables:
                    !Hveg = H-Hsrf
                    !LEveg = LE-LEsrf
                    !LWsrf = Tss^4*sb
                    !Rveg  = SWveg + LWveg      
! list of state variables that can be written to output bin files.
! at the moment, only 2d real variables are handled.
character(len=4), dimension(11):: &
  WRITE_STATE_VARS = (/&
                  'alse', &  ! albs
                ! 'hsnl', &  ! Ds(3d var)
                  'scfe', &  ! fsnow
                ! 'nsne', &  ! Nsnow
                ! 'sicl', &  ! Sice (3d var)
                ! 'slql', &  ! Sliq (3d var)
                  'tsfe', &  ! Tsrf
                ! 'tsnl', &  ! Tsnow (3d var)
                ! 'tsll', &  ! Tsoil (3d var)
                  'hsmn', &  ! snowdepthmin
                  'hsmx', &  ! snowdepthmax
                ! 'hshs', &  ! snowdepthhist (3d var, 1st dimension: time.)
                  'swmn', &  ! swemin
                  'swmx', &  ! swemax
                ! 'swhs', &  ! swehist (3d var, 1st dimension: time.)
                  'qcan', &  ! Qcan
                  'sveg', &  ! Sveg
                  'tcan', &  ! Tcan
                  'tveg'  &  ! Tveg
                  /)
character(len=4), dimension(35) :: &
  LIST_DIAG_RESULTS  ! List of diagnostic variables that the user wants to write into output bin files (subset of WRITE_DIAG_VARS)
character(len=4), dimension(11) :: &
  LIST_STATE_RESULTS ! List of result variables that the user wants to write into output bin files (subset of WRITE_STATE_VARS)
logical, dimension(35)  :: WRITE_DIAG_TABLE  ! table specifying whether the corresponding variables in WRITE_DIAG_VARS should be written or not.
logical, dimension(11) :: WRITE_STATE_TABLE ! table specifying whether the corresponding variables in WRITE_STATE_VARS should be written or not.
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
real*8 :: &
  hour                ! Hour of day
real*8 :: &
  dt,                &! Timestep (s)
  zT,                &! Temperature measurement height (m)
  zU                  ! Wind speed measurement height (m)
real*8, allocatable :: &
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
  z0P(:,:),          &! z0 perturbations
  wcP(:,:),          &! liquid water capacity perturbations
  fsP(:,:),          &! fresh snow density perturbations
  alP(:,:),          &! albedo perturbations
  slP(:,:)            ! settling perturbations
end module DRIVING

!-----------------------------------------------------------------------
! Grid parameters
!-----------------------------------------------------------------------
module GRID
integer :: &
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions
real*8, allocatable :: &
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
real*8 :: &
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
real*8 :: &
  a_eta,             &! Temperature factor for Crocus B92 compaction (K^-1)
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
  rgr0,              &! Fresh snow grain radius (m)
  rho0,              &! Fixed snow density (kg/m^3)
  rhob,              &! Temperature factor in fresh snow density (kg/m^3/K)
  rhoc,              &! Wind factor in fresh snow density (kg s^0.5/m^3.5)
  rhof,              &! Fresh snow density (kg/m^3)
  rcld,              &! Maximum density for cold snow (kg/m^3)
  rmlt,              &! Maximum density for melting snow (kg/m^3)
  snda,              &! Thermal metamorphism parameter (1/s)
  Talb,              &! Albedo decay temperature threshold (C)
  tcld,              &! Cold snow albedo decay time scale (s)
  tmlt,              &! Melting snow albedo decay time scale (s)
  trho,              &! Snow compaction time scale (s)
  Wirr,              &! Irreducible liquid water content of snow
  z0sn,              &! Snow roughness length (m)
  Sfmin               ! Minimum 24h snowfall to refresh albedo (kg/m^2)

! Surface parameters
real*8 :: &
  gsat                ! Surface conductance for saturated soil (m/s)
! Additional parameters used in forest snow process parametrizations
real*8 :: &      
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
real*8, allocatable :: &
  alb0(:,:),         &! Snow-free ground albedo
  canh(:,:),         &! Canopy heat capacity (J/K/m^2)
  fcly(:,:),         &! Soil clay fraction
  fsnd(:,:),         &! Soil sand fraction
  fsky(:,:),         &! Sky view fraction
  scap(:,:),         &! Canopy snow capacity (kg/m^2)
  trcn(:,:),         &! Canopy transmissivity
  VAI(:,:),          &! Vegetation area index
  z0sf(:,:)           ! Snow-free roughness length (m)
end module PARAMMAPS

!-----------------------------------------------------------------------
! Soil properties
!-----------------------------------------------------------------------
module SOILPARAMS
real*8, allocatable :: &
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
real*8, allocatable :: &
  Qcan(:,:),         &! Canopy air space humidity
  Tcan(:,:),         &! Canopy air space temperature (K)
  Sveg(:,:),         &! Snow mass on vegetation (kg/m^2)
  Tveg(:,:)           ! Vegetation temperature (K)
  
! Surface state variables
real*8, allocatable :: &
  Tsrf(:,:),         &! Surface skin temperature (K)
  fsnow(:,:)          ! Snow cover fraction terrain

! Snow state variables
integer, allocatable :: &
  Nsnow(:,:)          ! Number of snow layers
real*8, allocatable ::     &
  albs(:,:),             &! Snow albedo
  Ds(:,:,:),             &! Snow layer thicknesses (m)
  rgrn(:,:,:),           &! Snow layer grain radius (m)
  Sice(:,:,:),           &! Ice content of snow layers (kg/m^2)
  Sliq(:,:,:),           &! Liquid content of snow layers (kg/m^2)
  Tsnow(:,:,:),          &! Snow layer temperatures (K)
  swehist(:,:,:),        &! history of SWE during last 14 days (kg/m^2). Most recent entries first.
  swemin(:,:),           &! Minimum swe during the season (m)
  swemax(:,:),           &! Maximum swe during the season (m)
  snowdepthhist(:,:,:),  &! history of Snow depth during last 14 days (m). Most recent entries first.
  snowdepthmin(:,:),     &! Minimum Snow depth at time step of swemin (m)
  snowdepthmax(:,:)       ! Maximum Snow depth at time stemp of swemax(m)

! Soil state variables
real*8, allocatable :: &
  theta(:,:,:),      &! Volumetric moisture content of soil layers
  Tsoil(:,:,:)        ! Soil layer temperatures (K)
end module STATE_VARIABLES

!-----------------------------------------------------------------------
! Landuse information 
!-----------------------------------------------------------------------
module LANDUSE
! Canopy properties
real*8, allocatable :: &
  fveg(:,:),         &! Canopy cover fraction
  fves(:,:),         &! Stand-scale canopy cover fraction
  lai(:,:),          &! Leaf area index 
  vfhp(:,:),         &! Hemispherical sky-view fraction including canopy
  hcan(:,:)           ! Canopy height (m)

! Terrain properties
real*8, allocatable :: &
  fsky_terr(:,:),    &! Sky view fraction terrain
  slopemu(:,:),      &! slope parameter 
  xi(:,:),           &! terrain correlation length
  Ld(:,:),           &! grid cell size or domain size (m)
  lat(:,:),          &! latitude of each grid cell (center?)
  lon(:,:),          &! longitude of each grid cell (center?) 
  dem(:,:),          &! grid elevation
  pmultf(:,:)          ! precip multiplier to revert precip correction applied to open area

! Tile properties 
real*8, allocatable :: &
  tilefrac(:,:)       ! Tile fraction 
end module LANDUSE

!-----------------------------------------------------------------------
! Input / output unit numbers
!-----------------------------------------------------------------------
module IOUNITS
  integer, parameter :: &
    udmp = 11,         &! Dump file unit number
    umet = 21,         &! Driving file unit number
    uout = 31,         &! Output file unit number
    ustr = 41           ! Start file unit number
end module IOUNITS

module FILES
  character(len=200) :: &
    dump_file           ! Dump file name
end module FILES

module TEST
  integer :: &
    dump_data           ! Write debug info to file
end module TEST
