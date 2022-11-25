@with_kw mutable struct FSM{T}
  
  # Maximum snow and soil layer thicknesses
  Dzsnow = [0.1, 0.2, 0.4]
  Dzsoil = [0.1, 0.2, 0.4, 0.8]
  
  # Base variables
  Nx::Int = 1
  Ny::Int = 1
  Nsmax::Int = length(Dzsnow)
  Nsoil::Int = length(Dzsoil)

  # Driving data
  dt::T = 3600*1
  zT::T = 10
  zU::T = 10

  # Model configuration
  ALBEDO::Int = 2
  CANMOD::Int = 0
  CONDCT::Int = 1
  DENSTY::Int = 3
  EXCHNG::Int = 1
  HYDROL::Int = 2
  SNFRAC::Int = 3
  RADSBG::Int = 0
  ZOFFST::Int = 0
  OSHDTN::Int = 1
  HN_ON = false    # TODO add type
  FOR_HN = true    # TODO add type

  tthresh::T = 0.1    # TODO check if this is the best handling
  TILE = "open"    # TODO check if this is the best handling and type

  # Numerical solution parameters
  Nitr = 4

  # Defaults for canopy parameters
  avg0::T = 0.1
  avgs::T = 0.4
  cden::T = 0.004
  cvai::T = 6.6
  cveg::T = 20
  Gcn1::T = 0.5
  Gcn2::T = 0
  gsnf::T = 0
  kdif::T = 0.5
  kveg::T = 1
  rchd::T = 0.67
  rchz::T = 0.2          
  tcnc::T = 3600*240
  tcnm::T = 3600*48

  # Defaults for snow parameters
  a_eta::T = 0.1
  asmx::T = 0.8
  asmn::T = 0.5
  b_eta::T = 0.023
  bthr::T = 2
  c_eta::T = 250
  eta0::T = 3.7e7
  eta1::T = 7.62237e6
  hfsn::T = 0.1 
  kfix::T = 0.24
  rho0::T = 300
  rhob::T = 6
  rhoc::T = 26
  rhof::T = 109
  rcld::T = 300
  rgr0::T = 5e-5
  rmlt::T = 500
  Salb::T = 10
  snda::T = 2.8e-6
  Talb::T = -2
  tcld::T = 3600*1000
  tmlt::T = 3600*100
  trho::T = 3600*200
  Wirr::T = 0.03
  z0sn::T = 0.002
  Sfmin::T = 10

  # Defaults for ground surface parameters
  bstb::T = 5
  gsat::T = 0.01

  # Defaults for additional parameters required for forest snow process parametrization
  adfs::T = 3
  adfl::T = 2
  fsar::T = 0.1
  psf::T  = 1
  psr::T  = 0.1
  wcan::T = 2.5
  zsub::T = 2
  zgf::T = 5
  zgr::T = 5
  khcf::T = 3

  # Surface properties
  alb0::Array{T} = 0.2*ones(Nx,Ny)
  z0sf::Array{T} = 0.2*ones(Nx,Ny)
  fcly::Array{T} = 0.3*ones(Nx,Ny)
  fsnd::Array{T} = 0.6*ones(Nx,Ny)

  # Canopy parameters
  canh::Array{T} = zeros(Nx,Ny)   # todo may need to change
  fsky::Array{T} = zeros(Nx,Ny)   # todo may need to change
  fveg::Array{T} = zeros(Nx,Ny)   # todo may need to change
  fves::Array{T} = zeros(Nx,Ny)   # todo may need to change
  hcan::Array{T} = zeros(Nx,Ny)   # todo may need to change
  lai::Array{T} = zeros(Nx,Ny)   # todo may need to change
  pmultf::Array{T} = zeros(Nx,Ny)   # todo may need to change
  scap::Array{T} = zeros(Nx,Ny)   # todo may need to change
  trcn::Array{T} = zeros(Nx,Ny)   # todo may need to change
  VAI::Array{T} = zeros(Nx,Ny)   # todo may need to change
  vfhp::Array{T} = zeros(Nx,Ny)   # todo may need to change

  # Terrain properties
  slopemu::Array{T} = zeros(Nx,Ny)
  xi::Array{T} = zeros(Nx,Ny)
  Ld::Array{T} = zeros(Nx,Ny)
  lat::Array{T} = zeros(Nx,Ny)
  lon::Array{T} = zeros(Nx,Ny)
  dem::Array{T} = zeros(Nx,Ny)
  tilefrac::Array{T} = zeros(Nx,Ny)

  # Derived soil parameters
  b::Array{T} = zeros(Nx,Ny)
  hcap_soil::Array{T} = zeros(Nx,Ny)
  hcon_soil::Array{T} = zeros(Nx,Ny)
  sathh::Array{T} = zeros(Nx,Ny)
  Vsat::Array{T} = zeros(Nx,Ny)
  Vcrit::Array{T} = zeros(Nx,Ny)

  # State variables
  albs::Array{T} = 0.85.*ones(Nx,Ny)
  Ds::Array{T} = zeros(Nsmax,Nx,Ny)
  Nsnow::Array{Int} = zeros(Nx,Ny)
  Qcan::Array{T} = zeros(Nx,Ny)
  rgrn::Array{T} = zeros(Nsmax,Nx,Ny)
  Sice::Array{T} = zeros(Nsmax,Nx,Ny)
  Sliq::Array{T} = zeros(Nsmax,Nx,Ny)
  Sveg::Array{T} = zeros(Nx,Ny)
  Tcan::Array{T} = 273.15.*ones(Nx,Ny)
  theta::Array{T} = zeros(Nsoil,Nx,Ny)
  Tsnow::Array{T} = 273.15.*ones(Nsmax,Nx,Ny)
  Tsoil::Array{T} = 273.15.*ones(Nsoil,Nx,Ny)
  Tsrf::Array{T} = zeros(Nx,Ny)
  fsnow::Array{T} = zeros(Nx,Ny)
  Tveg::Array{T} = 273.15.*ones(Nx,Ny)
  snowdepthmin::Array{T} = zeros(Nx,Ny)
  snowdepthmax::Array{T} = zeros(Nx,Ny)
  snowdepthhist::Array{T} = zeros(14,Nx,Ny)
  swemin::Array{T} = zeros(Nx,Ny)
  swemax::Array{T} = zeros(Nx,Ny)
  swehist::Array{T} = zeros(14,Nx,Ny)
  fsky_terr::Array{T} = zeros(Nx,Ny)

  # Radiation - temporary arrays

  alb::Array{T} = zeros(Nx,Ny)
  asrf_out::Array{T} = zeros(Nx,Ny)
  Sdirt::Array{T} = zeros(Nx,Ny)
  Sdift::Array{T} = zeros(Nx,Ny)
  SWveg::Array{T} = zeros(Nx,Ny)
  SWsrf::Array{T} = zeros(Nx,Ny)
  SWsci::Array{T} = zeros(Nx,Ny)
  LWt::Array{T} = zeros(Nx,Ny)
  SWtopo_out::Array{T} = zeros(Nx,Ny)

  # Thermal - temporary arrays

  ksnow::Array{T} = zeros(Nsmax, Nx, Ny)
  csoil::Array{T} = zeros(Nsoil, Nx, Ny)
  ksoil::Array{T} = zeros(Nsoil, Nx, Ny)
  gs1::Array{T} = zeros(Nx, Ny)
  Ds1::Array{T} = zeros(Nx, Ny)
  Ts1::Array{T} = zeros(Nx, Ny)
  ks1::Array{T} = zeros(Nx, Ny)
  Tveg0::Array{T} = zeros(Nx, Ny)

  # Sfexch - temporary arrays

  KH::Array{T} = zeros(Nx, Ny)
  KHa::Array{T} = zeros(Nx, Ny)
  KHg::Array{T} = zeros(Nx, Ny)
  KHv::Array{T} = zeros(Nx, Ny)
  KWg::Array{T} = zeros(Nx, Ny)
  KWv::Array{T} = zeros(Nx, Ny)
  Usc::Array{T} = zeros(Nx, Ny)

  # Ebalsrf - temporary arrays

  Esrf::Array{T} = zeros(Nx,Ny)
  Eveg::Array{T} = zeros(Nx,Ny)
  G::Array{T} = zeros(Nx,Ny)
  H::Array{T} = zeros(Nx,Ny)
  Hsrf::Array{T} = zeros(Nx,Ny)
  LE::Array{T} = zeros(Nx,Ny)
  LEsrf::Array{T} = zeros(Nx,Ny)
  LWsci::Array{T} = zeros(Nx,Ny)
  LWveg::Array{T} = zeros(Nx,Ny)
  Melt::Array{T} = zeros(Nx,Ny)
  Rnet::Array{T} = zeros(Nx,Ny)
  Rsrf::Array{T} = zeros(Nx,Ny)

  # Snow - temporary arrays

  Gsoil::Array{T} = zeros(Nx, Ny)
  Roff::Array{T} = zeros(Nx, Ny)
  meltflux_out::Array{T} = zeros(Nx, Ny)
  Sbsrf::Array{T} = zeros(Nx, Ny)
  Roff_bare::Array{T} = zeros(Nx, Ny)
  Roff_snow::Array{T} = zeros(Nx, Ny)
  fsnow_thres::Array{T} = zeros(Nx, Ny)
  unload::Array{T} = zeros(Nx, Ny)

end