@with_kw mutable struct FSM{T}
  
  # Maximum snow and soil layer thicknesses
  Dzsnow::Vector{T} = [0.1, 0.2, 0.4]
  Dzsoil::Vector{T} = [0.1, 0.2, 0.4, 0.8]
  
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
  HN_ON::Bool = false
  FOR_HN::Bool = true

  tthresh::T = 0.1    # TODO check if this is the best handling
  TILE::String = "open"    # TODO check if this is the best handling and type

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
  alb0::Array{T,2} = 0.2*ones(Nx,Ny)
  z0sf::Array{T, 2} = 0.2*ones(Nx,Ny)
  fcly::Array{T, 2} = 0.3*ones(Nx,Ny)
  fsnd::Array{T, 2} = 0.6*ones(Nx,Ny)

  # Canopy parameters
  canh::Array{T, 2} = zeros(Nx,Ny)   # todo may need to change
  fsky::Array{T, 2} = zeros(Nx,Ny)   # todo may need to change
  fveg::Array{T, 2} = zeros(Nx,Ny)   # todo may need to change
  fves::Array{T, 2} = zeros(Nx,Ny)   # todo may need to change
  hcan::Array{T, 2} = zeros(Nx,Ny)   # todo may need to change
  lai::Array{T, 2} = zeros(Nx,Ny)   # todo may need to change
  pmultf::Array{T, 2} = zeros(Nx,Ny)   # todo may need to change
  scap::Array{T, 2} = zeros(Nx,Ny)   # todo may need to change
  trcn::Array{T, 2} = zeros(Nx,Ny)   # todo may need to change
  VAI::Array{T, 2} = zeros(Nx,Ny)   # todo may need to change
  vfhp::Array{T, 2} = zeros(Nx,Ny)   # todo may need to change

  # Terrain properties
  slopemu::Array{T, 2} = zeros(Nx,Ny)
  xi::Array{T, 2} = zeros(Nx,Ny)
  Ld::Array{T, 2} = zeros(Nx,Ny)
  lat::Array{T, 2} = zeros(Nx,Ny)
  lon::Array{T, 2} = zeros(Nx,Ny)
  dem::Array{T, 2} = zeros(Nx,Ny)
  tilefrac::Array{T,2} = zeros(Nx,Ny)

  # Derived soil parameters
  b::Array{T, 2} = zeros(Nx,Ny)
  hcap_soil::Array{T, 2} = zeros(Nx,Ny)
  hcon_soil::Array{T, 2} = zeros(Nx,Ny)
  sathh::Array{T, 2} = zeros(Nx,Ny)
  Vsat::Array{T, 2} = zeros(Nx,Ny)
  Vcrit::Array{T, 2} = zeros(Nx,Ny)

  # State variables
  albs::Array{T, 2} = 0.85.*ones(Nx,Ny)
  Ds::Array{T,3} = zeros(Nsmax,Nx,Ny)
  Nsnow::Array{Int,2} = zeros(Nx,Ny)
  Qcan::Array{T, 2} = zeros(Nx,Ny)
  rgrn::Array{T,3} = zeros(Nsmax,Nx,Ny)
  Sice::Array{T,3} = zeros(Nsmax,Nx,Ny)
  Sliq::Array{T,3} = zeros(Nsmax,Nx,Ny)
  Sveg::Array{T, 2} = zeros(Nx,Ny)
  Tcan::Array{T, 2} = 273.15.*ones(Nx,Ny)
  theta::Array{T,3} = zeros(Nsoil,Nx,Ny)
  Tsnow::Array{T,3} = 273.15.*ones(Nsmax,Nx,Ny)
  Tsoil::Array{T,3} = 273.15.*ones(Nsoil,Nx,Ny)
  Tsrf::Array{T, 2} = zeros(Nx,Ny)
  fsnow::Array{T, 2} = zeros(Nx,Ny)
  Tveg::Array{T, 2} = 273.15.*ones(Nx,Ny)
  snowdepthmin::Array{T, 2} = zeros(Nx,Ny)
  snowdepthmax::Array{T, 2} = zeros(Nx,Ny)
  snowdepthhist::Array{T,3} = zeros(14,Nx,Ny)
  swemin::Array{T, 2} = zeros(Nx,Ny)
  swemax::Array{T, 2} = zeros(Nx,Ny)
  swehist::Array{T,3} = zeros(14,Nx,Ny)
  fsky_terr::Array{T, 2} = zeros(Nx,Ny)

  # Radiation - temporary arrays

  afs::Array{T,2} = zeros(Nx, Ny)
  alb::Array{T,2} = zeros(Nx,Ny)
  asrf_out::Array{T,2} = zeros(Nx,Ny)
  Sdirt::Array{T,2} = zeros(Nx,Ny)
  Sdift::Array{T,2} = zeros(Nx,Ny)
  SWveg::Array{T,2} = zeros(Nx,Ny)
  SWsrf::Array{T,2} = zeros(Nx,Ny)
  SWsci::Array{T,2} = zeros(Nx,Ny)
  LWt::Array{T,2} = zeros(Nx,Ny)
  SWtopo_out::Array{T,2} = zeros(Nx,Ny)

  # Thermal - temporary arrays

  ksnow::Array{T,3} = zeros(Nsmax, Nx, Ny)
  csoil::Array{T,3} = zeros(Nsoil, Nx, Ny)
  ksoil::Array{T,3} = zeros(Nsoil, Nx, Ny)
  gs1::Array{T,2} = zeros(Nx, Ny)
  Ds1::Array{T,2} = zeros(Nx, Ny)
  Ts1::Array{T,2} = zeros(Nx, Ny)
  ks1::Array{T,2} = zeros(Nx, Ny)
  Tveg0::Array{T,2} = zeros(Nx, Ny)

  # Sfexch - temporary arrays

  KH::Array{T,2} = zeros(Nx, Ny)
  KHa::Array{T,2} = zeros(Nx, Ny)
  KHg::Array{T,2} = zeros(Nx, Ny)
  KHv::Array{T,2} = zeros(Nx, Ny)
  KWg::Array{T,2} = zeros(Nx, Ny)
  KWv::Array{T,2} = zeros(Nx, Ny)
  Usc::Array{T,2} = zeros(Nx, Ny)

  # Ebalsrf - temporary arrays
  dTs::Array{T, 2} = zeros(Nx,Ny)
  Esrf::Array{T,2} = zeros(Nx,Ny)
  Eveg::Array{T,2} = zeros(Nx,Ny)
  G::Array{T,2} = zeros(Nx,Ny)
  H::Array{T,2} = zeros(Nx,Ny)
  Hsrf::Array{T,2} = zeros(Nx,Ny)
  LE::Array{T,2} = zeros(Nx,Ny)
  LEsrf::Array{T,2} = zeros(Nx,Ny)
  LWsci::Array{T,2} = zeros(Nx,Ny)
  LWveg::Array{T,2} = zeros(Nx,Ny)
  Melt::Array{T,2} = zeros(Nx,Ny)
  Rnet::Array{T,2} = zeros(Nx,Ny)
  Rsrf::Array{T,2} = zeros(Nx,Ny)

  # Snow - temporary arrays

  Gsoil::Array{T,2} = zeros(Nx, Ny)
  Roff::Array{T,2} = zeros(Nx, Ny)
  meltflux_out::Array{T,2} = zeros(Nx, Ny)
  Sbsrf::Array{T,2} = zeros(Nx, Ny)
  Roff_bare::Array{T, 2} = zeros(Nx, Ny)
  Roff_snow::Array{T, 2} = zeros(Nx, Ny)
  fsnow_thres::Array{T, 2} = zeros(Nx, Ny)
  unload::Array{T, 2} = zeros(Nx, Ny)

  a::Vector{T} = zeros(Nsmax)
  bsnow::Vector{T} = zeros(Nsmax)
  c::Vector{T} = zeros(Nsmax)
  csnow::Vector{T} = zeros(Nsmax)
  dTssnow::Vector{T} = zeros(Nsmax)
  D::Vector{T} = zeros(Nsmax)
  E::Vector{T} = zeros(Nsmax)
  Gs::Vector{T} = zeros(Nsmax)
  rhs::Vector{T} = zeros(Nsmax)
  R::Vector{T} = zeros(Nsmax)
  S::Vector{T} = zeros(Nsmax)
  U::Vector{T} = zeros(Nsmax)
  W::Vector{T} = zeros(Nsmax)

end