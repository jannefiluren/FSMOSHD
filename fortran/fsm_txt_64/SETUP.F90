!-----------------------------------------------------------------------
! Set parameters, initialize prognostic variables and write metadata
!-----------------------------------------------------------------------
subroutine SETUP

use, intrinsic :: iso_fortran_env, only: dp=>real64 

use MODCONF, only: ALBEDO,CANMOD,CONDCT,DENSTY,EXCHNG,HYDROL,&
SNFRAC,RADSBG,ZOFFST,OSHDTN,HN_ON,FOR_HN

use MODPERT, only: Z0PERT,WCPERT,FSPERT,ALPERT,SLPERT

use MODTILE, only: TILE, tthresh 

use MODOUTPUT, only: LIST_DIAG_RESULTS, LIST_STATE_RESULTS

use CONSTANTS

use DIAGNOSTICS, only : Nave

use DRIVING, only: &
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  RH,                &! Relative humidity (%)
  Rf,                &! Rainfall rate (kg/m2/s)
  Sf,                &! Snowfall rate (kg/m2/s)
  Sf24h,             &! Snowfall 24hr (kg/m2)
  Sdir,              &! Incoming direct beam radiation on flat,unobstructed surface (W/m2)
  Sdif,              &! Incoming diffuse radiation on flat,unobstructed (W/m2)
  Ta,                &! Air temperature (K)
  Tv,                &! Time-varying canopy transmissivity for dSWR (-)
  Ua,                &! Wind speed (m/s)
  zT,                &! Temperature measurement height (m)
  zU,                &! Wind measurement height (m)		
  z0P,               &! z0 perturbations
  wcP,               &! liquid water capacity perturbations
  fsP,               &! fresh snow density perturbations
  alP,               &! albedo perturbations
  slP                 ! settling perturbations
  
use GRID

use PARAMETERS

use PARAMMAPS

use SOILPARAMS 

use STATE_VARIABLES

use LANDUSE

use IOUNITS, only : &
  uout,              &! Output file unit number
  umet,              &! Driving file unit number
  udmp,              &! Dump file unit number
  ustr                ! Start file unit number

use FILES, only : dump_file

use TEST

implicit none
 
integer :: & 
  i,j,               &! Point counters
  k,                 &! Level counter
  iresults_count
  
integer :: &
  NNsmax,NNsoil,NNx,NNy
  
integer :: &
  NALBEDO,NCANMOD,NCONDCT,NDENSTY,NEXCHNG,NHYDROL,NSNFRAC,NRADSBG,NZOFFST,NOSHDTN
  
real*8 :: &
  zzT,zzU

real*8 :: &
  rtthresh

real*8 :: &
  hcon_min            ! Thermal conductivity of soil minerals (W/m/K)
  
real*8, allocatable :: &
  fsat(:),           &! Initial moisture content of soil layers as fractions of saturation
  Tprof(:)            ! Initial soil layer temperatures (K)

character(len=200) :: &
  nlst_file,        & ! Namelist file name
  met_file,         & ! Drive file name
  out_file,         & ! Output file name
  start_file          ! Start file name


character(len=20) :: CTILE
  
character(len=4), dimension(36) :: CLIST_DIAG_RESULTS

character(len=4), dimension(12) :: CLIST_STATE_RESULTS

real*8, dimension(:), allocatable :: &
  DDzsnow,DDzsoil

logical :: LHN_ON  ! activate the HN model

logical :: LZ0PERT,LWCPERT,LFSPERT,LALPERT,LSLPERT,LFOR_HN

logical :: lexist

!-1- !!!!!!!!!!!!!!!!!!!!  READ THE NAMELIST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
namelist  /nam_grid/    NNx,NNy,NNsmax,NNsoil
namelist  /nam_layers/  DDzsnow,DDzsoil
namelist  /nam_driving/ zzT,zzU,met_file,out_file,start_file,dump_file,dump_data
namelist  /nam_modconf/ NALBEDO,NCANMOD,NCONDCT,NDENSTY,NEXCHNG,NHYDROL,NSNFRAC,NRADSBG,NZOFFST,NOSHDTN,LHN_ON,LFOR_HN
namelist  /nam_modtile/ CTILE, rtthresh
namelist  /nam_modpert/ LZ0PERT,LWCPERT,LFSPERT,LALPERT,LSLPERT
namelist  /nam_results/ CLIST_DIAG_RESULTS, CLIST_STATE_RESULTS
namelist  /nam_location/ fsky_terr,slopemu,xi,Ld,lat,lon,dem,pmultf

! get namelist path from first command argument.
call getarg(1, nlst_file)
INQUIRE (FILE=nlst_file, EXIST=lexist)
if (lexist) then
  open(5000, file = nlst_file)
else
  print*, '  namelist file: ', trim(nlst_file), ' does not exist'
  call exit(1)
endif

! Initialization of variables
NNsmax = 3
NNsoil = 4 ! Attention, 5 soil layers in JIM ... 
NNx = 1
NNy = 1
read(5000, nam_grid)
Nsmax = NNsmax
Nsoil = NNsoil
Nx = NNx
Ny = NNy

allocate(DDzsnow(Nsmax))
allocate(DDzsoil(Nsoil))
allocate(Dzsnow(Nsmax))
allocate(Dzsoil(Nsoil))
read(5000, nam_layers)
Dzsnow = DDzsnow
Dzsoil = DDzsoil

! Driving data
dt = 1_dp
zzT = 10_dp
zzU = 10_dp
start_file = 'none'
dump_file = 'none'
dump_data = 0
read(5000,nam_driving)
zT = zzT
zU = zzU

! Model configuration
! -1 for mandatory NLST arguments.
NALBEDO = -1
NCANMOD = -1
NCONDCT = -1
NDENSTY = -1
NEXCHNG = -1
NHYDROL = -1
NSNFRAC = -1
NRADSBG = -1
NZOFFST = -1
NOSHDTN = -1
LHN_ON = .FALSE.
LFOR_HN = .FALSE.
read(5000, nam_modconf)
ALBEDO = NALBEDO
CANMOD = NCANMOD
CONDCT = NCONDCT
DENSTY = NDENSTY
EXCHNG = NEXCHNG
HYDROL = NHYDROL
SNFRAC = NSNFRAC
RADSBG = NRADSBG
ZOFFST = NZOFFST
OSHDTN = NOSHDTN
HN_ON = LHN_ON
FOR_HN = LFOR_HN

if (ALBEDO==-1 .or. CANMOD==-1 .or. CONDCT==-1 .or. DENSTY==-1 .or. EXCHNG==-1 &
.or. HYDROL==-1 .or. SNFRAC==-1 .or. RADSBG ==-1 .or. ZOFFST ==-1 .or. OSHDTN ==-1) then
  print*, 'model configuration error:\n please specify all the fields of MODCONF in the namelist (&nam_modconf)'
  call exit(1)
endif

! Model perturbations
LZ0PERT = .FALSE.
LWCPERT = .FALSE.
LFSPERT = .FALSE.
LALPERT = .FALSE.
LSLPERT = .FALSE.
read(5000, nam_modpert)
Z0PERT = LZ0PERT
WCPERT = LWCPERT
FSPERT = LFSPERT
ALPERT = LALPERT
SLPERT = LSLPERT

! Modelled tile 
TILE = 'open'  ! set open as default if no tile has been defined 
tthresh = 0.1
read(5000, nam_modtile)
TILE = CTILE
tthresh = rtthresh

! List of output variables
CLIST_DIAG_RESULTS(:)  = '    ' ! BC 13-char length
CLIST_STATE_RESULTS(:) = '    ' ! BC 13-char length
read(5000, nam_results)
iresults_count = count(CLIST_DIAG_RESULTS /= '             ') ! BC 13-char length
LIST_DIAG_RESULTS(1:iresults_count) = CLIST_DIAG_RESULTS(1:iresults_count)
iresults_count = count(CLIST_STATE_RESULTS /= '             ') ! BC 13-char length
LIST_STATE_RESULTS(1:iresults_count) = CLIST_STATE_RESULTS(1:iresults_count)

! Outputs
Nave = 1 !24 Set to one, 24h averages created in OSHD Matlab wrapper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-2- !!!!!!!!!!!!!!!!!!!!    OPEN THE FILES   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open Files for OSHD states input and output and results output
!call OPEN_FILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(umet, file = met_file)
open(uout, file = out_file)


!-3- !!!!!!!!!!!!!!!!!!!!  ALLOCATE VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(LW(Nx,Ny))
allocate(Ps(Nx,Ny))
allocate(Qa(Nx,Ny))
allocate(Rf(Nx,Ny))
allocate(Sf(Nx,Ny))
allocate(Sdir(Nx,Ny))
allocate(Sdif(Nx,Ny))
allocate(Ta(Nx,Ny))
allocate(Ua(Nx,Ny))
allocate(RH(Nx,Ny))
allocate(Sf24h(Nx,Ny))
allocate(Tv(Nx,Ny))

if (Z0PERT) allocate(z0P(Nx,Ny))
if (WCPERT) allocate(wcP(Nx,Ny))
if (FSPERT) allocate(fsP(Nx,Ny))
if (ALPERT) allocate(alP(Nx,Ny))
if (SLPERT) allocate(slP(Nx,Ny))

! use Tv dummy in case of open simulations
if (CANMOD == 0) then
  Tv(:,:) = 1
endif

! Defaults for numerical solution parameters
Nitr = 4

! Defaults for canopy parameters
avg0 = 0.1_dp
avgs = 0.4_dp
cden = 0.004_dp
cvai = 6.6_dp
cveg = 20_dp
Gcn1 = 0.5_dp
Gcn2 = 0_dp
gsnf = 0_dp
kdif = 0.5_dp
kveg = 1_dp
rchd = 0.67_dp
rchz = 0.2_dp
tcnc = 240_dp
tcnm = 48_dp

! Defaults for snow parameters
a_eta = 0.1_dp
asmx = 0.8_dp       ! unused if OSHDTN = 1
asmn = 0.5_dp
b_eta = 0.023_dp
bstb = 5_dp
bthr = 2_dp
c_eta = 250_dp
eta0 = 3.7e7_dp
eta1 = 7.62237e6_dp
hfsn = 0.1_dp
kfix = 0.24_dp
rho0 = 300_dp
rhob = 6_dp
rhoc = 26_dp
rhof = 109_dp
rcld = 300_dp
rgr0 = 5e-5_dp
rmlt = 500_dp
snda = 2.8e-6_dp
Talb = -2_dp
tcld = 1000_dp
tmlt = 100_dp
trho = 200_dp
Wirr = 0.03_dp
z0sn = 0.002_dp
Sfmin = 10_dp

! some defaults different for forest tile - commented-out values based on FS-EBP runs, revisit during tuning
if (TILE == 'forest') then
  ! asmx = 0.88_dp
  hfsn = 0.3_dp    
  z0sn = 0.005_dp   
endif 

! Defaults for ground surface parameters
bstb = 5_dp
gsat = 0.01_dp

! Defaults for additional parameters required for forest snow process parametrization
adfs = 3_dp
adfl = 2_dp
fsar = 0.1_dp
psf  = 1_dp
psr  = 0.1_dp
wcan = 2.5_dp
zsub = 2_dp
zgf = 5_dp
zgr = 5_dp
khcf = 3_dp

if (DENSTY == 0) then
  rhof = rho0
endif

! Surface data from defaults or namelists  
! Surface properties 
allocate(alb0(Nx,Ny))
allocate(fcly(Nx,Ny))
allocate(fsnd(Nx,Ny))
allocate(z0sf(Nx,Ny))
if (TILE == 'glacier') then
  alb0(:,:) = 0.3_dp
  z0sf(:,:) = 0.04_dp
else
  alb0(:,:) = 0.2_dp
  z0sf(:,:) = 0.2_dp
endif
fcly(:,:) = 0.3_dp
fsnd(:,:) = 0.6_dp
if (TILE == 'forest') z0sf(:,:) = 0.2_dp

! Canopy parameters
allocate(canh(Nx,Ny))
allocate(fsky(Nx,Ny))
allocate(fveg(Nx,Ny)) 
allocate(fves(Nx,Ny))
allocate(hcan(Nx,Ny))
allocate(lai(Nx,Ny))
allocate(pmultf(Nx,Ny))
allocate(scap(Nx,Ny))
allocate(trcn(Nx,Ny))
allocate(VAI(Nx,Ny))
allocate(vfhp(Nx,Ny))
canh(:,:)  = undef
fsky(:,:)  = undef
fveg(:,:)  = undef
fves(:,:)  = undef
hcan(:,:)  = undef
pmultf(:,:) = undef
scap(:,:)  = undef
trcn(:,:)  = undef
VAI(:,:)   = undef

!Terrain properties
allocate(slopemu(Nx,Ny))
allocate(xi(Nx,Ny))
allocate(Ld(Nx,Ny))
allocate(lat(Nx,Ny))
allocate(lon(Nx,Ny))
allocate(dem(Nx,Ny))
allocate(tilefrac(Nx,Ny))
slopemu(:,:) = undef
xi(:,:) = undef
Ld(:,:) = undef
lat(:,:) = undef
lon(:,:) = undef
dem(:,:) = undef
tilefrac(:,:) = undef

! Derived soil parameters
allocate(b(Nx,Ny))
allocate(hcap_soil(Nx,Ny))
allocate(hcon_soil(Nx,Ny))
allocate(sathh(Nx,Ny))
allocate(Vsat(Nx,Ny))
allocate(Vcrit(Nx,Ny))
do j = 1, Ny
  do i = 1, Nx
    if (fcly(i,j) + fsnd(i,j) > 1) then
      fcly(i,j) = 1_dp - fsnd(i,j)
    endif
    b(i,j) = 3.1_dp + 15.7_dp*fcly(i,j) - 0.3_dp*fsnd(i,j)
    hcap_soil(i,j) = (2.128_dp*fcly(i,j) + 2.385_dp*fsnd(i,j))*1e6_dp / (fcly(i,j) + fsnd(i,j))
    sathh(i,j) = 10_dp**(0.17_dp - 0.63_dp*fcly(i,j) - 1.58_dp*fsnd(i,j))
    Vsat(i,j) = 0.505_dp - 0.037_dp*fcly(i,j) - 0.142_dp*fsnd(i,j)
    Vcrit(i,j) = Vsat(i,j)*(sathh(i,j)/3.364_dp)**(1_dp/b(i,j))
    hcon_min = (hcon_clay**fcly(i,j)) * (hcon_sand**(1_dp - fcly(i,j)))
    hcon_soil(i,j) = (hcon_air**Vsat(i,j)) * (hcon_min**(1_dp - Vsat(i,j)))
  end do
end do

! Convert time scales from hours to seconds
dt = 3600_dp*dt
tcnc = 3600_dp*tcnc
tcnm = 3600_dp*tcnm
tcld = 3600_dp*tcld
tmlt = 3600_dp*tmlt
trho = 3600_dp*trho

! Allocate state variables
allocate(albs(Nx,Ny))
allocate(Ds(Nsmax,Nx,Ny))
allocate(Nsnow(Nx,Ny))
allocate(Qcan(Nx,Ny))
allocate(rgrn(Nsmax,Nx,Ny))
allocate(Sice(Nsmax,Nx,Ny))
allocate(Sliq(Nsmax,Nx,Ny))
allocate(Sveg(Nx,Ny))
allocate(Tcan(Nx,Ny))
allocate(theta(Nsoil,Nx,Ny))
allocate(Tsnow(Nsmax,Nx,Ny))
allocate(Tsoil(Nsoil,Nx,Ny))
allocate(Tsrf(Nx,Ny))
allocate(fsnow(Nx,Ny))
allocate(Tveg(Nx,Ny))
allocate(snowdepthmin(Nx,Ny))
allocate(snowdepthmax(Nx,Ny))
allocate(snowdepthhist(14,Nx,Ny))
allocate(swemin(Nx,Ny))
allocate(swemax(Nx,Ny))
allocate(swehist(14,Nx,Ny))
allocate(fsky_terr(Nx,Ny))

! Default initialization of state variables 
albs(:,:)    = 0.8500_dp
Ds(:,:,:)    = 0_dp
fsnow(:,:)   = 0_dp
Nsnow(:,:)   = 0_dp
Qcan(:,:)    = 0_dp
rgrn(:,:,:)  = undef !*GM watch out: rgrn currently not tracked
Sice(:,:,:)  = 0_dp
Sliq(:,:,:)  = 0_dp
Sveg(:,:)    = 0_dp
Tcan(:,:)    = 273.1500_dp
Tsnow(:,:,:) = 273.1500_dp
Tsoil(:,:,:) = 273.1500_dp
Tveg(:,:)    = 273.1500_dp
snowdepthmin(:,:) = 0_dp
snowdepthmax(:,:) = 0_dp
snowdepthhist(:,:,:) = 0_dp
swemin(:,:) = 0_dp
swemax(:,:) = 0_dp
swehist(:,:,:) = 0_dp

! Initial soil profiles from namelist
allocate(fsat(Nsoil))
allocate(Tprof(Nsoil))
fsat(:)  = 0.5_dp
Tprof(:) = 285_dp
do k = 1, Nsoil
  theta(k,:,:) = fsat(k)*Vsat(:,:)
  Tsoil(k,:,:) = Tprof(k)
end do
Tsrf(:,:) = Tsoil(1,:,:)

! Initialize state variables from a named start file
if (start_file /= 'none') then
  open(ustr, file = start_file)
  read(ustr,*) albs(:,:)
  read(ustr,*) Ds(:,:,:)
  read(ustr,*) Nsnow(:,:)
  read(ustr,*) Qcan(:,:)
  read(ustr,*) Sice(:,:,:)
  read(ustr,*) Sliq(:,:,:)
  read(ustr,*) Sveg(:,:)
  read(ustr,*) Tcan(:,:)
  read(ustr,*) theta(:,:,:)
  read(ustr,*) Tsnow(:,:,:)
  read(ustr,*) Tsoil(:,:,:)
  read(ustr,*) Tsrf(:,:)
  read(ustr,*) fsnow(:,:)
  read(ustr,*) Tveg(:,:)
  read(ustr,*) snowdepthmin(:,:)
  read(ustr,*) snowdepthmax(:,:)
  read(ustr,*) snowdepthhist(:,:,:)
  read(ustr,*) swemin(:,:)
  read(ustr,*) swemax(:,:)
  read(ustr,*) swehist(:,:,:)

  close(ustr)
end if

! print *,  albs
! print *,  Ds
! print *,  Nsnow
! print *,  Qcan
! print *,  Sice
! print *,  Sliq
! print *,  Sveg
! print *,  Tcan
! print *,  theta
! print *,  Tsnow
! print *,  Tsoil
! print *,  Tsrf
! print *,  fsnow
! print *,  Tveg
! print *,  snowdepthmin
! print *,  snowdepthmax
! print *,  snowdepthhist
! print *,  swemin
! print *,  swemax
! print *,  swehist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !-4- !!!!!!!!!!!!!!!!!!!! READ DRIVING/STATES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! states relevant to both open and forest simulation
! read(1101) albs
! read(1102) Ds
! read(1103) fsnow
! read(1104) Nsnow
! read(1106) Sice
! read(1107) Sliq
! read(1116) Tsrf
! read(1119) Tsnow
! read(1120) Tsoil
! read(1123) fsky_terr
! read(1127) lat
! read(1128) lon
! read(1129) dem

read(5000, nam_location)

! Cap glacier temperatures to 0°C
if (TILE == 'glacier') then
  do j = 1, Ny
  do i = 1, Nx
    Tsrf(i,j) = min(Tsrf(i,j),Tm)
    do k = 1, Nsoil
      Tsoil(k,i,j) = min(Tsoil(k,i,j),Tm)
    end do
  end do
  end do
endif

! model tile fractions 
! if (TILE == 'open') then 
tilefrac = dem/dem   ! temporary fix to get ones within our entire domain, assuming we always want to run an open tile. may have to be revisited
! else 
!   read(1139) tilefrac
! endif 

if (SNFRAC == 0 .or. SNFRAC == 2) then
  ! read(1110) snowdepthmax
endif

if (SNFRAC == 0 ) then
  ! ! states specific to open runs
  ! read(1109) snowdepthmin
  ! read(1111) snowdepthhist
  ! read(1113) swemin
  ! read(1114) swemax
  ! read(1115) swehist
  ! read(1124) slopemu
  ! read(1125) xi
  ! read(1126) Ld
endif

if (TILE /= 'forest') then
  ! canopy properties (no canopy)
  VAI(:,:)  = 0_dp
  hcan(:,:) = 0_dp
  fsky(:,:) = 1_dp
  trcn(:,:) = exp(-kdif*VAI(:,:))
  fveg(:,:) = 1_dp - exp(-kveg*VAI(:,:))
  fves(:,:) = 1_dp - exp(-kveg*VAI(:,:))
else ! TILE == 'forest'
  ! ! lus fields specific to forest runs
  ! read(1130) Qcan
  ! read(1131) Sveg
  ! read(1132) Tcan
  ! read(1133) Tveg
  ! read(1134) fveg
  ! read(1135) hcan
  ! read(1136) lai
  ! read(1137) vfhp
  ! read(1138) fves
  ! read(1140) pmultf

  ! derived canopy properties 
  VAI(:,:) = lai(:,:) 
  trcn(:,:) = 1_dp-0.9_dp*fveg(:,:)  
  do j = 1, Ny
    do i = 1, Nx
      fsky(i,j) = vfhp(i,j)/trcn(i,j)
      if ( fsky(i,j) > 1 ) trcn(i,j) = vfhp(i,j)
      if ( fsky(i,j) > 1 ) fsky(i,j) = 1_dp
    end do
  end do 
endif

! derived canopy parameters
canh(:,:) = 12500_dp*VAI(:,:)
scap(:,:) = cvai*VAI(:,:)

close(5000)

if (dump_data == 1) then
  open(1, file="temp/test_setup.txt")
  
  write(1,*)  albs
  write(1,*)  Ds
  write(1,*)  Nsnow
  write(1,*)  Qcan
  write(1,*)  Sice
  write(1,*)  Sliq
  write(1,*)  Sveg
  write(1,*)  Tcan
  write(1,*)  theta
  write(1,*)  Tsnow
  write(1,*)  Tsoil
  write(1,*)  Tsrf
  write(1,*)  fsnow
  write(1,*)  Tveg
  write(1,*)  snowdepthmin
  write(1,*)  snowdepthmax
  write(1,*)  snowdepthhist
  write(1,*)  swemin
  write(1,*)  swemax
  write(1,*)  swehist
  write(1,*)  fcly
  write(1,*)  b
  write(1,*)  hcap_soil
  write(1,*)  sathh
  write(1,*)  Vsat
  write(1,*)  Vcrit
  close(1) 
end if

end subroutine SETUP
