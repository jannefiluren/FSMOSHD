# Initialization of variables
Nx = 1
Ny = 1
Nsmax = 3
Nsoil = 4

Dzsnow = [0.1, 0.2, 0.4]
Dzsoil = [0.1, 0.2, 0.4, 0.8]

# Driving data
dt = 1
zT = 10
zU = 10

# Model configuration
ALBEDO = 2
CANMOD = 0
CONDCT = 1
DENSTY = 3
EXCHNG = 1
HYDROL = 2
SNFRAC = 3
RADSBG = 0
ZOFFST = 0
OSHDTN = 1
HN_ON = false
FOR_HN = true


tthresh = 0.1    ### hack


# Model perturbations
Z0PERT = false
WCPERT = false
FSPERT = false
ALPERT = false
SLPERT = false

# Modelled tile
TILE = "open"
rtthresh = 0.1

# Outputs
Nave = 1

# Defaults for numerical solution parameters
Nitr = 4


# Defaults for canopy parameters
avg0 = 0.1
avgs = 0.4
cden = 0.004
cvai = 6.6
cveg = 20
Gcn1 = 0.5
Gcn2 = 0
gsnf = 0
kdif = 0.5
kveg = 1
rchd = 0.67
rchz = 0.2          
tcnc = 240
tcnm = 48

# Defaults for snow parameters
a_eta = 0.1
asmx = 0.8
asmn = 0.5
b_eta = 0.023
bstb = 5
bthr = 2
c_eta = 250
eta0 = 3.7e7
eta1 = 7.62237e6
hfsn = 0.1
kfix = 0.24
rho0 = 300
rhob = 6
rhoc = 26
rhof = 109
rcld = 300
rgr0 = 5e-5
rmlt = 500
Salb = 10
snda = 2.8e-6
Talb = -2
tcld = 1000
tmlt = 100
trho = 200
Wirr = 0.03
z0sn = 0.002
Sfmin = 10

# some defaults different for forest tile - commented-out values based on FS-EBP runs, revisit during tuning
if (TILE == "forest")
  hfsn = 0.3    
  z0sn = 0.005
end

# Defaults for ground surface parameters
bstb = 5
gsat = 0.01

# Defaults for additional parameters required for forest snow process parametrization
adfs = 3
adfl = 2
fsar = 0.1
psf  = 1
psr  = 0.1
wcan = 2.5
zsub = 2
zgf = 5
zgr = 5
khcf = 3

if (DENSTY == 0)
  rhof = rho0
end

# Surface properties
if (TILE == "glacier")
  alb0 = 0.3*ones(Nx,Ny)
  z0sf = 0.04*ones(Nx,Ny)
else
  alb0 = 0.2*ones(Nx,Ny)
  z0sf = 0.2*ones(Nx,Ny)
end

fcly = 0.3*ones(Nx,Ny)
fsnd = 0.6*ones(Nx,Ny)

if (TILE == "forest")
  z0sf = 0.2*ones(Nx,Ny)
end

# Canopy parameters
canh  = zeros(Nx,Ny)   # todo may need to change
fsky  = zeros(Nx,Ny)   # todo may need to change
fveg  = zeros(Nx,Ny)   # todo may need to change
fves  = zeros(Nx,Ny)   # todo may need to change
hcan  = zeros(Nx,Ny)   # todo may need to change
lai  = zeros(Nx,Ny)   # todo may need to change
pmultf = zeros(Nx,Ny)   # todo may need to change
scap  = zeros(Nx,Ny)   # todo may need to change
trcn  = zeros(Nx,Ny)   # todo may need to change
VAI   = zeros(Nx,Ny)   # todo may need to change
vfhp   = zeros(Nx,Ny)   # todo may need to change

# Allocate terrain properties
slopemu = zeros(Nx,Ny)
xi = zeros(Nx,Ny)
Ld = zeros(Nx,Ny)
lat = zeros(Nx,Ny)
lon = zeros(Nx,Ny)
dem = zeros(Nx,Ny)
tilefrac = zeros(Nx,Ny)

# Derived soil parameters
b = zeros(Nx,Ny)
hcap_soil = zeros(Nx,Ny)
hcon_soil = zeros(Nx,Ny)
sathh = zeros(Nx,Ny)
Vsat = zeros(Nx,Ny)
Vcrit = zeros(Nx,Ny)

for j = 1:Ny
  for i = 1:Nx
    if (fcly[i,j] + fsnd[i,j] > 1)
      fcly[i,j] = 1 - fsnd[i,j]
    end
    b[i,j] = 3.1 + 15.7*fcly[i,j] - 0.3*fsnd[i,j]
    hcap_soil[i,j] = (2.128*fcly[i,j] + 2.385*fsnd[i,j])*1e6 / (fcly[i,j] + fsnd[i,j])
    sathh[i,j] = 10^(0.17 - 0.63*fcly[i,j] - 1.58*fsnd[i,j])
    Vsat[i,j] = 0.505 - 0.037*fcly[i,j] - 0.142*fsnd[i,j]
    Vcrit[i,j] = Vsat[i,j]*(sathh[i,j]/3.364)^(1/b[i,j])
    hcon_min = (hcon_clay^fcly[i,j]) * (hcon_sand^(1 - fcly[i,j]))
    hcon_soil[i,j] = (hcon_air^Vsat[i,j]) * (hcon_min^(1 - Vsat[i,j]))
  end
end

# Convert time scales from hours to seconds
dt = 3600*dt
tcnc = 3600*tcnc
tcnm = 3600*tcnm
tcld = 3600*tcld
tmlt = 3600*tmlt
trho = 3600*trho

# Allocate state variables
albs = 0.85.*ones(Nx,Ny)
Ds = zeros(Nsmax,Nx,Ny)
Nsnow = zeros(Int64,Nx,Ny)
Qcan = zeros(Nx,Ny)
rgrn = zeros(Nsmax,Nx,Ny)
Sice = zeros(Nsmax,Nx,Ny)
Sliq = zeros(Nsmax,Nx,Ny)
Sveg = zeros(Nx,Ny)
Tcan = 273.15.*ones(Nx,Ny)
theta = zeros(Nsoil,Nx,Ny)
Tsnow = 273.15.*ones(Nsmax,Nx,Ny)
Tsoil = 273.15.*ones(Nsoil,Nx,Ny)
Tsrf = zeros(Nx,Ny)
fsnow = zeros(Nx,Ny)
Tveg = 273.15.*ones(Nx,Ny)
snowdepthmin = zeros(Nx,Ny)
snowdepthmax = zeros(Nx,Ny)
snowdepthhist = zeros(14,Nx,Ny)
swemin = zeros(Nx,Ny)
swemax = zeros(Nx,Ny)
swehist = zeros(14,Nx,Ny)
fsky_terr = zeros(Nx,Ny)

# Initial soil profiles

fsat = zeros(Nsoil)
Tprof = zeros(Nsoil)
fsat[:] .= 0.5
Tprof[:] .= 285
for k = 1:Nsoil
  theta[k,:,:] .= fsat[k]*Vsat[:,:]
  Tsoil[k,:,:] .= Tprof[k]
end
Tsrf[:,:] .= Tsoil[1,:,:]

# Read state variables

if isfile(state_file)
    tmp = readlines(state_file)

    albs[1,1] = parse.(Float64,split(tmp[1]))[1]
    Ds[:,1,1] .= parse.(Float64,split(tmp[2]))
    Nsnow[1,1] = parse.(Float64,split(tmp[3]))[1]
    Qcan[1,1] = parse.(Float64,split(tmp[4]))[1]
    Sice[:,1,1] .= parse.(Float64,split(tmp[5]))
    Sliq[:,1,1] .= parse.(Float64,split(tmp[6]))
    Sveg[1,1] = parse.(Float64,split(tmp[7]))[1]
    Tcan[1,1] = parse.(Float64,split(tmp[8]))[1]
    theta[:,1,1] .= parse.(Float64,split(tmp[9]))
    Tsnow[:,1,1] .= parse.(Float64,split(tmp[10]))
    Tsoil[:,1,1] .= parse.(Float64,split(tmp[11]))
    Tsrf[1,1] = parse.(Float64,split(tmp[12]))[1]
    fsnow[1,1] = parse.(Float64,split(tmp[13]))[1]
    Tveg[1,1] = parse.(Float64,split(tmp[14]))[1]
    snowdepthmin[1,1] = parse.(Float64,split(tmp[15]))[1]
    snowdepthmax[1,1] = parse.(Float64,split(tmp[16]))[1]
    snowdepthhist[:,1,1] .= parse.(Float64,split(tmp[17]))
    swemin[1,1] = parse.(Float64,split(tmp[18]))[1]
    swemax[1,1] = parse.(Float64,split(tmp[19]))[1]
    swehist[:,1,1] .= parse.(Float64,split(tmp[20]))

end

# println(albs)
# println(Ds)
# println(Nsnow)
# println(Qcan)
# println(Sice)
# println(Sliq)
# println(Sveg)
# println(Tcan)
# println(theta)
# println(Tsnow)
# println(Tsoil)
# println(Tsrf)
# println(fsnow)
# println(Tveg)
# println(snowdepthmin)
# println(snowdepthmax)
# println(snowdepthhist)
# println(swemin)
# println(swemax)
# println(swehist)


# ! Cap glacier temperatures to 0°C
# if (TILE == 'glacier') then
#   do j = 1, Ny
#   do i = 1, Nx
#     Tsrf(i,j) = min(Tsrf(i,j),Tm)
#     do k = 1, Nsoil
#       Tsoil(k,i,j) = min(Tsoil(k,i,j),Tm)
#     end do
#   end do
#   end do
# endif

# ! model tile fractions 
# if (TILE == 'open') then 
#   tilefrac = dem/dem   ! temporary fix to get ones within our entire domain, assuming we always want to run an open tile. may have to be revisited
# else 
#   read(1139) tilefrac
# endif 

# if (TILE == "open")
#   tilefrac .= dem./dem
# end
  
  
  # if (SNFRAC == 0 .or. SNFRAC == 2) then
  #   read(1110) snowdepthmax
  # endif
  
  # if (SNFRAC == 0 ) then
  #   ! states specific to open runs
  #   read(1109) snowdepthmin
  #   read(1111) snowdepthhist
  #   read(1113) swemin
  #   read(1114) swemax
  #   read(1115) swehist
  #   read(1124) slopemu
  #   read(1125) xi
  #   read(1126) Ld
  # endif
  
  if (TILE != "forest") 
    # canopy properties (no canopy)
    VAI[:,:]  .= 0
    hcan[:,:] .= 0
    fsky[:,:] .= 1
    trcn[:,:] .= exp(-kdif.*VAI[:,:])
    fveg[:,:] .= 1 .- exp(-kveg.*VAI[:,:])
    fves[:,:] .= 1 .- exp(-kveg.*VAI[:,:])
  else
    # # lus fields specific to forest runs
    # read(1130) Qcan
    # read(1131) Sveg
    # read(1132) Tcan
    # read(1133) Tveg
    # read(1134) fveg
    # read(1135) hcan
    # read(1136) lai
    # read(1137) vfhp
    # read(1138) fves
    # read(1140) pmultf
  
    # ! derived canopy properties 
    # VAI(:,:) = lai(:,:) 
    # trcn(:,:) = 1-0.9*fveg(:,:)  
    # do j = 1, Ny
    #   do i = 1, Nx
    #     fsky(i,j) = vfhp(i,j)/trcn(i,j)
    #     if ( fsky(i,j) > 1 ) trcn(i,j) = vfhp(i,j)
    #     if ( fsky(i,j) > 1 ) fsky(i,j) = 1
    #   end do
    # end do 
  end
  
  # derived canopy parameters
  canh[:,:] .= 12500*VAI[:,:]
  scap[:,:] .= cvai*VAI[:,:]

  # hacks...
  fsky_terr[:,:] .= 0.9483
  slopemu[:,:] .= 0.0122
  xi[:,:] .= 0
  Ld[:,:] .= 1
  lat[:,:] .= 46.8296
  lon[:,:] .= 9.8092
  dem[:,:] .= 2540
  pmultf[:,:] .= 1
  tilefrac[:,:] .= 1

  year = 0
  month = 0
  day = 0
  hour = 0
  Sdir = zeros(Nx,Ny)
  Sdif = zeros(Nx,Ny)
  LW = zeros(Nx,Ny)
  Sf = zeros(Nx,Ny)
  Rf = zeros(Nx,Ny)
  Ta = zeros(Nx,Ny)
  RH = zeros(Nx,Ny)
  Ua = zeros(Nx,Ny)
  Ps = zeros(Nx,Ny)
  Sf24h = zeros(Nx,Ny)
  Tv = zeros(Nx,Ny)
  
  # # prints
  
  # println("Nsmax: ", Nsmax)
  # println("Nsoil: ", Nsoil)
  # println("Nx: ", Nx)
  # println("Ny: ", Ny)
  
  # println("Dzsnow: ", Dzsnow)
  # println("Dzsoil: ", Dzsoil)
  
  # println("dt: ", dt)
  # println("zT: ", zT)
  # println("zU: ", zU)
  
  # println("ALBEDO: ", ALBEDO)
  # println("CANMOD: ", CANMOD)
  # println("CONDCT: ", CONDCT)
  # println("DENSTY: ", DENSTY)
  # println("EXCHNG: ", EXCHNG)
  # println("HYDROL: ", HYDROL)
  # println("SNFRAC: ", SNFRAC)
  # println("RADSBG: ", RADSBG)
  # println("ZOFFST: ", ZOFFST)
  # println("OSHDTN: ", OSHDTN)
  # println("HN_ON: ", HN_ON)
  # println("FOR_HN: ", FOR_HN)
  
  # println("fcly: ", fcly)
  # println("b: ", b)
  # println("hcap_soil: ", hcap_soil)
  # println("sathh: ", sathh)
  # println("Vsat: ", Vsat)
  # println("Vcrit: ", Vcrit)
  # println("hcon_soil: ", hcon_soil)
  
  # println("dt: ", dt)
  # println("tcnc: ", tcnc)
  # println("tcnm: ", tcnm)
  # println("tcld: ", tcld)
  # println("tmlt: ", tmlt)
  # println("trho: ", trho)
  
  # println("albs: ", albs)
  # println("Ds: ", Ds)
  # println("fsnow: ", fsnow)
  # println("Nsnow: ", Nsnow)
  # println("Sice: ", Sice)
  # println("Sliq: ", Sliq)
  # println("Tsrf: ", Tsrf)
  # println("Tsnow: ", Tsnow)
  # println("Tsoil: ", Tsoil)
  # println("fsky_terr: ", fsky_terr)
  # println("lat: ", lat)
  # println("lon: ", lon)
  # println("dem: ", dem)
  
  # println("VAI: ", VAI)
  # println("hcan: ", hcan)
  # println("fsky: ", fsky)
  # println("trcn: ", trcn)
  # println("fveg: ", fveg)
  # println("fves: ", fves)
  
  # println("fsat: ", fsat)
  # println("Tprof: ", Tprof)
  # println("theta: ", theta)
  # println("Tsoil: ", Tsoil)
  # println("Tsrf: ", Tsrf)
  
  # exit()