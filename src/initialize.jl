@unpack Nx, Ny, Nsmax, Nsoil, albs = fsm

# Radiation

alb = similar(albs)   ### hack
asrf_out = similar(albs)   ### hack
Sdirt = similar(albs)   ### hack
Sdift = similar(albs)   ### hack
SWveg = similar(albs)   ### hack
SWsrf = similar(albs)   ### hack
SWsci = similar(albs)   ### hack
LWt = similar(albs)   ### hack
SWtopo_out = similar(albs)   ### hack

# Thermal

ksnow = zeros(Nsmax, Nx, Ny)
csoil = zeros(Nsoil, Nx, Ny)
ksoil = zeros(Nsoil, Nx, Ny)
gs1 = zeros(Nx, Ny)
Ds1 = zeros(Nx, Ny)
Ts1 = zeros(Nx, Ny)
ks1 = zeros(Nx, Ny)
Tveg0 = zeros(Nx, Ny)

# Sfexch

KH = zeros(Nx, Ny)
KHa = zeros(Nx, Ny)
KHg = zeros(Nx, Ny)
KHv = zeros(Nx, Ny)
KWg = zeros(Nx, Ny)
KWv = zeros(Nx, Ny)
Usc = zeros(Nx, Ny)

# EBALSRF

Esrf = zeros(Nx,Ny)
Eveg = zeros(Nx,Ny)
G = zeros(Nx,Ny)
H = zeros(Nx,Ny)
Hsrf = zeros(Nx,Ny)
LE = zeros(Nx,Ny)
LEsrf = zeros(Nx,Ny)
LWsci = zeros(Nx,Ny)
LWveg = zeros(Nx,Ny)
Melt = zeros(Nx,Ny)
Rnet = zeros(Nx,Ny)
Rsrf = zeros(Nx,Ny)

# Snow

Gsoil = zeros(Nx, Ny)
Roff = zeros(Nx, Ny)
meltflux_out = zeros(Nx, Ny)
Sbsrf = zeros(Nx, Ny)
Roff_bare = zeros(Nx, Ny)
Roff_snow = zeros(Nx, Ny)
fsnow_thres = zeros(Nx, Ny)
unload = zeros(Nx, Ny)

# Drive variables
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
