using FSMOSHD

station = "SLF.5WJ"

projdir = dirname(dirname(@__FILE__))

terrain_file = joinpath(projdir, "fortran", "input", "terrain_SLF_5WJ.txt")
state_file = joinpath(projdir, "fortran", "temp", "states_end_64.txt")

fsm = FSM{Float64}()
setup!(fsm, terrain_file) #, state_file=state_file)

year = 2022
month = 3
day = 10
hour = 12
Sdir = 50*ones(fsm.Nx, fsm.Ny)
Sdif = 50*ones(fsm.Nx, fsm.Ny)
LW = 300*ones(fsm.Nx, fsm.Ny)
Sf = 1*ones(fsm.Nx, fsm.Ny)
Rf = 1*ones(fsm.Nx, fsm.Ny)
Ta = 5*ones(fsm.Nx, fsm.Ny)
RH = 50*ones(fsm.Nx, fsm.Ny)
Ua = 3*ones(fsm.Nx, fsm.Ny)
Ps = 73238*ones(fsm.Nx, fsm.Ny)
Sf24h = 2*ones(fsm.Nx, fsm.Ny)
Tc = zeros(fsm.Nx, fsm.Ny)
es = zeros(fsm.Nx, fsm.Ny)
Qa = zeros(fsm.Nx, fsm.Ny)
Tv = zeros(fsm.Nx, fsm.Ny)

drive!(fsm, Tc, es, Qa, Ua, Sf, Rf, Ta, RH, Ps)

@code_warntype ebalsrf(fsm, LW, Ps, Qa, Ta)

using BenchmarkTools
@btime ebalsrf(fsm, LW, Ps, Qa, Ta)

using Profile, PProf
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=1 ebalsrf(fsm, LW, Ps, Qa, Ta)
results = Profile.Allocs.fetch();
PProf.Allocs.pprof(results; from_c=false)