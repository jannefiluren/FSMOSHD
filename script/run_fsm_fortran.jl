function run_fsm_fortran(station)

write_namelist(station)

projdir = dirname(dirname(@__FILE__))

cd(joinpath(projdir, "fortran"))
if Sys.iswindows()
  run(`FSM2_txt_64.exe nlst_from_julia.nam`)
else
  run(`./FSM2_TXT_64 nlst_from_julia.nam`)
end
cd(joinpath(projdir, "script"))

end


function write_namelist(station; NALBEDO=2, NCANMOD=0, NCONDCT=1, NDENSTY=3, NEXCHNG=1, NHYDROL=2, NSNFRAC=3, NRADSBG=0, NZOFFST=0, NOSHDTN=1)

  projdir = dirname(dirname(@__FILE__))
  terrain = readline(joinpath(projdir, "fortran", "input", "terrain_") * replace(station,"." => "_") * ".txt")
  terrain = parse.(Float64,split(terrain,","))
  
  fsky_terr = terrain[1]
  slopemu = terrain[2]
  xi = terrain[3]
  Ld = terrain[4]
  lat = terrain[5]
  lon = terrain[6]
  dem = terrain[7]
  
  in_file = "input/input_" * replace(station,"." => "_") * ".txt"
  out_file = "output_64/output_" * replace(station,"." => "_") * "_run_from_julia.txt"
  
  nlst = """
  &nam_grid
    NNx = 1,
    NNy = 1,
    NNsmax = 3,
    NNsoil = 4,
  /
  &nam_layers
    DDzsnow = 0.1, 0.2, 0.4,
    DDzsoil = 0.1, 0.2, 0.4, 0.8,
  /
  &nam_driving
    zzT = 10,
    zzU = 10,
    met_file = '$(in_file)',
    out_file = '$(out_file)',
  /
  &nam_modconf
    NALBEDO = $(NALBEDO),
    NCANMOD = $(NCANMOD),
    NCONDCT = $(NCONDCT),
    NDENSTY = $(NDENSTY),
    NEXCHNG = $(NEXCHNG),
    NHYDROL = $(NHYDROL),
    NSNFRAC = $(NSNFRAC),
    NRADSBG = $(NRADSBG),
    NZOFFST = $(NZOFFST),
    NOSHDTN = $(NOSHDTN),
    LHN_ON  = .FALSE.,
    LFOR_HN = .TRUE.,
  /
  &nam_modpert
    LZ0PERT = .FALSE.,
  /
  &nam_modtile
    CTILE = 'open',
    rtthresh = 0.1,
  /
  &nam_results
     CLIST_DIAG_RESULTS = 'rotc', 'hsnt', 'swet', 'slqt', 'swtb', 'swtd', 'lwtr', 'romc', 'sbsc',
     CLIST_STATE_RESULTS = 'tsfe', 'scfe',
  /
  &nam_location
    fsky_terr = $(fsky_terr),
    slopemu = $(slopemu),
    xi = $(xi),
    Ld = $(Ld),
    lat = $(lat),
    lon = $(lon),
    dem = $(dem),
    pmultf = 1,
  /
  """

  open(joinpath(projdir, "fortran", "nlst_from_julia.nam"), "w") do io
    write(io,nlst)
  end;
  
  end

