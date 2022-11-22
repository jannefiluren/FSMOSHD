state_file = ""

include("../src/parameters.jl")
include("../src/setup.jl")
include("../src/initialize.jl")
include("../src/qsat.jl")
include("../src/tridiag.jl")
include("../src/drive.jl")
include("../src/radiation.jl")
include("../src/thermal.jl")
include("../src/sfexch.jl")
include("../src/ebalsrf.jl")
include("../src/snow.jl")
include("../src/soil.jl")

function run_fsm_julia(station)

  drive_file = "../fortran/input/input_" * replace(station, "." => "_") * ".txt"
  output_file = "../fortran/output_julia/output_" * replace(station, "." => "_") * "_test.txt"

  terrain = readline("../fortran/input/terrain_" * replace(station,"." => "_") * ".txt")
  terrain = parse.(Float64,split(terrain,","))
  
  fsky_terr[:,:] .= terrain[1]
  slopemu[:,:] .= terrain[2]
  xi[:,:] .= terrain[3]
  Ld[:,:] .= terrain[4]
  lat[:,:] .= terrain[5]
  lon[:,:] .= terrain[6]
  dem[:,:] .= terrain[7]

  global Qa

  Qa = similar(Ta)

  fout = open(output_file, "w")

  for (index, data) in enumerate(readlines(drive_file))

    println("Time step: ", index)

    ### Run drive

    global year, month, day, hour

    year, month, day, hour = drive(data)

    ### Run radiation

    radiation()

    ### Run thermal

    thermal()

    ### Run sfexch and ebalsrf

    for i in 1:Nitr
      sfexch()
      ebalsrf()
    end

    # Run snow

    snow()

    # Run soil

    soil()

    println(fout, "$(year) $(month) $(day) $(hour) $(sum(Ds[:,1,1])) $(fsnow[1,1]) $(sum(Sice[:,1,1]+Sliq[:,1,1])) $(Tsrf[1,1]) $(Nsnow[1,1])")

  end

  close(fout)

end