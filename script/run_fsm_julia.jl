using Pkg
Pkg.activate("../.")

using FSMOSHD
using Parameters

include("../src/parameters.jl")
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

Qa = similar(Ta)


function run_fsm_julia(station)

  drive_file = "../fortran/input/input_" * replace(station, "." => "_") * ".txt"
  output_file = "../fortran/output_julia/output_" * replace(station, "." => "_") * "_test.txt"
  terrain_file = "../fortran/input/terrain_" * replace(station,"." => "_") * ".txt"

  fsm = FSM{Float64}()
  setup!(fsm, terrain_file)
  
  fout = open(output_file, "w")

  for (index, data) in enumerate(readlines(drive_file))

    println("Time step: ", index)

    ### Run drive

    global year, month, day, hour

    year, month, day, hour = drive(fsm, data)

    ### Run radiation

    radiation(fsm)

    ### Run thermal

    thermal(fsm)

    ### Run sfexch and ebalsrf

    for i in 1:fsm.Nitr
      sfexch(fsm)
      ebalsrf(fsm)
    end

    # Run snow

    snow(fsm)

    # Run soil

    soil(fsm)

    println(fout, "$(year) $(month) $(day) $(hour) $(sum(fsm.Ds[:,1,1])) $(fsm.fsnow[1,1]) $(sum(fsm.Sice[:,1,1]+fsm.Sliq[:,1,1])) $(fsm.Tsrf[1,1]) $(fsm.Nsnow[1,1])")

  end

  close(fout)

end