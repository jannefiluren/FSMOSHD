# ./run_all_whole_season.bat

station = ARGS[1]

drive_file = "../fortran/input/input_" * replace(station, "." => "_") * ".txt"
output_file = "../fortran/output_julia/output_" * replace(station, "." => "_") * "_test.txt"
terrain_file = "../fortran/input/terrain_" * replace(station, "." => "_") * ".txt"
state_file = ""
check_final_vals = false

include("run_fsm_oshd.jl")