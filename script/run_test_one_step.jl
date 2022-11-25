# ./run_test_one_step.bat

station = "SLF.5WJ"

drive_file = "../fortran/input/input_fake_one_timestep_" * replace(station, "." => "_") * ".txt"
terrain_file = "../fortran/input/terrain_" * replace(station, "." => "_") * ".txt"
state_file = "../fortran/temp/states_end_64.txt"
check_final_vals = true
output_file = ""

include("run_fsm_oshd.jl")