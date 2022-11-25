# ./FSM2_TXT_64.exe nlst_one_step_no_snow_64.nam

station = "SLF.5WJ"

drive_file = "../fortran/input/input_fake_one_step_no_snow.txt"
terrain_file = "../fortran/input/terrain_" * replace(station, "." => "_") * ".txt"
output_file = ""
state_file = ""
check_final_vals = true

include("run_fsm_oshd.jl")