# ./FSM2_TXT_64.exe nlst_one_step_no_snow_64.nam

station = "MCH.OTE2"

drive_file = "../fortran/input/input_fake_one_step_no_snow.txt"
output_file = ""
state_file = ""
check_final_vals = true

include("run_fsm_oshd.jl")