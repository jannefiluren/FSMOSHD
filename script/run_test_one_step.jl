# ./run_test_one_step.bat

drive_file = "../fortran/data/input_fake_one_timestep_5wj.txt"
state_file = "../fortran/data/states_end_64.txt"
check_final_vals = true
output_file = ""

include("run_fsm_oshd.jl")