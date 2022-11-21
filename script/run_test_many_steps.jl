# ./run_test_many_steps.bat

drive_file = "../fortran/input/input_fake_long_SLF_5WJ.txt"
state_file = "../fortran/temp/states_end_64.txt"
check_final_vals = true
output_file = ""

include("run_fsm_oshd.jl")