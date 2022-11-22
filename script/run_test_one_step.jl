# ./run_test_one_step.bat

station = "MCH.OTE2"

drive_file = "../fortran/input/input_" * replace(station, "." => "_") * "_final_step.txt"
state_file = "../fortran/temp/states_end_64.txt"
check_final_vals = true
output_file = ""

include("run_fsm_oshd.jl")