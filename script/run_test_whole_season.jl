# ./run_all_whole_season.bat

drive_file = "../fortran/input/input_SLF_5WJ.txt"
state_file = ""
check_final_vals = false
output_file = "../fortran/output_julia/output_SLF_5WJ.txt"

include("run_fsm_oshd.jl")