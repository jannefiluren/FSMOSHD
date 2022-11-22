# FSMOSHD

Short instructions on how to find errors between fortran and julia:

1) Use compare_julia.m to identify the exact time step where the difference arises.

2) Create two input files that runs the model until the difference arises (input_station_first_half.txt and input_station_final_step.txt)

3) Adjust the two namelist files: nlst_first_half_64.nam and nlst_input_fake_one_timestep_64.nam

4) Do not forget to update the terrain characteristics in the namelist files

5) Run the fortran code using run_test_one_step.bat

6) Run the julia script run_test_one_step.jl

7) Run the last timestep with both codes and use prints and stops to debug.

8) Run fortran using ./FSM2_TXT_64.exe nlst_input_fake_one_timestep_64.nam
