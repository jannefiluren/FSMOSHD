using Test
using FSMOSHD
using DelimitedFiles

include("../script/run_fsm_fortran.jl")
include("../script/compile_fsm.jl")

@testset "complete_season" begin
  
  path = pwd()

  cd("../fortran/fsm_txt_64/")

  compile()

  cd(path)

  for station in ["SLF.5WJ", "MCH.BLS2", "MCH.MAG2", "MCH.OTE2", "MCH.SCD2", "MCH.LUN2", "MCH.JUN2"]

    run_fsm_fortran(station)

    julia_run = run_fsm(station)

    fortran_run = readdlm("../fortran/output_64/output_" * replace(station, "." => "_") * "_run_from_julia.txt")

    @test maximum(abs.(julia_run .- fortran_run)) < 10e-8

  end

end
