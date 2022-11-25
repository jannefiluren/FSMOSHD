using Test
using FSMOSHD
using DelimitedFiles

@testset "complete_season" begin
  
  for station in ["SLF.5WJ", "MCH.BLS2", "MCH.OTE2"]

    julia_run = run_fsm(station)

    fortran_run = readdlm("../fortran/output_64/output_" * replace(station, "." => "_") * "_run_from_julia.txt")

    @test maximum(abs.(julia_run .- fortran_run)) < 10e-8

  end

end
