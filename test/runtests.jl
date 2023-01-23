using Test
using FSMOSHD
using DelimitedFiles

#absolute path to project folder
projdir = dirname(dirname(@__FILE__))

include(joinpath(projdir, "script", "run_fsm_fortran.jl"))
include(joinpath(projdir, "script", "compile_fsm.jl"))

@testset "complete_season" begin

  cd(joinpath(projdir, "fortran", "fsm_txt_64"))

  compile()

  mkpath(joinpath(projdir, "fortran", "output_64/"))

  for station in ["SLF.5WJ", "MCH.BLS2", "MCH.MAG2", "MCH.OTE2", "MCH.SCD2", "MCH.LUN2", "MCH.JUN2"]

    run_fsm_fortran(station)

    julia_run = run_fsm_point(station)

    fortran_run = readdlm(joinpath(projdir, "fortran", "output_64", "output_") * replace(station, "." => "_") * "_run_from_julia.txt")

    @test maximum(abs.(julia_run .- fortran_run)) < 10e-8

  end

end
