using Test
using FSMOSHD


@testset "setup" begin
  
  fsm = FSM{Float64}()

  setup!(fsm, "../fortran/input/terrain_SLF_5WJ.txt")


end
