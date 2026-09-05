using Adapt
using FlexibleSnowModelOSHD
using Test

# Stand-in for a device array type (CuArray and friends). Adapt performs the same
# reconstruction at every GPU kernel launch, so this exercises it without a GPU.
struct DeviceArray{T, N} <: AbstractArray{T, N}
    data::Array{T, N}
end
Base.size(a::DeviceArray) = size(a.data)
Base.getindex(a::DeviceArray, i...) = a.data[i...]

struct DeviceAdaptor end
Adapt.adapt_storage(::DeviceAdaptor, a::Array) = DeviceArray(a)

@testset "on_architecture round trip" begin

    fsm = FSM{Float32, Int32}(; Nx = 3, Ny = 2)
    met = MET{Float32, Int32}(Nx = 3, Ny = 2)

    moved = on_architecture(FlexibleSnowModelOSHD.CPU(), fsm)
    @test typeof(moved) === typeof(fsm)
    @test moved.grid.Dzsnow == fsm.grid.Dzsnow
    @test moved.landuse.alb0 == fsm.landuse.alb0
    @test moved.state.Ds == fsm.state.Ds
    @test moved.diag.ksnow == fsm.diag.ksnow
    @test moved.physics == fsm.physics

    movedmet = on_architecture(FlexibleSnowModelOSHD.CPU(), met)
    @test typeof(movedmet) === typeof(met)
    @test movedmet.Nx == met.Nx
    @test movedmet.Sf_history_f64 == met.Sf_history_f64

end

@testset "adapt rebuilds sub-structs on a device array type" begin

    fsm = FSM{Float32, Int32}(; Nx = 3, Ny = 2)
    met = MET{Float32, Int32}(Nx = 3, Ny = 2)
    to = DeviceAdaptor()

    grid = Adapt.adapt(to, fsm.grid)
    @test grid.Dzsnow isa DeviceArray{Float32, 1}
    @test grid.Nsmax === fsm.grid.Nsmax

    # Tf is the declared type of no field on these three, so the reconstruction has
    # to recover it from a representative array
    landuse = Adapt.adapt(to, fsm.landuse)
    @test landuse isa FlexibleSnowModelOSHD.Landuse
    @test landuse.alb0 isa DeviceArray{Float32, 2}
    @test landuse.prec_multi isa DeviceArray{Float64, 2}

    state = Adapt.adapt(to, fsm.state)
    @test state isa FlexibleSnowModelOSHD.State
    @test state.albs isa DeviceArray{Float32, 2}
    @test state.Nsnow isa DeviceArray{Int32, 2}
    @test state.Ds isa DeviceArray{Float32, 3}

    diag = Adapt.adapt(to, fsm.diag)
    @test diag isa FlexibleSnowModelOSHD.Diagnostics
    @test diag.es isa DeviceArray{Float32, 2}
    @test diag.ksnow isa DeviceArray{Float32, 3}

    m = Adapt.adapt(to, met)
    @test m isa MET{Float32, Int32}
    @test m.Sdir isa DeviceArray{Float32, 2}
    @test m.Sf24h_f64 isa DeviceArray{Float64, 2}
    @test m.Sf_history_f64 isa DeviceArray{Float64, 3}

end
