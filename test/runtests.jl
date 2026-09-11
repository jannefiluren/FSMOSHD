using Test
using FlexibleSnowModelOSHD

@testset "Mass Balance" begin
    include("test_mass_balance.jl")
end

@testset "MET Immutability" begin
    include("test_met_immutability.jl")
end

@testset "Soil Energy Balance" begin
    include("test_soil_energy_balance.jl")
end

@testset "Architectures" begin
    include("test_architectures.jl")
end

@testset "SnowSlide" begin
    include("test_snowslide.jl")
end

@testset "SnowTran3D" begin
    include("test_snowtran3d.jl")
end

@testset "Transport step" begin
    include("test_transport_step.jl")
end

@testset "Regression Tests" begin
    include("test_regression.jl")
end

# GPU test — runs a full CPU-vs-GPU comparison only if CUDA is functional in the (stacked)
# default environment; otherwise it logs a skip and adds no tests. See test_gpu.jl.
include("test_gpu.jl")
