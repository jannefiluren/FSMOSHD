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

@testset "Regression Tests" begin
    include("test_regression.jl")
end
