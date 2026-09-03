module FlexibleSnowModelOSHD

using Parameters
using Dates
using Adapt: Adapt, @adapt_structure

# Every physics parameterization derives from this, which is what lets
# on_architecture (architectures.jl) move any of them in one generic method.
abstract type AbstractParameterization{Tf <: Real} end
abstract type AbstractConductivity{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractAlbedo{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractCanopy{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractSubstrate{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractReferenceHeight{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractSurfaceLayer{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractStabilityCorrection{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractFreshSnowDensity{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractCompaction{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractHydrology{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractLayering{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractSnowFraction{Tf} <: AbstractParameterization{Tf} end

# KernelAbstractions is imported qualified because it exports its own CPU/GPU
# names, which would clash with the architecture types defined here.
# NOTE: @Const is deliberately NOT used on kernel arguments - on the CPU
# backend it wraps the kernel body in an aliasscope, which miscompiles large
# kernel bodies on Julia >= 1.11 when combined with inbounds (silently wrong
# results; see JuliaGPU/KernelAbstractions.jl#652 for the aliasscope issue).
import KernelAbstractions
using KernelAbstractions: @kernel, @index, get_backend
using StaticArrays: MVector, MMatrix

include("parameters.jl")
include("schemes.jl")
include("types.jl")
include("architectures.jl")
include("kernel_utils.jl")
include("setup.jl")
include("qsat.jl")
include("tridiag.jl")
include("ludcmp.jl")
include("fresh_snow_density.jl")
include("snow_compaction.jl")
include("snow_hydrology.jl")
include("snow_relayering.jl")
include("snow_layering.jl")
include("drive.jl")
include("canopy.jl")
include("radiation.jl")
include("thermal.jl")
include("sfexch.jl")
include("ebalsrf.jl")
include("ebalfor.jl")
include("snow.jl")
include("soil.jl")
include("step.jl")
include("snowcoverfraction.jl")

export FSM, MET
export AbstractParameterization, grid_array, check_grid, build_scheme
export AbstractConductivity, FixedConductivity, DensityConductivity, snow_conductivity!
export AbstractAlbedo, DiagnosticAlbedo, DecayAlbedo, PrognosticAlbedo, snow_albedo!
export AbstractCanopy, NoCanopy, OneLayerCanopy, surface_energy_balance!
export AbstractSubstrate, SoilSubstrate, IceSubstrate
export AbstractReferenceHeight, AboveGround, AboveCanopy
export AbstractSurfaceLayer, OpenSurfaceLayer, ForestSurfaceLayer
export AbstractStabilityCorrection, NoStabilityCorrection, LouisStabilityCorrection
export AbstractFreshSnowDensity, FixedFreshSnowDensity, ClimateFreshSnowDensity, ElevationFreshSnowDensity, fresh_snow_density
export AbstractCompaction, AgeCompaction, OverburdenCompaction, CrocusCompaction, compact_snow!
export AbstractHydrology, FreeDrainingHydrology, BucketHydrology, DensityBucketHydrology, snow_hydrology!
export AbstractLayering, OriginalLayering, DensityLayering, relayer_snow!
export AbstractSnowFraction, SeasonalSnowFraction, HelbigSnowFraction, HelbigMaxSnowFraction, PointSnowFraction, TanhSnowFraction, snow_covered_fraction!
export AbstractArchitecture, CPU, GPU, on_architecture
export canopy!, radiation!, thermal!, sfexch!, ebalsrf!, ebalfor!, snow!, soil!, snowcoverfraction!
export qsat, tridiag!, ludcmp!
export drive!, step!, setup
export @unpack_constants

end
