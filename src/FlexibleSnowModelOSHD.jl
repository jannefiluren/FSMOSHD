module FlexibleSnowModelOSHD

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
import Libdl

# Core: physical constants, precision-parametric structs, architecture, parameter overrides, setup
include("parameters.jl")
include("schemes.jl")
include("types.jl")
include("architectures.jl")
include("kernel_utils.jl")
include("setup.jl")

# Numerical utilities (grid-agnostic solvers and math)
include("numerics/qsat.jl")
include("numerics/tridiag.jl")
include("numerics/ludcmp.jl")

# Parameterizations: one AbstractParameterization scheme family per file. Included before the
# processes, whose kernels dispatch on / call these schemes.
include("parameterizations/albedo.jl")
include("parameterizations/conductivity.jl")
include("parameterizations/canopy.jl")
include("parameterizations/substrate.jl")
include("parameterizations/reference_height.jl")
include("parameterizations/stability.jl")
include("parameterizations/surface_layer.jl")
include("parameterizations/fresh_snow_density.jl")
include("parameterizations/snow_compaction.jl")
include("parameterizations/snow_hydrology.jl")
include("parameterizations/snowcoverfraction.jl")
include("parameterizations/snow_relayering.jl")

# Processes: each defines a launcher foo!(fsm, ...) + its kernel, called from step!
include("processes/drive.jl")
include("processes/canopy.jl")
include("processes/radiation.jl")
include("processes/thermal.jl")
include("processes/surface_exchange_coefficients.jl")
include("processes/surface_energy_balance.jl")
include("processes/snow.jl")
include("processes/soil.jl")
include("step.jl")

# Snow transport (SnowSlide + SnowTran3D). Standalone, CPU-only operators - NOT part of
# step!. Both a Fortran (ccall) and a pure-Julia implementation are kept so the two can be
# cross-validated as the Fortran evolves. The Fortran libraries are built by deps/build.jl
# (Pkg.build); the ccall wrappers resolve them lazily, so the module loads without them.
include("transport/transport_types.jl")
include("transport/transport_setup.jl")
include("transport/snowslide.jl")
include("transport/snowslide_julia.jl")
include("transport/snowtran3d.jl")
include("transport/snowtran3d_julia.jl")
include("transport/relayer.jl")
include("transport/transport.jl")

export FSM, MET, Grid
export AbstractParameterization, grid_array, check_grid, instantiate
export AbstractConductivity, FixedConductivity, DensityConductivity, snow_conductivity!
export AbstractAlbedo, DiagnosticAlbedo, DecayAlbedo, PrognosticAlbedo, snow_albedo!
export AbstractCanopy, NoCanopy, OneLayerCanopy, surface_energy_balance!, energy_balance!, canopy_snow!
export solar_radiation!, thermal_radiation!
export AbstractSubstrate, SoilSubstrate, IceSubstrate, soil_properties!
export AbstractReferenceHeight, AboveGround, AboveCanopy
export AbstractSurfaceLayer, OpenSurfaceLayer, ForestSurfaceLayer
export AbstractStabilityCorrection, NoStabilityCorrection, LouisStabilityCorrection
export AbstractFreshSnowDensity, FixedFreshSnowDensity, ClimateFreshSnowDensity, ElevationFreshSnowDensity, fresh_snow_density
export AbstractCompaction, AgeCompaction, OverburdenCompaction, CrocusCompaction, compact_snow!
export AbstractHydrology, FreeDrainingHydrology, BucketHydrology, DensityBucketHydrology, snow_hydrology!
export AbstractLayering, OriginalLayering, DensityLayering, relayer_snow!
export AbstractSnowFraction, SeasonalSnowFraction, HelbigSnowFraction, HelbigMaxSnowFraction, PointSnowFraction, TanhSnowFraction, snow_covered_fraction!
export AbstractArchitecture, CPU, GPU, on_architecture
export canopy!, radiation!, thermal!, surface_exchange_coefficients!, snow!, soil!, snowcoverfraction!
export qsat, tridiag!, ludcmp!
export drive!, step!, setup
export SnowTransport, setup_transport, transport!, relayer!
export snowslide!, snowslide_julia!, snowtran3d!, snowtran3d_julia!
export @unpack_constants

end
