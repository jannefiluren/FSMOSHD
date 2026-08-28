module FlexibleSnowModelOSHD

using Parameters
using Dates

abstract type AbstractConductivity{Tf <: Real} end

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
include("types.jl")
include("architectures.jl")
include("kernel_utils.jl")
include("setup.jl")
include("qsat.jl")
include("tridiag.jl")
include("ludcmp.jl")
include("fresh_snow_density.jl")
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
export AbstractConductivity, FixedConductivity, DensityConductivity, snow_conductivity!
export AbstractArchitecture, CPU, GPU, on_architecture
export canopy!, radiation!, thermal!, sfexch!, ebalsrf!, ebalfor!, snow!, soil!, snowcoverfraction!
export qsat, tridiag!, ludcmp!
export drive!, step!, setup
export @unpack_constants

end
