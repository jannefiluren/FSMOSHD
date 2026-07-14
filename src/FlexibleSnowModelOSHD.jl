module FlexibleSnowModelOSHD

using Parameters
using Dates

# KernelAbstractions is imported qualified because it exports its own CPU/GPU
# names, which would clash with the architecture types defined here
import KernelAbstractions
using KernelAbstractions: @kernel, @index, @Const, get_backend

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
export AbstractArchitecture, CPU, GPU, on_architecture
export canopy!, radiation!, thermal!, sfexch!, ebalsrf!, ebalfor!, snow!, soil!, snowcoverfraction!
export qsat, tridiag!, ludcmp!
export drive!, step!, setup
export @unpack_constants

end
