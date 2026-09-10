"""
Architecture abstraction for CPU/GPU portability (Oceananigans-style).

An architecture decides where the model arrays live. `CPU()` is the default
and keeps everything as plain Julia `Array`s. `GPU(backend)` wraps a
KernelAbstractions GPU backend, e.g. with CUDA.jl loaded:

    using CUDA
    arch = GPU(CUDABackend())
    fsm = setup(arch, Grid(Float32; Nx = Nx, Ny = Ny), surface, settings)
    met = on_architecture(arch, MET{Float32}(Nx = Nx, Ny = Ny))

Physics routines pick their compute backend from the arrays themselves (via
`KernelAbstractions.get_backend`).
"""
abstract type AbstractArchitecture end

"""
    CPU()

Architecture for CPU runs: all model arrays are plain Julia `Array`s.
KernelAbstractions kernels are partitioned across Julia threads.
"""
struct CPU <: AbstractArchitecture end

"""
    GPU(backend)

Architecture for GPU runs on a KernelAbstractions GPU backend
(e.g. `CUDABackend()` from CUDA.jl, `ROCBackend()` from AMDGPU.jl).
"""
struct GPU{B} <: AbstractArchitecture
    backend::B
end

backend(::CPU) = KernelAbstractions.CPU()
backend(arch::GPU) = arch.backend

"""
    on_architecture(arch, x)

Move `x` to the architecture `arch`: arrays are converted to the
architecture's array type (a no-op for `Array`s on `CPU()`), all
structures are rebuilt with all their array fields converted, and everything
else is passed through unchanged. Note that unconverted fields (and arrays
already on the right architecture) are aliased, not copied.
"""
on_architecture(::AbstractArchitecture, x) = x

on_architecture(::CPU, a::Array) = a
on_architecture(::CPU, a::AbstractArray) = Array(a)

function on_architecture(arch::GPU, a::AbstractArray)
    out = KernelAbstractions.allocate(backend(arch), eltype(a), size(a))
    copyto!(out, a)
    return out
end

function on_architecture(arch::AbstractArchitecture, s::AbstractParameterization)
    isbitstype(typeof(s)) && return s
    T = typeof(s).name.wrapper
    return T(map(f -> on_architecture(arch, getfield(s, f)), fieldnames(typeof(s)))...)
end

on_architecture(::AbstractArchitecture, p::Parameters) = p

function on_architecture(arch::AbstractArchitecture, x::Union{Grid, Surface, State, Diagnostics, MET})
    T = typeof(x).name.wrapper
    return T(map(f -> on_architecture(arch, getfield(x, f)), fieldnames(typeof(x)))...)
end

on_architecture(arch::AbstractArchitecture, nt::NamedTuple) = map(x -> on_architecture(arch, x), nt)

function on_architecture(arch::AbstractArchitecture, fsm::FSM)
    values = map(name -> on_architecture(arch, getfield(fsm, name)), fieldnames(typeof(fsm)))
    return FSM(values...)
end
