#
# Shared machinery for physics parameterizations.
#
# Every scheme derives from AbstractParameterization (FlexibleSnowModelOSHD.jl), which is what
# lets on_architecture (architectures.jl) move any of them with one generic method.
#

"""
    grid_array(Tf, x, Nx, Ny)

Materialize a parameterization parameter as an `Nx` by `Ny` array of element type `Tf`.
A scalar is broadcast over the whole grid; an array is converted element-wise.
"""
grid_array(::Type{Tf}, x::Number, Nx, Ny) where {Tf} = fill(Tf(x), Nx, Ny)
grid_array(::Type{Tf}, x::AbstractArray, Nx, Ny) where {Tf} = convert(Array{Tf, 2}, x)

"""
    check_grid(scheme, Nx, Ny)

Assert that any grid-shaped parameter held by `scheme` matches the `Nx` by `Ny` model grid.
The fallback accepts anything, so a parameterization built only from scalars needs no method.

Without this a mismatch is not caught at setup; it surfaces later as a `BoundsError` from
inside a kernel, which says nothing about the actual cause.
"""
check_grid(scheme, Nx, Ny) = nothing

function check_grid(scheme::AbstractParameterization, Nx, Ny)
    for name in fieldnames(typeof(scheme))
        value = getfield(scheme, name)
        value isa AbstractArray || continue
        size(value) == (Nx, Ny) || throw(
            DimensionMismatch(
                "$(nameof(typeof(scheme))) field `$name` is $(size(value)) but the model grid is ($Nx, $Ny)"
            )
        )
    end
    return nothing
end

"""
    instantiate(scheme, grid)

Return a physics parameterization ready for the model: a scheme **type** is default-constructed at
the grid's precision (`Scheme{eltype(grid)}(grid)`), while a ready-made **instance** is returned
unchanged. Non-default schemes are constructed by the caller, e.g.
`PrognosticAlbedo{Float32}(grid; adm = 200, adc = my_array)`.
"""
instantiate(scheme::Type, grid) = scheme{eltype(grid)}(grid)
instantiate(scheme, grid) = scheme

"""
    reconstruct(x; kwargs...)

Copy the immutable struct `x` with the named fields replaced. `Parameters` is rebuilt
rather than mutated so that it stays isbits and can cross into a kernel by value.
"""
function reconstruct(x::T; kwargs...) where {T}
    names = fieldnames(T)
    fields = NamedTuple{names}(map(f -> getfield(x, f), names))
    return T(; merge(fields, NamedTuple(kwargs))...)
end

# Route a scalar into the immutable Parameters via a functional update.
@inline function set_param!(fsm, sym::Symbol, value)
    v = convert(fieldtype(typeof(fsm.params), sym), value)
    fsm.params = reconstruct(fsm.params; NamedTuple{(sym,)}((v,))...)
    return fsm
end

"""
    apply_params!(fsm, params)

Apply parameter overrides: a `Parameters` scalar is reconstructed onto `fsm.params`;
a `Surface` per-cell array is filled (scalar) or copied (array) in place.
"""
function apply_params!(fsm, params)
    for (key, value) in params
        sym = Symbol(key)
        if hasfield(typeof(fsm.params), sym)
            set_param!(fsm, sym, value)
        elseif hasfield(typeof(fsm.surface), sym)
            arr = getfield(fsm.surface, sym)
            value isa AbstractArray ? (arr .= eltype(arr).(value)) : fill!(arr, eltype(arr)(value))
        else
            throw(ArgumentError("unknown parameter override \"$key\""))
        end
    end
    return fsm
end
