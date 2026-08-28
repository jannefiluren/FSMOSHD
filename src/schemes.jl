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
    build_scheme(Tf, requested, Nx, Ny, params)

Instantiate a parameterization. `requested` is either a type - constructed at precision `Tf` on an
`Nx` by `Ny` grid, consuming any `params` entry named like one of its fields - or a ready-made
instance, returned unchanged.

Routing matters: a scheme's parameters live on the scheme, so an entry such as "adm" would
otherwise be set on `FSM`, where nothing reads it any more.
"""
function build_scheme(Tf, requested, Nx, Ny, params)

    requested isa Type || return requested

    kwargs = Dict{Symbol, Any}()
    for name in fieldnames(requested)
        key = string(name)
        haskey(params, key) && (kwargs[name] = pop!(params, key))
    end

    return requested{Tf}(Nx, Ny; kwargs...)

end
