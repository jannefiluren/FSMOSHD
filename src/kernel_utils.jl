"""
    column_sum(A, i, j)

Sum of the layer dimension of a `(Nlayer, Nx, Ny)` array at pixel (i, j).
Equivalent to `sum(@view A[:, i, j])` (same order of operations, so results
are bit-identical), written as a plain loop that is safe and fast inside
KernelAbstractions kernels on both CPU and GPU.
"""
@inline function column_sum(A::AbstractArray{Tf, 3}, i::Integer, j::Integer) where {Tf}
    s = zero(Tf)
    @inbounds for k in axes(A, 1)
        s += A[k, i, j]
    end
    return s
end

"""
    first_argmin(v, n)

Index of the first minimum of `v[1:n]`. Equivalent to `argmin(v[1:n])` for
data without NaNs, written as a plain loop that is safe inside kernels
(range indexing of an `MVector` would allocate).
"""
@inline function first_argmin(v, n::Integer)
    im = 1
    @inbounds for k in 2:n
        if v[k] < v[im]
            im = k
        end
    end
    return im
end

"""
    first_argmax(v, n)

Index of the first maximum of `v[1:n]`; see [`first_argmin`](@ref).
"""
@inline function first_argmax(v, n::Integer)
    im = 1
    @inbounds for k in 2:n
        if v[k] > v[im]
            im = k
        end
    end
    return im
end
