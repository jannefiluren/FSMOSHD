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
