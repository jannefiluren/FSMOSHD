"""
    tridiag!(x, Nvec, gamma, Nmax, a, b, c, r)

Tridiagonal matrix solver using Thomas algorithm.

# Arguments
- `x::Vector`: Solution vector (output)
- `Nvec`: Number of equations in the system
- `gamma`: Workspace vector for elimination coefficients
- `Nmax`: Maximum system size (for array bounds)
- `a`: Sub-diagonal coefficients
- `b`: Main diagonal coefficients  
- `c`: Super-diagonal coefficients
- `r`: Right-hand side vector
"""
# @inline so that kernel-local MVector arguments do not escape (escaping
# would force them onto the heap, allocating once per grid cell)
@inline function tridiag!(x::AbstractVector{Tf}, Nvec, gamma, Nmax, a, b, c, r) where {Tf <: Real}

    # @inbounds (callers guarantee Nvec <= length): without it, the
    # bounds-check error paths would capture the kernel-local MVector
    # arguments and force them onto the heap
    @inbounds begin
        fill!(gamma, zero(Tf))

        beta = b[1]
        x[1] = r[1] / beta

        for n in 2:Nvec
            gamma[n] = c[n - 1] / beta
            beta = b[n] - a[n] * gamma[n]
            x[n] = (r[n] - a[n] * x[n - 1]) / beta
        end

        for n in (Nvec - 1):-1:1
            x[n] = x[n] - gamma[n + 1] * x[n + 1]
        end
    end

    return nothing
end
