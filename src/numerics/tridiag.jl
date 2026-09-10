"""
    tridiag!(x, Nvec, a, b, c, r)

Tridiagonal matrix solver using the Thomas algorithm: solve the `Nvec`-equation system with
sub/main/super-diagonals `a`/`b`/`c` and right-hand side `r`, writing the result into `x`.

# Arguments
- `x::AbstractVector`: Solution vector (output)
- `Nvec`: Number of equations in the system (`<= length(x)`)
- `a`: Sub-diagonal coefficients
- `b`: Main diagonal coefficients
- `c`: Super-diagonal coefficients
- `r`: Right-hand side vector
"""
# Base.@propagate_inbounds (not @inline): carries the caller kernel's `inbounds = true` into this
# function so the internally-allocated `gamma` MVector stays on the stack rather than heap-allocating
# once per grid cell. Callers guarantee Nvec <= length(x); the @inbounds block relies on it.
Base.@propagate_inbounds function tridiag!(x::AbstractVector{Tf}, Nvec, a, b, c, r) where {Tf <: Real}

    # Elimination-coefficient workspace, sized from x (a stack MVector when x is one).
    gamma = similar(x)

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
