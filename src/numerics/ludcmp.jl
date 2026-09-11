"""
    ludcmp!(N, A, b, x)

Solve the `N`-by-`N` linear system `A x = b` by LU decomposition with partial pivoting, writing
the solution into `x`. `A` is left unchanged; the decomposition and pivoting workspace is
allocated internally.

Numerical Recipes' `ludcmp`/`lubksb`: Crout ordering with implicit (scaled) partial pivoting; the
factors are Doolittle-form (unit-diagonal `L`, `U` carries the diagonal).
"""
# @propagate_inbounds (not @inline): the kernel's inbounds context reaches the internal scratch, keeping it off the heap.
Base.@propagate_inbounds function ludcmp!(N::Integer, A::AbstractMatrix{Tf}, b::AbstractVector{Tf}, x::AbstractVector{Tf}) where {Tf <: Real}

    Acp = similar(A)
    vv = similar(x)
    indx = similar(x, Int32)    # TODO check if this should be Int...

    @inbounds begin

        Acp .= A
        x .= b
        vv .= 0
        indx .= 0

        # Scaling
        for i in 1:N
            aamax = abs(Acp[i, 1])
            for k in 2:N
                aamax = max(aamax, abs(Acp[i, k]))
            end
            vv[i] = 1 / aamax
        end

        # LU decomposition with partial pivoting
        for j in 1:N
            for i in 1:(j - 1)
                sum = Acp[i, j]
                if i > 1
                    for k in 1:(i - 1)
                        sum -= Acp[i, k] * Acp[k, j]
                    end
                end
                Acp[i, j] = sum
            end
            aamax = 0.0
            imax = j
            for i in j:N
                sum = Acp[i, j]
                for k in 1:(j - 1)
                    sum -= Acp[i, k] * Acp[k, j]
                end
                Acp[i, j] = sum
                dum = vv[i] * abs(sum)
                if dum >= aamax
                    imax = i
                    aamax = dum
                end
            end
            if (j != imax)
                for k in 1:N
                    dum = Acp[imax, k]
                    Acp[imax, k] = Acp[j, k]
                    Acp[j, k] = dum
                end
                vv[imax] = vv[j]
            end
            indx[j] = imax
            if Acp[j, j] == 0.0
                Acp[j, j] = 1.0e-20
            end
            if j != N
                dum = 1 / Acp[j, j]
                for i in (j + 1):N
                    Acp[i, j] *= dum
                end
            end
        end

        # Forward substitution
        ii = 0
        for i in 1:N
            ll = indx[i]
            sum = x[ll]
            x[ll] = x[i]
            if ii != 0
                for j in ii:(i - 1)
                    sum -= Acp[i, j] * x[j]
                end
            elseif sum != 0.0
                ii = i
            end
            x[i] = sum
        end

        # Backward substitution
        for i in N:-1:1
            sum = x[i]
            for j in (i + 1):N
                sum -= Acp[i, j] * x[j]
            end
            x[i] = sum / Acp[i, i]
        end

    end

    return nothing
end
