"""
    soil!(fsm)

Soil thermal processes and heat conduction calculations.

The per-cell physics lives in `soil_kernel!`, a KernelAbstractions kernel
launched over the whole grid (see `ebalsrf!` for the pattern). Each cell
solves its own tridiagonal system using kernel-local `MVector` scratch, so
the routine is thread-safe per cell (the former shared scratch vectors in
`FSM` are no longer used). The number of soil layers is passed as
`Val(Nsoil)` because the scratch size must be known at compile time.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
"""
function soil!(fsm::FSM{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    @unpack TILE, tthresh = fsm

    @unpack dt = fsm

    @unpack Dzsoil, Nsoil, Nx, Ny = fsm

    @unpack Tsoil = fsm

    @unpack tilefrac = fsm

    @unpack csoil, ksoil = fsm

    @unpack Gsoil = fsm

    # Strings cannot cross into kernels: resolve the tile test here
    glacier_tile = TILE == "glacier"

    backend = get_backend(Tsoil)
    kernel! = soil_kernel!(backend)
    kernel!(
        Tsoil,
        Dzsoil, tilefrac, csoil, ksoil, Gsoil,
        dt, tthresh, glacier_tile, Val(Int(Nsoil));
        ndrange = (Int(Nx), Int(Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

# inbounds = true keeps the kernel-local MVector scratch off the heap (the
# bounds-check error paths would otherwise capture it, forcing a heap
# allocation per grid cell). NOTE: a raw @inbounds block written directly in
# a @kernel body must NOT be used instead - it interferes with the
# KernelAbstractions CPU code transformation and silently corrupts results.
# Tests run with --check-bounds=yes, which overrides this, so all indexing
# stays validated in CI.
@kernel inbounds = true function soil_kernel!(
        Tsoil,
        Dzsoil, tilefrac,
        csoil, ksoil, Gsoil,
        dt::Tf, tthresh::Tf, glacier_tile::Bool,
        ::Val{Nsoil},
    ) where {Tf, Nsoil}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        # Kernel-local scratch (one set per grid cell)
        asoil = zero(MVector{Nsoil, Tf})
        bsoil = zero(MVector{Nsoil, Tf})
        cssoil = zero(MVector{Nsoil, Tf})
        dTssoil = zero(MVector{Nsoil, Tf})
        Gssoil = zero(MVector{Nsoil, Tf})
        rhssoil = zero(MVector{Nsoil, Tf})
        gammasoil = zero(MVector{Nsoil, Tf})

        for k in 1:(Nsoil - 1)
            Gssoil[k] = Tf(2) / (Dzsoil[k] / ksoil[k, i, j] + Dzsoil[k + 1] / ksoil[k + 1, i, j])
        end
        asoil[1] = Tf(0)
        bsoil[1] = csoil[1, i, j] + Gssoil[1] * dt
        cssoil[1] = -Gssoil[1] * dt
        rhssoil[1] = (Gsoil[i, j] - Gssoil[1] * (Tsoil[1, i, j] - Tsoil[2, i, j])) * dt
        for k in 2:(Nsoil - 1)
            asoil[k] = cssoil[k - 1]
            bsoil[k] = csoil[k, i, j] + (Gssoil[k - 1] + Gssoil[k]) * dt
            cssoil[k] = -Gssoil[k] * dt
            rhssoil[k] = Gssoil[k - 1] * (Tsoil[k - 1, i, j] - Tsoil[k, i, j]) * dt + Gssoil[k] * (Tsoil[k + 1, i, j] - Tsoil[k, i, j]) * dt
        end
        k = Nsoil
        Gssoil[k] = ksoil[k, i, j] / Dzsoil[k]
        asoil[k] = cssoil[k - 1]
        bsoil[k] = csoil[k, i, j] + (Gssoil[k - 1] + Gssoil[k]) * dt
        cssoil[k] = Tf(0)
        rhssoil[k] = Gssoil[k - 1] * (Tsoil[k - 1, i, j] - Tsoil[k, i, j]) * dt
        tridiag!(dTssoil, Nsoil, gammasoil, Nsoil, asoil, bsoil, cssoil, rhssoil)
        for k in 1:Nsoil
            Tsoil[k, i, j] = Tsoil[k, i, j] + dTssoil[k]
        end

        # Cap glacier temperatures to 0°C
        # This does not conserve energy.
        # The excess energy would correspond to glacier melting, which we don't track.
        if glacier_tile
            for k in 1:Nsoil
                Tsoil[k, i, j] = min(Tsoil[k, i, j], Tm)
            end
        end

    end
end
