"""
    canopy!(fsm, meteo)

Snow interception, sublimation, and unloading from vegetation canopy.

The per-cell physics lives in `canopy_kernel!`, a KernelAbstractions kernel
launched over the whole grid (see `ebalsrf!` for the pattern).

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions (read-only)
"""
function canopy!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    @unpack tthresh = fsm

    @unpack Sfeff = fsm

    @unpack Nx, Ny, dt = fsm

    @unpack tcnc, tcnm, psf, psr = fsm

    @unpack scap = fsm

    @unpack Sveg, Tveg = fsm

    @unpack fveg, pmultf, tilefrac = fsm

    @unpack Eveg = fsm

    @unpack intcpt, Sbveg, unload = fsm

    backend = get_backend(Sveg)
    kernel! = canopy_kernel!(backend)
    kernel!(
        unload, intcpt, Sbveg, Sveg, Sfeff,
        scap, Tveg, fveg, pmultf, tilefrac, Eveg,
        dt, tthresh, tcnc, tcnm, psf, psr;
        ndrange = (Int(Nx), Int(Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function canopy_kernel!(
        unload, intcpt, Sbveg, Sveg, Sfeff,
        scap, Tveg, fveg, pmultf,
        tilefrac, Eveg,
        dt::Tf, tthresh::Tf, tcnc::Tf, tcnm::Tf, psf::Tf, psr::Tf,
    ) where {Tf}

    i, j = @index(Global, NTuple)

    @unpack_constants(Tf)

    unload[i, j] = Tf(0)
    intcpt[i, j] = Tf(0)
    Sbveg[i, j] = Tf(0)

    if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        # rescale precipitation to correct back precip multiplier applied to open area
        Sfeff[i, j] = pmultf[i, j] * Sfeff[i, j]

        # interception
        intcpt[i, j] = (scap[i, j] - Sveg[i, j]) * (Tf(1) - exp(-fveg[i, j] * Sfeff[i, j] * dt / scap[i, j]))
        Sveg[i, j] = Sveg[i, j] + intcpt[i, j]
        Sfeff[i, j] = Sfeff[i, j] - intcpt[i, j] / dt
        Sfeff[i, j] = (psf - psr * fveg[i, j]) * Sfeff[i, j] # including preferential deposition in canopy gaps; might have to be revisited to ensure mass conservation, potentially integrate with pmultf

        # sublimation
        Evegs = Tf(0)
        if (Sveg[i, j] > eps(Tf) || Tveg[i, j] < Tm)
            Evegs = Eveg[i, j]
        end
        Sveg[i, j] = Sveg[i, j] - Evegs * dt
        Sbveg[i, j] = Evegs * dt
        if (Sveg[i, j] < Tf(0))
            Sbveg[i, j] = Sbveg[i, j] + Sveg[i, j]
        end
        Sveg[i, j] = max(Sveg[i, j], Tf(0))

        # unloading
        tunl = tcnc
        if (Tveg[i, j] >= Tm)
            tunl = tcnm
        end
        tunl = max(tunl, dt)
        unload[i, j] = Sveg[i, j] * dt / tunl
        Sveg[i, j] = Sveg[i, j] - unload[i, j]


    end
end
