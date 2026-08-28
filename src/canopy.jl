# ---------------------------------------------------------------------------
# Canopy parameterizations
#
# Whether a tile has canopy is a property of the tile, not of a cell: after the
# Stage 4 mask narrowing, every active cell of a canopy tile has fveg > 0 and
# every active cell of a NoCanopy tile has fveg == 0. That is what makes this a
# dispatch axis rather than a per-cell branch.
#
# Both schemes hold scalars only, so they are isbits and cross into a kernel by
# value.
# ---------------------------------------------------------------------------

struct NoCanopy{Tf} <: AbstractCanopy{Tf} end

@kwdef struct OneLayerCanopy{Tf} <: AbstractCanopy{Tf}
    fsar::Tf = 0.1                       # Albedo adjustment range vs vegetation fraction (-)
    avg0::Tf = 0.1                       # Snow-free vegetation albedo (-)
    avgs::Tf = 0.4                       # Snow-covered vegetation albedo (-)
    psf::Tf = 1                          # Solid precipitation multiplier at min canopy cover (-)
    psr::Tf = 0.1                        # Solid precipitation multiplier range (-)
end

NoCanopy{Tf}(Nx, Ny; kwargs...) where {Tf} = NoCanopy{Tf}()
OneLayerCanopy{Tf}(Nx, Ny; kwargs...) where {Tf} = OneLayerCanopy{Tf}(; kwargs...)


# Neutral values let the shared radiation prologue run without branching: with no
# canopy, fveg == 0 makes fsar's contribution vanish and aveg is multiplied by
# acan == 0, so any finite value is correct.
canopy_fsar(c::OneLayerCanopy) = c.fsar
canopy_fsar(::NoCanopy{Tf}) where {Tf} = zero(Tf)
canopy_avg0(c::OneLayerCanopy) = c.avg0
canopy_avg0(::NoCanopy{Tf}) where {Tf} = zero(Tf)
canopy_avgs(c::OneLayerCanopy) = c.avgs
canopy_avgs(::NoCanopy{Tf}) where {Tf} = zero(Tf)

"""
    surface_balance!(canopy, fsm, met)

Solve the surface energy balance. `NoCanopy` uses the surface-only solver; `OneLayerCanopy`
uses the joint surface+canopy solver. Replaces the `TILE == "forest"` test in `step!`.
"""
surface_balance!(::NoCanopy, fsm, met) = ebalsrf!(fsm, met)
surface_balance!(::OneLayerCanopy, fsm, met) = ebalfor!(fsm, met)

"""
    canopy!(canopy, fsm, met)

Canopy interception, sublimation and unloading. A no-op without canopy.
"""
canopy!(::NoCanopy, fsm, met) = nothing

"""
    canopy!(fsm, meteo)

Snow interception, sublimation, and unloading from vegetation canopy.

The per-cell physics lives in `canopy_kernel!`, a KernelAbstractions kernel
launched over the whole grid (see `ebalsrf!` for the pattern).

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions (read-only)
"""
function canopy!(::OneLayerCanopy, fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    @unpack tthresh = fsm

    @unpack Sfeff = fsm

    @unpack Nx, Ny, dt = fsm

    @unpack tcnc, tcnm, CANOPY = fsm

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
        dt, tthresh, tcnc, tcnm, CANOPY;
        ndrange = (Int(Nx), Int(Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function canopy_kernel!(
        unload, intcpt, Sbveg, Sveg, Sfeff,
        scap, Tveg, fveg, pmultf,
        tilefrac, Eveg,
        dt::Tf, tthresh::Tf, tcnc::Tf, tcnm::Tf, CANOPY::OneLayerCanopy{Tf},
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
        Sfeff[i, j] = (CANOPY.psf - CANOPY.psr * fveg[i, j]) * Sfeff[i, j] # including preferential deposition in canopy gaps; might have to be revisited to ensure mass conservation, potentially integrate with pmultf

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
