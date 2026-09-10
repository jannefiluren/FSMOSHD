# Canopy mass-balance parameterizations. radiation.jl and surface_energy_balance.jl add further
# AbstractCanopy-dispatched methods (solar_radiation!/thermal_radiation!, energy_balance!).

struct NoCanopy{Tf} <: AbstractCanopy{Tf} end

@kwdef struct OneLayerCanopy{Tf} <: AbstractCanopy{Tf}
    fsar::Tf = 0.1                       # Albedo adjustment range vs vegetation fraction (-)
    avg0::Tf = 0.1                       # Snow-free vegetation albedo (-)
    avgs::Tf = 0.4                       # Snow-covered vegetation albedo (-)
    psf::Tf = 1                          # Solid precipitation multiplier at min canopy cover (-)
    psr::Tf = 0.1                        # Solid precipitation multiplier range (-)
    tcnc::Tf = 3600 * 240                # Canopy unloading time scale for cold snow (s)
    tcnm::Tf = 3600 * 48                 # Canopy unloading time scale for melting snow (s)
end

NoCanopy{Tf}(grid::Grid; kwargs...) where {Tf} = NoCanopy{Tf}()
OneLayerCanopy{Tf}(grid::Grid; kwargs...) where {Tf} = OneLayerCanopy{Tf}(; kwargs...)

canopy_fsar(c::OneLayerCanopy) = c.fsar
canopy_avg0(c::OneLayerCanopy) = c.avg0
canopy_avgs(c::OneLayerCanopy) = c.avgs

"""
    canopy_snow!(canopy, i, j, state, diag, surface, params)

Snow on the canopy at cell `(i, j)`: interception from the throughfall `diag.Sfeff`,
sublimation and unloading, updating `state.Sveg` and `diag.intcpt`/`Sbveg`/`unload`.
"""
function canopy_snow! end

@inline canopy_snow!(::NoCanopy, i, j, state, diag, surface, params) = nothing

@inline function canopy_snow!(canopy::OneLayerCanopy{Tf}, i, j, state, diag, surface, params) where {Tf}

    @unpack_constants(Tf)

    (; dt) = params
    (; tcnc, tcnm) = canopy
    (; scap, fveg, pmultf) = surface
    (; Sveg, Tveg) = state
    (; unload, intcpt, Sbveg, Sfeff, Eveg) = diag

    unload[i, j] = Tf(0)
    intcpt[i, j] = Tf(0)
    Sbveg[i, j] = Tf(0)

    # Remove precipitation scaling applied to forcing data
    Sfeff[i, j] = pmultf[i, j] * Sfeff[i, j]

    # Interception of snow on canopies
    intcpt[i, j] = (scap[i, j] - Sveg[i, j]) * (Tf(1) - exp(-fveg[i, j] * Sfeff[i, j] * dt / scap[i, j]))
    Sveg[i, j] = Sveg[i, j] + intcpt[i, j]
    Sfeff[i, j] = Sfeff[i, j] - intcpt[i, j] / dt

    # Preferential deposition of snowfall in canopy gaps (not mass conserving)
    Sfeff[i, j] = (canopy.psf - canopy.psr * fveg[i, j]) * Sfeff[i, j]

    # Sublimation of intercepted snow
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

    # Unloading of intercepted snow
    tunl = tcnc
    if (Tveg[i, j] >= Tm)
        tunl = tcnm
    end
    tunl = max(tunl, dt)
    unload[i, j] = Sveg[i, j] * dt / tunl
    Sveg[i, j] = Sveg[i, j] - unload[i, j]

    return nothing
end
