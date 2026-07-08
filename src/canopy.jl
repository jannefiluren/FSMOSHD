"""
    canopy!(fsm, meteo)

Snow interception, sublimation, and unloading from vegetation canopy.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
- `meteo::MET`: Current meteorological conditions (may be modified)
"""
function canopy!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    @unpack_constants(Tf)

    @unpack tthresh = fsm

    @unpack Sf = meteo

    @unpack Nx, Ny, dt = fsm

    @unpack tcnc, tcnm, psf, psr = fsm

    @unpack scap = fsm

    @unpack Sveg, Tveg = fsm

    @unpack fveg, pmultf, tilefrac = fsm

    @unpack Eveg = fsm

    @unpack intcpt, Sbveg, unload = fsm

    for j in 1:Ny
        for i in 1:Nx

            unload[i, j] = Tf(0)
            intcpt[i, j] = Tf(0)
            Sbveg[i, j] = Tf(0)

            if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

                if (fveg[i, j] > eps(Tf))
                    # rescale precipitation to correct back precip multiplier applied to open area
                    Sf[i, j] = pmultf[i, j] * Sf[i, j]

                    # interception
                    intcpt[i, j] = (scap[i, j] - Sveg[i, j]) * (Tf(1) - exp(-fveg[i, j] * Sf[i, j] * dt / scap[i, j]))
                    Sveg[i, j] = Sveg[i, j] + intcpt[i, j]
                    Sf[i, j] = Sf[i, j] - intcpt[i, j] / dt
                    Sf[i, j] = (psf - psr * fveg[i, j]) * Sf[i, j] # including preferential deposition in canopy gaps; might have to be revisited to ensure mass conservation, potentially integrate with pmultf

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

        end
    end

    return nothing
end
