"""
    soil!(fsm)

Soil thermal processes and heat conduction calculations.

# Arguments
- `fsm::FSM`: Model state structure (modified in-place)
"""
function soil!(fsm::FSM{Tf, Ti}) where {Tf <: Real, Ti <: Integer}

    @unpack_constants(Tf)

    @unpack TILE, tthresh, glacierfrac = fsm

    @unpack dt = fsm

    @unpack Dzsoil, Nsoil, Nx, Ny = fsm

    @unpack Tsoil = fsm

    @unpack tilefrac = fsm

    @unpack csoil, ksoil = fsm

    @unpack Gsoil = fsm

    @unpack asoil, bsoil, cssoil, dTssoil, Gssoil, rhssoil = fsm

    @unpack gammasoil = fsm

    asoil .= Tf(0)
    bsoil .= Tf(0)
    cssoil .= Tf(0)
    dTssoil .= Tf(0)
    Gssoil .= Tf(0)
    rhssoil .= Tf(0)

    for j in 1:Ny
        for i in 1:Nx

            if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

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
                if (TILE == "glacier")
                    for k in 1:Nsoil
                        Tsoil[k, i, j] = min(Tsoil[k, i, j], Tm)
                    end
                end

            end

        end
    end

    return nothing
end
