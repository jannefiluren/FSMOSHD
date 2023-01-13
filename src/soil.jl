function soil(fsm::FSM)

  @unpack TILE, tthresh = fsm

  @unpack dt = fsm

  @unpack Dzsoil, Nsoil, Nx, Ny = fsm

  @unpack Tsoil = fsm

  @unpack tilefrac = fsm

  @unpack csoil, ksoil = fsm

  @unpack Gsoil = fsm

  @unpack asoil, bsoil, cssoil, dTssoil, Gssoil, rhssoil = fsm

  @unpack gammasoil = fsm

  asoil .= 0.0
  bsoil .= 0.0
  cssoil .= 0.0
  dTssoil .= 0.0
  Gssoil .= 0.0
  rhssoil .= 0.0

  for j = 1:Ny
    for i = 1:Nx

      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        for k = 1:Nsoil-1
          Gssoil[k] = 2 / (Dzsoil[k] / ksoil[k, i, j] + Dzsoil[k+1] / ksoil[k+1, i, j])
        end
        asoil[1] = 0
        bsoil[1] = csoil[1, i, j] + Gssoil[1] * dt
        cssoil[1] = -Gssoil[1] * dt
        rhssoil[1] = (Gsoil[i, j] - Gssoil[1] * (Tsoil[1, i, j] - Tsoil[2, i, j])) * dt
        for k = 2:Nsoil-1
          asoil[k] = cssoil[k-1]
          bsoil[k] = csoil[k, i, j] + (Gssoil[k-1] + Gssoil[k]) * dt
          cssoil[k] = -Gssoil[k] * dt
          rhssoil[k] = Gssoil[k-1] * (Tsoil[k-1, i, j] - Tsoil[k, i, j]) * dt + Gssoil[k] * (Tsoil[k+1, i, j] - Tsoil[k, i, j]) * dt
        end
        k = Nsoil
        Gssoil[k] = ksoil[k, i, j] / Dzsoil[k]
        asoil[k] = cssoil[k-1]
        bsoil[k] = csoil[k, i, j] + (Gssoil[k-1] + Gssoil[k]) * dt
        cssoil[k] = 0
        rhssoil[k] = Gssoil[k-1] * (Tsoil[k-1, i, j] - Tsoil[k, i, j]) * dt
        tridiag!(dTssoil, Nsoil, gammasoil, Nsoil, asoil, bsoil, cssoil, rhssoil)
        ###call TRIDIAG(Nsoil,Nsoil,asoil,bsoil,cssoil,rhssoil,dTssoil)
        for k = 1:Nsoil
          Tsoil[k, i, j] = Tsoil[k, i, j] + dTssoil[k]
        end

        # Cap glacier temperatures to 0°C
        # This does not conserve energy.
        # The excess energy would correspond to glacier melting, which we don't track.
        if (TILE == "glacier")
          for k = 1:Nsoil
            Tsoil[k, i, j] = min(Tsoil[k, i, j], Tm)
          end
        end

      end

    end
  end

end