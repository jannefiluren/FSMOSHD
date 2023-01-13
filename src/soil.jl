function soil(fsm::FSM)

  @unpack TILE, tthresh = fsm

  @unpack dt = fsm

  @unpack Dzsoil, Nsoil, Nx, Ny = fsm

  @unpack Tsoil = fsm

  @unpack tilefrac = fsm

  @unpack csoil, ksoil = fsm

  @unpack Gsoil = fsm

  @unpack gammasoil = fsm

  a = zeros(Nsoil)
  b = zeros(Nsoil)
  c = zeros(Nsoil)
  dTs = zeros(Nsoil)
  Gs = zeros(Nsoil)
  rhs = zeros(Nsoil)

  for j = 1:Ny
    for i = 1:Nx

      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        for k = 1:Nsoil-1
          Gs[k] = 2 / (Dzsoil[k] / ksoil[k, i, j] + Dzsoil[k+1] / ksoil[k+1, i, j])
        end
        a[1] = 0
        b[1] = csoil[1, i, j] + Gs[1] * dt
        c[1] = -Gs[1] * dt
        rhs[1] = (Gsoil[i, j] - Gs[1] * (Tsoil[1, i, j] - Tsoil[2, i, j])) * dt
        for k = 2:Nsoil-1
          a[k] = c[k-1]
          b[k] = csoil[k, i, j] + (Gs[k-1] + Gs[k]) * dt
          c[k] = -Gs[k] * dt
          rhs[k] = Gs[k-1] * (Tsoil[k-1, i, j] - Tsoil[k, i, j]) * dt + Gs[k] * (Tsoil[k+1, i, j] - Tsoil[k, i, j]) * dt
        end
        k = Nsoil
        Gs[k] = ksoil[k, i, j] / Dzsoil[k]
        a[k] = c[k-1]
        b[k] = csoil[k, i, j] + (Gs[k-1] + Gs[k]) * dt
        c[k] = 0
        rhs[k] = Gs[k-1] * (Tsoil[k-1, i, j] - Tsoil[k, i, j]) * dt
        tridiag!(dTs, Nsoil, gammasoil, Nsoil, a, b, c, rhs)
        ###call TRIDIAG(Nsoil,Nsoil,a,b,c,rhs,dTs)
        for k = 1:Nsoil
          Tsoil[k, i, j] = Tsoil[k, i, j] + dTs[k]
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