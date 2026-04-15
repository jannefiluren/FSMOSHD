function fresh_snow_density!(fsm::FSM{Tf, Ti}, Ta, Ua, dem) where {Tf<:Real, Ti<:Integer}
  
  @unpack_constants(Tf)

  @unpack rho0, rhob, rhoc, rhof, rhos_min = fsm

  @unpack FSNRHO = fsm

  if (FSNRHO == 0)
    # Fixed fresh snow density
    rhonew = rho0
  elseif (FSNRHO == 1)
    # Climate-dependent fresh snow density
    rhonew = max(rhof + rhob * (Ta - Tm) + rhoc * Ua^Tf(0.5), rhos_min)
  else # FSNRHO == 2
    # Climate-dependent with elevation-dependent decompaction
    rhonew = rhof + rhob * (Ta - Tm) + rhoc * Ua^Tf(0.5)
    if (dem <= Tf(1000))
      t_decompaction = Tf(24.0)
    elseif (dem > Tf(4000))
      t_decompaction = Tf(0.0)
    else
      t_decompaction = Tf(24) + (dem - Tf(1000)) / (Tf(4000) - Tf(1000)) * (Tf(0) - Tf(24))
    end
    rhonew = Tf(300) + (rhonew - Tf(300)) * exp(t_decompaction / Tf(100))
    rhonew = max(rhonew, rhos_min)
  end
  
  return rhonew
  
end