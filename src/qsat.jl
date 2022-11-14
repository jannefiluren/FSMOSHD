function qsat(P, T)

  eps_fsm = 0.622   # Ratio of molecular weights of water and dry air
  Tm = 273.15       # Melting point (K)
  e0 = 610.78       # Saturation vapour pressure at Tm (Pa)

  Tc = T - Tm
  if (Tc > 0)
    es = e0 * exp(17.5043 * Tc / (241.3 + Tc))
  else
    es = e0 * exp(22.4422 * Tc / (272.186 + Tc))
  end
  Qs = eps_fsm * es / P

  return Qs

end