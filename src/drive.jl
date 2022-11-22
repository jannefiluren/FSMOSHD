function drive(data)

  tmp = parse.(Float64, split(data))

  year = tmp[1]
  month = tmp[2]
  day = tmp[3]
  hour = tmp[4]
  Sdir[:, :] .= tmp[5]
  Sdif[:, :] .= tmp[6]
  LW[:, :] .= tmp[7]
  Sf[:, :] .= tmp[8]
  Rf[:, :] .= tmp[9]
  Ta[:, :] .= tmp[10]
  RH[:, :] .= tmp[11]
  Ua[:, :] .= tmp[12]
  Ps[:, :] .= tmp[13]
  Sf24h[:, :] .= tmp[14]

  Tc = similar(Ta)   #### hack
  es = similar(Ta)   #### hack
  # Qa = similar(Ta)   #### hack

  Ua .= max.(Ua, 0.1)

  Sf .= Sf ./ dt
  Rf .= Rf ./ dt
  Tc .= Ta .- Tm
  es .= e0 * exp.(17.5043 * Tc ./ (241.3 .+ Tc))
  Qa .= (RH ./ 100) .* eps_fsm .* es ./ Ps

  return year, month, day, hour

end
