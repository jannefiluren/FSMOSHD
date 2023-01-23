function drive(fsm, data)

  @unpack dt = fsm

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


function drive!(fsm::FSM, meteo::MET)

  @unpack dt = fsm

  @unpack Tc, es, Qa, Ua, Sf, Rf, Ta, RH, Ps = meteo

  Ua .= max.(Ua, 0.1)

  Sf .= Sf ./ dt
  Rf .= Rf ./ dt
  Tc .= Ta .- Tm
  es .= e0 .* exp.(17.5043 .* Tc ./ (241.3 .+ Tc))
  Qa .= (RH ./ 100) .* eps_fsm .* es ./ Ps

end

function compute_psolid(ptot, ta, thres_prec=274.19, m_prec=0.1500)

  p_corr = 1.0
  tp = @. (ta - thres_prec) / m_prec
  p_multi = @. p_corr / (1 + exp(tp))

  return @. p_multi * ptot

end

function compute_pliquid(ptot, ta, thres_prec=274.19, m_prec=0.1500)

  p_corr = 1.0

  Tp = @. (ta - thres_prec) / m_prec
  p_multi = @. p_corr * exp(Tp) / (1 + exp(Tp))
  return @. p_multi * ptot

end

function read_meteo!(Sdir, Sdif, LW, Sf, Rf, Ta, RH, Ua, Ps, Sf24h, year, month, day, hour)

  folder = joinpath("K:/DATA_COSMO/OUTPUT_GRID_OSHD_0250/PROCESSED_ANALYSIS/COSMO_1EFA", Dates.format(t, "yyyy.mm"))
  filename = searchdir(folder, "COSMODATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_")

  meteo = matread(joinpath(folder, filename[1]))

  Ta = meteo["tais"]["data"]
  RH = meteo["rhus"]["data"]
  Ua = meteo["wnsc"]["data"]
  SW = meteo["sdrd"]["data"] .+ meteo["sdfd"]["data"]
  LW = meteo["lwrc"]["data"]
  Sf = compute_psolid(meteo["prcs"]["data"], meteo["tais"]["data"]) ./ 3600
  Rf = compute_pliquid(meteo["prcs"]["data"], meteo["tais"]["data"]) ./ 3600
  Ps = meteo["pail"]["data"]

  return Ta, RH, Ua, SW, LW, Sf, Rf, Ps

end
