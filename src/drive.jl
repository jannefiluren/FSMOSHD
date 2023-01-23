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

searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))

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

function drive_grid!(meteo::MET, fsm::FSM, t::DateTime)

  folder = joinpath("home", "haugened", "Documents", "data", "FSM_input", "grid", "DATA_COSMO", "OUTPUT_GRID_OSHD_0250", "PROCESSED_ANALYSIS", "COSMO_1EFA", Dates.format(t, "yyyy.mm"))
  filename = searchdir(folder, "COSMODATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_")

  meteo_in = matread(joinpath(folder, filename[1]))

  meteo.Ta = meteo_in["tais"]["data"]
  meteo.RH = meteo_in["rhus"]["data"]
  meteo.Ua = maximum.(meteo_in["wnsc"]["data"], 0.1)
  meteo.Sdir = meteo_in["sdrd"]["data"]
  meteo.Sdif = meteo_in["sdfd"]["data"]
  meteo.LW = meteo_in["lwrc"]["data"]
  meteo.Sf = compute_psolid(meteo_in["prcs"]["data"], meteo_in["tais"]["data"]) ./ fsm.dt
  meteo.Rf = compute_pliquid(meteo_in["prcs"]["data"], meteo_in["tais"]["data"]) ./ fsm.dt
  meteo.Tc .= meteo.Ta .- Tm
  meteo.Ps = meteo_in["pail"]["data"]

  meteo.es .= e0 .* exp.(17.5043 .* meteo.Tc ./ (241.3 .+ meteo.Tc))
  meteo.Qa .= (meteo.RH ./ 100) .* eps_fsm .* meteo.es ./ meteo.Ps

  #computation of Sf24h assuming hourly input
  curr_hour = Dates.value(Hour(t))
  Sf_sum .-= Sf_history[:,:,curr_hour]
  Sf_sum .+= meteo.Sf
  Sf_history[:,:,curr_hour] = meteo.Sf

  #TODO: missing: Tv

end
