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

function compute_psolid!(Sf, ptot, ta, thres_prec=274.19, m_prec=0.1500)

  p_corr = 1.0
  Sf .= ptot .* p_corr ./ (1 .+ exp.((ta .- thres_prec) ./ m_prec))

end

function compute_pliquid!(Rf, ptot, ta, thres_prec=274.19, m_prec=0.1500)

  p_corr = 1.0 
  Rf .= ptot .* p_corr .* exp.((ta .- thres_prec) ./ m_prec) ./ (1 .+ exp.((ta .- thres_prec) ./ m_prec))

end

function drive_grid!(meteo::MET, fsm::FSM, t::DateTime)

  folder = joinpath("K:/DATA_COSMO/OUTPUT_GRID_OSHD_0250/PROCESSED_ANALYSIS/COSMO_1EFA", Dates.format(t, "yyyy.mm"))
  filename = searchdir(folder, "COSMODATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_")

  meteo_in = matread(joinpath(folder, filename[1]))

  meteo.Ta = meteo_in["tais"]["data"]
  meteo.RH = meteo_in["rhus"]["data"]
  meteo.Ua = meteo_in["wnsc"]["data"]
  meteo.Ua[meteo.Ua .< 0.1] .= 0.1
  meteo.Sdir = meteo_in["sdri"]["data"]
  meteo.Sdif = meteo_in["sdfd"]["data"]
  meteo.LW = meteo_in["lwrc"]["data"]
  ptot = meteo_in["prcs"]["data"] ./ fsm.dt
  compute_psolid!(meteo.Sf, ptot, meteo.Ta)
  compute_pliquid!(meteo.Rf, ptot, meteo.Ta)
  meteo.Tc .= meteo.Ta .- Tm
  meteo.Ps = meteo_in["pail"]["data"]

  meteo.es .= e0 .* exp.(17.5043 .* meteo.Tc ./ (241.3 .+ meteo.Tc))
  meteo.Qa .= (meteo.RH ./ 100) .* eps_fsm .* meteo.es ./ meteo.Ps

  #computation of Sf24h assuming hourly input
  curr_hour = Dates.value(Hour(t)) + 1 #0h -> 1; 23h -> 24
  meteo.Sf24h .+= meteo.Sf .* fsm.dt
  meteo.Sf24h .-= meteo.Sf_history[:,:,curr_hour] .* fsm.dt
  meteo.Sf_history[:,:,curr_hour] = meteo.Sf

  #TODO: missing: Tv

end
