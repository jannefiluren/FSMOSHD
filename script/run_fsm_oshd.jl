using Pkg
Pkg.activate("../.")

using FSMOSHD
using Parameters
using Dates
using DelimitedFiles

drive_data = readdlm(drive_file)

fsm = FSM{Float64}()
setup_point!(fsm, terrain_file, state_file=state_file)

meteo = MET{Float64}()

include("../src/parameters.jl")
include("../src/initialize.jl")
include("../src/qsat.jl")
include("../src/tridiag.jl")
include("../src/drive.jl")
include("../src/radiation.jl")
include("../src/thermal.jl")
include("../src/sfexch.jl")
include("../src/ebalsrf.jl")
include("../src/snow.jl")
include("../src/soil.jl")

fortran = readlines("../fortran/temp/test_setup.txt")

albs_fortran = parse.(Float64, split(fortran[1]))
Ds_fortran = parse.(Float64, split(fortran[2]))
Nsnow_fortran = parse.(Float64, split(fortran[3]))
Qcan_fortran = parse.(Float64, split(fortran[4]))
Sice_fortran = parse.(Float64, split(fortran[5]))
Sliq_fortran = parse.(Float64, split(fortran[6]))
Sveg_fortran = parse.(Float64, split(fortran[7]))
Tcan_fortran = parse.(Float64, split(fortran[8]))
theta_fortran = parse.(Float64, split(fortran[9]))
Tsnow_fortran = parse.(Float64, split(fortran[10]))
Tsoil_fortran = parse.(Float64, split(fortran[11]))
Tsrf_fortran = parse.(Float64, split(fortran[12]))
fsnow_fortran = parse.(Float64, split(fortran[13]))
Tveg_fortran = parse.(Float64, split(fortran[14]))
snowdepthmin_fortran = parse.(Float64, split(fortran[15]))
snowdepthmax_fortran = parse.(Float64, split(fortran[16]))
snowdepthhist_fortran = parse.(Float64, split(fortran[17]))
swemin_fortran = parse.(Float64, split(fortran[18]))
swemax_fortran = parse.(Float64, split(fortran[19]))
swehist_fortran = parse.(Float64, split(fortran[20]))
fcly_fortran = parse.(Float64, split(fortran[21]))
b_fortran = parse.(Float64, split(fortran[22]))
hcap_soil_fortran = parse.(Float64, split(fortran[23]))
sathh_fortran = parse.(Float64, split(fortran[24]))
Vsat_fortran = parse.(Float64, split(fortran[25]))
Vcrit_fortran = parse.(Float64, split(fortran[26]))

if check_final_vals

    ### Test setup

    println("\nTest setup")

    println(maximum(abs.(albs_fortran - fsm.albs)))
    println(maximum(abs.(Ds_fortran - fsm.Ds)))
    println(maximum(abs.(Nsnow_fortran - fsm.Nsnow)))
    println(maximum(abs.(Qcan_fortran - fsm.Qcan)))
    println(maximum(abs.(Sice_fortran - fsm.Sice)))
    println(maximum(abs.(Sliq_fortran - fsm.Sliq)))
    println(maximum(abs.(Sveg_fortran - fsm.Sveg)))
    println(maximum(abs.(Tcan_fortran - fsm.Tcan)))
    println(maximum(abs.(theta_fortran - fsm.theta)))
    println(maximum(abs.(Tsnow_fortran - fsm.Tsnow)))
    println(maximum(abs.(Tsoil_fortran - fsm.Tsoil)))
    println(maximum(abs.(Tsrf_fortran - fsm.Tsrf)))
    println(maximum(abs.(fsnow_fortran - fsm.fsnow)))
    println(maximum(abs.(Tveg_fortran - fsm.Tveg)))
    println(maximum(abs.(snowdepthmin_fortran - fsm.snowdepthmin)))
    println(maximum(abs.(snowdepthmax_fortran - fsm.snowdepthmax)))
    println(maximum(abs.(snowdepthhist_fortran - fsm.snowdepthhist)))
    println(maximum(abs.(swemin_fortran - fsm.swemin)))
    println(maximum(abs.(swemax_fortran - fsm.swemax)))
    println(maximum(abs.(swehist_fortran - fsm.swehist)))
    println(maximum(abs.(fcly_fortran - fsm.fcly)))
    println(maximum(abs.(b_fortran - fsm.b)))
    println(maximum(abs.(hcap_soil_fortran - fsm.hcap_soil)))
    println(maximum(abs.(sathh_fortran - fsm.sathh)))
    println(maximum(abs.(Vsat_fortran - fsm.Vsat)))
    println(maximum(abs.(Vcrit_fortran - fsm.Vcrit)))

end


Qa = similar(Ta)

if length(output_file) > 0
  fout = open(output_file, "w") 
end

for (istep, indata) in enumerate(eachrow(drive_data))

  global fortran

  println("Time step: ", istep)

  # Time information

  global year, month, day, hour

  year = indata[1]
  month = indata[2]
  day = indata[3]
  hour = indata[4]

  # Forcing data
  t = DateTime(indata[1], indata[2], indata[3], indata[4], 00, 00)
  meteo.Sdir[:, :] .= indata[5]
  meteo.Sdif[:, :] .= indata[6]
  meteo.LW[:, :] .= indata[7]
  meteo.Sf[:, :] .= indata[8]
  meteo.Rf[:, :] .= indata[9]
  meteo.Ta[:, :] .= indata[10]
  meteo.RH[:, :] .= indata[11]
  meteo.Ua[:, :] .= indata[12]
  meteo.Ps[:, :] .= indata[13]
  meteo.Sf24h[:, :] .= indata[14]

  ### Run drive

  drive!(fsm, meteo)

  ### Run radiation

  radiation(fsm, meteo, t)

  fortran = readlines("../fortran/temp/test_radiation.txt")

  global alb_fortran
  global asrf_out_fortran
  global Sdirt_fortran
  global Sdift_fortran
  global SWveg_fortran
  global SWsrf_fortran
  global SWsci_fortran
  global LWt_fortran

  alb_fortran = parse.(Float64, split(fortran[1]))
  asrf_out_fortran = parse.(Float64, split(fortran[2]))
  Sdirt_fortran = parse.(Float64, split(fortran[3]))
  Sdift_fortran = parse.(Float64, split(fortran[4]))
  SWveg_fortran = parse.(Float64, split(fortran[5]))
  SWsrf_fortran = parse.(Float64, split(fortran[6]))
  SWsci_fortran = parse.(Float64, split(fortran[7]))
  LWt_fortran = parse.(Float64, split(fortran[8]))

  ### Run thermal

  thermal(fsm)

  fortran = readlines("../fortran/temp/test_thermal.txt")

  global Ds1_fortran
  global gs1_fortran
  global ks1_fortran
  global Ts1_fortran
  global Tveg0_fortran
  global csoil_fortran
  global ksnow_fortran
  global ksoil_fortran

  Ds1_fortran = parse.(Float64, split(fortran[1]))
  gs1_fortran = parse.(Float64, split(fortran[2]))
  ks1_fortran = parse.(Float64, split(fortran[3]))
  Ts1_fortran = parse.(Float64, split(fortran[4]))
  Tveg0_fortran = parse.(Float64, split(fortran[5]))
  csoil_fortran = parse.(Float64, split(fortran[6]))
  ksnow_fortran = parse.(Float64, split(fortran[7]))
  ksoil_fortran = parse.(Float64, split(fortran[8]))

  ### Run sfexch and ebalsrf

  for i in 1:fsm.Nitr
    sfexch(fsm, meteo)
    ebalsrf(fsm, meteo)
  end

  fortran = readlines("../fortran/temp/test_sfexch.txt")

  global KH_fortran
  global KHa_fortran
  global KHg_fortran
  global KHv_fortran
  global KWg_fortran
  global KWv_fortran
  global Usc_fortran

  KH_fortran = parse.(Float64, split(fortran[1]))
  KHa_fortran = parse.(Float64, split(fortran[2]))
  KHg_fortran = parse.(Float64, split(fortran[3]))
  KHv_fortran = parse.(Float64, split(fortran[4]))
  KWg_fortran = parse.(Float64, split(fortran[5]))
  KWv_fortran = parse.(Float64, split(fortran[6]))
  Usc_fortran = parse.(Float64, split(fortran[7]))

  fortran = readlines("../fortran/temp/test_ebalsrf.txt")

  global Esrf_fortran
  global Eveg_fortran
  global G_fortran
  global H_fortran
  global Hsrf_fortran
  global LE_fortran
  global LEsrf_fortran
  global LWsci_fortran
  global LWveg_fortran
  global Melt_fortran
  global Rnet_fortran
  global Rsrf_fortran

  Esrf_fortran = parse.(Float64, split(fortran[1]))
  Eveg_fortran = parse.(Float64, split(fortran[2]))
  G_fortran = parse.(Float64, split(fortran[3]))
  H_fortran = parse.(Float64, split(fortran[4]))
  Hsrf_fortran = parse.(Float64, split(fortran[5]))
  LE_fortran = parse.(Float64, split(fortran[6]))
  LEsrf_fortran = parse.(Float64, split(fortran[7]))
  LWsci_fortran = parse.(Float64, split(fortran[8]))
  LWveg_fortran = parse.(Float64, split(fortran[9]))
  Melt_fortran = parse.(Float64, split(fortran[10]))
  Rnet_fortran = parse.(Float64, split(fortran[11]))
  Rsrf_fortran = parse.(Float64, split(fortran[12]))

  # Run snow

  snow(fsm, meteo)

  fortran = readlines("../fortran/temp/test_snow.txt")

  global Gsoil_fortran
  global Roff_fortran
  global meltflux_out_fortran
  global Sbsrf_fortran
  global Roff_bare_fortran
  global Roff_snow_fortran
  global fsnow_thres_fortran
  global unload_fortran
  global Tsnow_fortran
  global Sice_fortran
  global Sliq_fortran
  global Ds_fortran

  Gsoil_fortran = parse.(Float64, split(fortran[1]))
  Roff_fortran = parse.(Float64, split(fortran[2]))
  meltflux_out_fortran = parse.(Float64, split(fortran[3]))
  Sbsrf_fortran = parse.(Float64, split(fortran[4]))
  Roff_bare_fortran = parse.(Float64, split(fortran[5]))
  Roff_snow_fortran = parse.(Float64, split(fortran[6]))
  fsnow_thres_fortran = parse.(Float64, split(fortran[7]))
  unload_fortran = parse.(Float64, split(fortran[8]))
  Tsnow_fortran = parse.(Float64, split(fortran[9]))
  Sice_fortran = parse.(Float64, split(fortran[10]))
  Sliq_fortran = parse.(Float64, split(fortran[11]))
  Ds_fortran = parse.(Float64, split(fortran[12]))

  # Run soil

  soil(fsm)

  fortran = readlines("../fortran/temp/test_soil.txt")

  global Tsoil_fortran

  Tsoil_fortran = parse.(Float64, split(fortran[1]))

  if length(output_file) > 0
    println(fout,"$(year) $(month) $(day) $(hour) $(sum(fsm.Ds[:,1,1])) $(fsm.fsnow[1,1]) $(sum(fsm.Sice[:,1,1]+fsm.Sliq[:,1,1])) $(fsm.Tsrf[1,1]) $(fsm.Nsnow[1,1])")
  end

end

if length(output_file) > 0
  close(fout)
end

if check_final_vals

  ### Test radiation

  println("\nTest radiation")

  println(maximum(abs.(alb_fortran - fsm.alb)))
  println(maximum(abs.(asrf_out_fortran - fsm.asrf_out)))
  println(maximum(abs.(Sdirt_fortran - fsm.Sdirt)))
  println(maximum(abs.(Sdift_fortran - fsm.Sdift)))
  println(maximum(abs.(SWveg_fortran - fsm.SWveg)))
  println(maximum(abs.(SWsrf_fortran - fsm.SWsrf)))
  println(maximum(abs.(SWsci_fortran - fsm.SWsci)))
  println(maximum(abs.(LWt_fortran - fsm.LWt)))

  ### Test thermal

  println("\nTest thermal")

  println(maximum(abs.(Ds1_fortran - fsm.Ds1)))
  println(maximum(abs.(gs1_fortran - fsm.gs1)))
  println(maximum(abs.(ks1_fortran - fsm.ks1)))
  println(maximum(abs.(Ts1_fortran - fsm.Ts1)))
  println(maximum(abs.(Tveg0_fortran - fsm.Tveg0)))
  println(maximum(abs.(csoil_fortran - fsm.csoil[:, 1, 1])))
  println(maximum(abs.(ksnow_fortran - fsm.ksnow[:, 1, 1])))
  println(maximum(abs.(ksoil_fortran - fsm.ksoil[:, 1, 1])))

  ### Test sfexch and ebalsrf

  println("\nTest sfexch and ebalsrf")

  println(maximum(abs.(KH_fortran - fsm.KH)))
  println(maximum(abs.(KHa_fortran - fsm.KHa)))
  println(maximum(abs.(KHg_fortran - fsm.KHg)))
  println(maximum(abs.(KHv_fortran - fsm.KHv)))
  println(maximum(abs.(KWg_fortran - fsm.KWg)))
  println(maximum(abs.(KWv_fortran - fsm.KWv)))
  println(maximum(abs.(Usc_fortran - fsm.Usc)))

  println(maximum(abs.(Esrf_fortran - fsm.Esrf)))
  println(maximum(abs.(Eveg_fortran - fsm.Eveg)))
  println(maximum(abs.(G_fortran - fsm.G)))
  println(maximum(abs.(H_fortran - fsm.H)))
  println(maximum(abs.(Hsrf_fortran - fsm.Hsrf)))
  println(maximum(abs.(LE_fortran - fsm.LE)))
  println(maximum(abs.(LEsrf_fortran - fsm.LEsrf)))
  println(maximum(abs.(LWsci_fortran - fsm.LWsci)))
  println(maximum(abs.(LWveg_fortran - fsm.LWveg)))
  println(maximum(abs.(Melt_fortran - fsm.Melt)))
  println(maximum(abs.(Rnet_fortran - fsm.Rnet)))
  println(maximum(abs.(Rsrf_fortran - fsm.Rsrf)))

  # Test snow

  println("\nTest snow")

  println(maximum(abs.(Gsoil_fortran - fsm.Gsoil)))
  println(maximum(abs.(Roff_fortran - fsm.Roff)))
  println(maximum(abs.(meltflux_out_fortran - fsm.meltflux_out)))
  println(maximum(abs.(Sbsrf_fortran - fsm.Sbsrf)))
  println(maximum(abs.(Roff_bare_fortran - fsm.Roff_bare)))
  println(maximum(abs.(Roff_snow_fortran - fsm.Roff_snow)))
  println(maximum(abs.(fsnow_thres_fortran - fsm.fsnow_thres)))
  println(maximum(abs.(unload_fortran - fsm.unload)))
  println(maximum(abs.(Tsnow_fortran - fsm.Tsnow[:, 1, 1])))
  println(maximum(abs.(Sice_fortran - fsm.Sice[:, 1, 1])))
  println(maximum(abs.(Sliq_fortran - fsm.Sliq[:, 1, 1])))
  println(maximum(abs.(Ds_fortran - fsm.Ds[:, 1, 1])))

  # Test soil

  println("\nTest soil")

  println(maximum(abs.(Tsoil_fortran - fsm.Tsoil)))

end
