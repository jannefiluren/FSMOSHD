include("../src/parameters.jl")
include("../src/setup.jl")
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

terrain = readline("../fortran/input/terrain_" * replace(station,"." => "_") * ".txt")
terrain = parse.(Float64,split(terrain,","))

fsky_terr[:,:] .= terrain[1]
slopemu[:,:] .= terrain[2]
xi[:,:] .= terrain[3]
Ld[:,:] .= terrain[4]
lat[:,:] .= terrain[5]
lon[:,:] .= terrain[6]
dem[:,:] .= terrain[7]

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

    println(maximum(abs.(albs_fortran - albs)))
    println(maximum(abs.(Ds_fortran - Ds)))
    println(maximum(abs.(Nsnow_fortran - Nsnow)))
    println(maximum(abs.(Qcan_fortran - Qcan)))
    println(maximum(abs.(Sice_fortran - Sice)))
    println(maximum(abs.(Sliq_fortran - Sliq)))
    println(maximum(abs.(Sveg_fortran - Sveg)))
    println(maximum(abs.(Tcan_fortran - Tcan)))
    println(maximum(abs.(theta_fortran - theta)))
    println(maximum(abs.(Tsnow_fortran - Tsnow)))
    println(maximum(abs.(Tsoil_fortran - Tsoil)))
    println(maximum(abs.(Tsrf_fortran - Tsrf)))
    println(maximum(abs.(fsnow_fortran - fsnow)))
    println(maximum(abs.(Tveg_fortran - Tveg)))
    println(maximum(abs.(snowdepthmin_fortran - snowdepthmin)))
    println(maximum(abs.(snowdepthmax_fortran - snowdepthmax)))
    println(maximum(abs.(snowdepthhist_fortran - snowdepthhist)))
    println(maximum(abs.(swemin_fortran - swemin)))
    println(maximum(abs.(swemax_fortran - swemax)))
    println(maximum(abs.(swehist_fortran - swehist)))
    println(maximum(abs.(fcly_fortran - fcly)))
    println(maximum(abs.(b_fortran - b)))
    println(maximum(abs.(hcap_soil_fortran - hcap_soil)))
    println(maximum(abs.(sathh_fortran - sathh)))
    println(maximum(abs.(Vsat_fortran - Vsat)))
    println(maximum(abs.(Vcrit_fortran - Vcrit)))

end


Qa = similar(Ta)

if length(output_file) > 0
  fout = open(output_file, "w") 
end


for (index,data) in enumerate(readlines(drive_file))

  global fortran

  println("Time step: ", index)

  ### Run drive

  global year, month, day, hour

  year, month, day, hour = drive(data)

  ### Run radiation

  radiation()

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

  thermal()

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

  for i in 1:Nitr
    sfexch()
    ebalsrf()
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

  snow()

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

  soil()

  fortran = readlines("../fortran/temp/test_soil.txt")

  global Tsoil_fortran

  Tsoil_fortran = parse.(Float64, split(fortran[1]))

  if length(output_file) > 0
    println(fout,"$(year) $(month) $(day) $(hour) $(sum(Ds[:,1,1])) $(fsnow[1,1]) $(sum(Sice[:,1,1]+Sliq[:,1,1])) $(Tsrf[1,1]) $(Nsnow[1,1])")
  end

end

if length(output_file) > 0
  close(fout)
end

if check_final_vals

  ### Test radiation

  println("\nTest radiation")

  println(maximum(abs.(alb_fortran - alb)))
  println(maximum(abs.(asrf_out_fortran - asrf_out)))
  println(maximum(abs.(Sdirt_fortran - Sdirt)))
  println(maximum(abs.(Sdift_fortran - Sdift)))
  println(maximum(abs.(SWveg_fortran - SWveg)))
  println(maximum(abs.(SWsrf_fortran - SWsrf)))
  println(maximum(abs.(SWsci_fortran - SWsci)))
  println(maximum(abs.(LWt_fortran - LWt)))

  ### Test thermal

  println("\nTest thermal")

  println(maximum(abs.(Ds1_fortran - Ds1)))
  println(maximum(abs.(gs1_fortran - gs1)))
  println(maximum(abs.(ks1_fortran - ks1)))
  println(maximum(abs.(Ts1_fortran - Ts1)))
  println(maximum(abs.(Tveg0_fortran - Tveg0)))
  println(maximum(abs.(csoil_fortran - csoil[:, 1, 1])))
  println(maximum(abs.(ksnow_fortran - ksnow[:, 1, 1])))
  println(maximum(abs.(ksoil_fortran - ksoil[:, 1, 1])))

  ### Test sfexch and ebalsrf

  println("\nTest sfexch and ebalsrf")

  println(maximum(abs.(KH_fortran - KH)))
  println(maximum(abs.(KHa_fortran - KHa)))
  println(maximum(abs.(KHg_fortran - KHg)))
  println(maximum(abs.(KHv_fortran - KHv)))
  println(maximum(abs.(KWg_fortran - KWg)))
  println(maximum(abs.(KWv_fortran - KWv)))
  println(maximum(abs.(Usc_fortran - Usc)))

  println(maximum(abs.(Esrf_fortran - Esrf)))
  println(maximum(abs.(Eveg_fortran - Eveg)))
  println(maximum(abs.(G_fortran - G)))
  println(maximum(abs.(H_fortran - H)))
  println(maximum(abs.(Hsrf_fortran - Hsrf)))
  println(maximum(abs.(LE_fortran - LE)))
  println(maximum(abs.(LEsrf_fortran - LEsrf)))
  println(maximum(abs.(LWsci_fortran - LWsci)))
  println(maximum(abs.(LWveg_fortran - LWveg)))
  println(maximum(abs.(Melt_fortran - Melt)))
  println(maximum(abs.(Rnet_fortran - Rnet)))
  println(maximum(abs.(Rsrf_fortran - Rsrf)))

  # Test snow

  println("\nTest snow")

  println(maximum(abs.(Gsoil_fortran - Gsoil)))
  println(maximum(abs.(Roff_fortran - Roff)))
  println(maximum(abs.(meltflux_out_fortran - meltflux_out)))
  println(maximum(abs.(Sbsrf_fortran - Sbsrf)))
  println(maximum(abs.(Roff_bare_fortran - Roff_bare)))
  println(maximum(abs.(Roff_snow_fortran - Roff_snow)))
  println(maximum(abs.(fsnow_thres_fortran - fsnow_thres)))
  println(maximum(abs.(unload_fortran - unload)))
  println(maximum(abs.(Tsnow_fortran - Tsnow[:, 1, 1])))
  println(maximum(abs.(Sice_fortran - Sice[:, 1, 1])))
  println(maximum(abs.(Sliq_fortran - Sliq[:, 1, 1])))
  println(maximum(abs.(Ds_fortran - Ds[:, 1, 1])))

  # Test soil

  println("\nTest soil")

  println(maximum(abs.(Tsoil_fortran - Tsoil)))

end
