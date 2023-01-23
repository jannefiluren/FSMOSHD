function run_fsm_point(station)

  projdir = dirname(dirname(@__FILE__))

  drive_file = joinpath(projdir, "fortran", "input", "input_") * replace(station, "." => "_") * ".txt"
  terrain_file = joinpath(projdir, "fortran", "input", "terrain_") * replace(station,"." => "_") * ".txt"

  drive_data = readdlm(drive_file)

  fsm = FSM{Float64}()
  setup_point!(fsm, terrain_file)

  meteo = MET{Float64}()

  output_data = zeros(size(drive_data, 1), 9)

  for (istep, indata) in enumerate(eachrow(drive_data))

    # Forcing data

    year = indata[1]
    month = indata[2]
    day = indata[3]
    hour = indata[4]
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

    # Run model

    drive!(fsm, meteo)

    radiation(fsm, year, month, day, hour, meteo)

    thermal(fsm)
 
    for i in 1:fsm.Nitr
      sfexch(fsm, meteo)
      ebalsrf(fsm, meteo)
    end

    snow(fsm, meteo)
    
    soil(fsm)

    # Output data

    output_data[istep, 1] = year
    output_data[istep, 2] = month
    output_data[istep, 3] = day
    output_data[istep, 4] = hour
    tmpsum = 0.0
    for si in 1:size(fsm.Ds, 1)
        tmpsum += fsm.Ds[si,1,1]
    end
    output_data[istep, 5] = tmpsum
    output_data[istep, 6] = fsm.fsnow[1,1]
    tmpsum = 0.0
    for si in 1:size(fsm.Sice, 1)
        tmpsum += fsm.Sice[si,1,1]+fsm.Sliq[si,1,1]
    end
    output_data[istep, 7] = tmpsum
    output_data[istep, 8] = fsm.Tsrf[1,1]
    output_data[istep, 9] = fsm.Nsnow[1,1]
  end
  
  return output_data

end

function run_fsm_grid()

  #projdir = dirname(dirname(@__FILE__))

  landuse_file_loc = string("/home/haugened/Documents/data/FSM_input/grid/", "BAFU_LUS_0250_2023a.mat")

  #drive_data = readdlm(drive_file)

  #read landuse from .mat-file
  landuse_file = matopen(landuse_file_loc)
  landuse = read(landuse_file, "landuse")
  close(landuse_file)

  Nx = round(Int, landuse["nrows"])
  Ny = round(Int, landuse["ncols"])

  fsm = FSM{Float64}(Nx = Nx, Ny = Ny)
  setup_grid!(fsm, landuse)

  Sdir = zeros(fsm.Nx, fsm.Ny)
  Sdif = zeros(fsm.Nx, fsm.Ny)
  LW = zeros(fsm.Nx, fsm.Ny)
  Sf = zeros(fsm.Nx, fsm.Ny)
  Rf = zeros(fsm.Nx, fsm.Ny)
  Ta = zeros(fsm.Nx, fsm.Ny)
  RH = zeros(fsm.Nx, fsm.Ny)
  Ua = zeros(fsm.Nx, fsm.Ny)
  Ps = zeros(fsm.Nx, fsm.Ny)
  Sf24h = zeros(fsm.Nx, fsm.Ny)
  Tc = zeros(fsm.Nx, fsm.Ny)
  es = zeros(fsm.Nx, fsm.Ny)
  Qa = zeros(fsm.Nx, fsm.Ny)
  Tv = zeros(fsm.Nx, fsm.Ny)

  output_data = zeros(size(drive_data, 1), 9)

  for (istep, indata) in enumerate(eachrow(drive_data))

    # Forcing data

    year = indata[1]
    month = indata[2]
    day = indata[3]
    hour = indata[4]
    Sdir[:, :] .= indata[5]
    Sdif[:, :] .= indata[6]
    LW[:, :] .= indata[7]
    Sf[:, :] .= indata[8]
    Rf[:, :] .= indata[9]
    Ta[:, :] .= indata[10]
    RH[:, :] .= indata[11]
    Ua[:, :] .= indata[12]
    Ps[:, :] .= indata[13]
    Sf24h[:, :] .= indata[14]

    # Run model

    drive!(fsm, Tc, es, Qa, Ua, Sf, Rf, Ta, RH, Ps)

    radiation(fsm, year, month, day, hour, LW, Sdif, Sdir, Sf, Sf24h, Ta, Tv)

    thermal(fsm)
 
    for i in 1:fsm.Nitr
      sfexch(fsm, Ta, Ps, Qa, Ua)
      ebalsrf(fsm, LW, Ps, Qa, Ta)
    end

    snow(fsm, Rf, Sf, Ta, Ua)
    
    soil(fsm)

    # Output data

    output_data[istep, 1] = year
    output_data[istep, 2] = month
    output_data[istep, 3] = day
    output_data[istep, 4] = hour
    tmpsum = 0.0
    for si in 1:size(fsm.Ds, 1)
        tmpsum += fsm.Ds[si,1,1]
    end
    output_data[istep, 5] = tmpsum
    output_data[istep, 6] = fsm.fsnow[1,1]
    tmpsum = 0.0
    for si in 1:size(fsm.Sice, 1)
        tmpsum += fsm.Sice[si,1,1]+fsm.Sliq[si,1,1]
    end
    output_data[istep, 7] = tmpsum
    output_data[istep, 8] = fsm.Tsrf[1,1]
    output_data[istep, 9] = fsm.Nsnow[1,1]
  end
  
  return output_data

end