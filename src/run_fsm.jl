function run_fsm_point(station)

  projdir = dirname(dirname(@__FILE__))

  drive_file = joinpath(projdir, "fortran", "input", "input_") * replace(station, "." => "_") * ".txt"
  terrain_file = joinpath(projdir, "fortran", "input", "terrain_") * replace(station, "." => "_") * ".txt"

  drive_data = readdlm(drive_file)

  fsm = FSM{Float64}()
  setup_point!(fsm, terrain_file)

  meteo = MET{Float64}()

  output_data = zeros(size(drive_data, 1), 9)

  for (istep, indata) in enumerate(eachrow(drive_data))

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

    # Run model

    drive!(fsm, meteo)

    radiation(fsm, meteo, t)

    thermal(fsm)

    for i in 1:fsm.Nitr
      sfexch(fsm, meteo)
      ebalsrf(fsm, meteo)
    end

    snow(fsm, meteo)

    soil(fsm)

    # Output data

    output_data[istep, 1] = Dates.value(Year(t))
    output_data[istep, 2] = Dates.value(Month(t))
    output_data[istep, 3] = Dates.value(Day(t))
    output_data[istep, 4] = Dates.value(Hour(t))
    tmpsum = 0.0
    for si in 1:size(fsm.Ds, 1)
      tmpsum += fsm.Ds[si, 1, 1]
    end
    output_data[istep, 5] = tmpsum
    output_data[istep, 6] = fsm.fsnow[1, 1]
    tmpsum = 0.0
    for si in 1:size(fsm.Sice, 1)
      tmpsum += fsm.Sice[si, 1, 1] + fsm.Sliq[si, 1, 1]
    end
    output_data[istep, 7] = tmpsum
    output_data[istep, 8] = fsm.Tsrf[1, 1]
    output_data[istep, 9] = fsm.Nsnow[1, 1]
  end

  return output_data

end

function run_fsm_grid(starttime::DateTime=DateTime(2021,10,01,00,00,00), endtime::DateTime=DateTime(2022,09,30,23,00,00))

  landuse_file_loc = string("/home/haugened/Documents/data/FSM_input/grid/", "BAFU_LUS_2022_250.mat")

  #read landuse from .mat-file
  landuse_file = matopen(landuse_file_loc)
  landuse = read(landuse_file, "landuse")
  close(landuse_file)

  Nx = round(Int, landuse["nrows"])
  Ny = round(Int, landuse["ncols"])

  fsm = FSM{Float64}(Nx=Nx, Ny=Ny)
  setup_grid!(fsm, landuse)

  meteo = MET{Float64}(Nx=Nx, Ny=Ny)

  times = collect(starttime:Hour(1):endtime)

  output_time = zeros(DateTime, size(times, 1))
  output_data = zeros(size(times, 1), fsm.Nx, fsm.Ny, 9)

  for (istep, t) in enumerate(times)

    @show t

    # Run model

    drive_grid!(meteo, fsm, t)

    radiation(fsm, meteo, t)

    thermal(fsm)

    for i in 1:fsm.Nitr
      sfexch(fsm, meteo)
      ebalsrf(fsm, meteo)
    end

    snow(fsm, meteo)

    soil(fsm)

    # Output data

    output_time[istep] = t
    tmpsum = zeros(Float64, fsm.Nx, fsm.Ny)
    for si in 1:size(fsm.Ds, 1)
      tmpsum[:,:] .+= fsm.Ds[si, :, :]
    end
    output_data[istep, :, :, 5] = tmpsum
    output_data[istep, :, :, 6] = fsm.fsnow
    tmpsum = zeros(Float64, fsm.Nx, fsm.Ny)
    for si in 1:size(fsm.Sice, 1)
      tmpsum[:,:] .+= fsm.Sice[si, :, :] + fsm.Sliq[si, :, :]
    end
    output_data[istep, :, :, 7] = tmpsum
    output_data[istep, :, :, 8] = fsm.Tsrf
    output_data[istep, :, :, 9] = fsm.Nsnow
  end

  return output_data

end