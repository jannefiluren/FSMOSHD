function setup_point!(fsm::FSM, terrain_file; state_file="")

  # Adjust parameters depending on tile
  if (fsm.TILE == "forest")
    fsm.hfsn = 0.3    
    fsm.z0sn = 0.005
    fsm.z0sf = 0.2*ones(Nx,Ny)
  end

  if (fsm.TILE == "glacier")
    fsm.alb0 = 0.3*ones(Nx,Ny)
    fsm.z0sf = 0.04*ones(Nx,Ny)
  end

  if (fsm.TILE != "forest") 
    fsm.VAI[:,:]  .= 0
    fsm.hcan[:,:] .= 0
    fsm.fsky[:,:] .= 1
    fsm.trcn[:,:] .= exp(-fsm.kdif.*fsm.VAI[:,:])
    fsm.fveg[:,:] .= 1 .- exp(-fsm.kveg.*fsm.VAI[:,:])
    fsm.fves[:,:] .= 1 .- exp(-fsm.kveg.*fsm.VAI[:,:])
  end

  # Adjust parameters depending on configuration
  if (fsm.DENSTY == 0)
    fsm.rhof = rho0
  end

  # Derived soil parameters
  for j = 1:fsm.Ny
    for i = 1:fsm.Nx
      if (fsm.fcly[i,j] + fsm.fsnd[i,j] > 1)
        fsm.fcly[i,j] = 1 - fsm.fsnd[i,j]
      end
      fsm.b[i,j] = 3.1 + 15.7*fsm.fcly[i,j] - 0.3*fsm.fsnd[i,j]
      fsm.hcap_soil[i,j] = (2.128*fsm.fcly[i,j] + 2.385*fsm.fsnd[i,j])*1e6 / (fsm.fcly[i,j] + fsm.fsnd[i,j])
      fsm.sathh[i,j] = 10^(0.17 - 0.63*fsm.fcly[i,j] - 1.58*fsm.fsnd[i,j])
      fsm.Vsat[i,j] = 0.505 - 0.037*fsm.fcly[i,j] - 0.142*fsm.fsnd[i,j]
      fsm.Vcrit[i,j] = fsm.Vsat[i,j]*(fsm.sathh[i,j]/3.364)^(1/fsm.b[i,j])
      hcon_min = (hcon_clay^fsm.fcly[i,j]) * (hcon_sand^(1 - fsm.fcly[i,j]))
      fsm.hcon_soil[i,j] = (hcon_air^fsm.Vsat[i,j]) * (hcon_min^(1 - fsm.Vsat[i,j]))
    end
  end

  # Initial soil profiles
  fsat = 0.5*ones(fsm.Nsoil)
  Tprof = 285*ones(fsm.Nsoil)
  for k = 1:fsm.Nsoil
    fsm.theta[k,:,:] .= fsat[k]*fsm.Vsat[:,:]
    fsm.Tsoil[k,:,:] .= Tprof[k]
  end
  fsm.Tsrf[:,:] .= fsm.Tsoil[1,:,:]

  # Derived canopy parameters
  fsm.canh[:,:] .= 12500*fsm.VAI[:,:]
  fsm.scap[:,:] .= fsm.cvai*fsm.VAI[:,:]

  # Read terrain parameters from file
  terrain = readline(terrain_file)
  terrain = parse.(Float64,split(terrain,","))

  fsm.fsky_terr[:,:] .= terrain[1]
  fsm.slopemu[:,:] .= terrain[2]
  fsm.xi[:,:] .= terrain[3]
  fsm.Ld[:,:] .= terrain[4]
  fsm.lat[:,:] .= terrain[5]
  fsm.lon[:,:] .= terrain[6]
  fsm.dem[:,:] .= terrain[7]

  fsm.tilefrac[:,:] .= 1   # TODO read this variable from file...

  # Read states from file
  if !isempty(state_file)    # TODO this only work is Nx and Ny == 1
    tmp = readlines(state_file)

    fsm.albs[1,1] = parse.(Float64,split(tmp[1]))[1]
    fsm.Ds[:,1,1] .= parse.(Float64,split(tmp[2]))
    fsm.Nsnow[1,1] = parse.(Float64,split(tmp[3]))[1]
    fsm.Qcan[1,1] = parse.(Float64,split(tmp[4]))[1]
    fsm.Sice[:,1,1] .= parse.(Float64,split(tmp[5]))
    fsm.Sliq[:,1,1] .= parse.(Float64,split(tmp[6]))
    fsm.Sveg[1,1] = parse.(Float64,split(tmp[7]))[1]
    fsm.Tcan[1,1] = parse.(Float64,split(tmp[8]))[1]
    fsm.theta[:,1,1] .= parse.(Float64,split(tmp[9]))
    fsm.Tsnow[:,1,1] .= parse.(Float64,split(tmp[10]))
    fsm.Tsoil[:,1,1] .= parse.(Float64,split(tmp[11]))
    fsm.Tsrf[1,1] = parse.(Float64,split(tmp[12]))[1]
    fsm.fsnow[1,1] = parse.(Float64,split(tmp[13]))[1]
    fsm.Tveg[1,1] = parse.(Float64,split(tmp[14]))[1]
    fsm.snowdepthmin[1,1] = parse.(Float64,split(tmp[15]))[1]
    fsm.snowdepthmax[1,1] = parse.(Float64,split(tmp[16]))[1]
    fsm.snowdepthhist[:,1,1] .= parse.(Float64,split(tmp[17]))
    fsm.swemin[1,1] = parse.(Float64,split(tmp[18]))[1]
    fsm.swemax[1,1] = parse.(Float64,split(tmp[19]))[1]
    fsm.swehist[:,1,1] .= parse.(Float64,split(tmp[20]))
  end

end

function setup_grid!(fsm::FSM, landuse; state_file="")

  # Adjust parameters depending on tile
  if (fsm.TILE == "forest")
    fsm.hfsn = 0.3    
    fsm.z0sn = 0.005
    fsm.z0sf = 0.2*ones(Nx,Ny)
  end

  if (fsm.TILE == "glacier")
    fsm.alb0 = 0.3*ones(Nx,Ny)
    fsm.z0sf = 0.04*ones(Nx,Ny)
  end

  if (fsm.TILE != "forest") 
    fsm.VAI[:,:]  .= 0
    fsm.hcan[:,:] .= 0
    fsm.fsky[:,:] .= 1
    fsm.trcn[:,:] .= exp.(-fsm.kdif.*fsm.VAI[:,:])
    fsm.fveg[:,:] .= 1 .- exp.(-fsm.kveg.*fsm.VAI[:,:])
    fsm.fves[:,:] .= 1 .- exp.(-fsm.kveg.*fsm.VAI[:,:])
  end

  # Adjust parameters depending on configuration
  if (fsm.DENSTY == 0)
    fsm.rhof = rho0
  end

  # Derived soil parameters
  for j = 1:fsm.Ny
    for i = 1:fsm.Nx
      if (fsm.fcly[i,j] + fsm.fsnd[i,j] > 1)
        fsm.fcly[i,j] = 1 - fsm.fsnd[i,j]
      end
      fsm.b[i,j] = 3.1 + 15.7*fsm.fcly[i,j] - 0.3*fsm.fsnd[i,j]
      fsm.hcap_soil[i,j] = (2.128*fsm.fcly[i,j] + 2.385*fsm.fsnd[i,j])*1e6 / (fsm.fcly[i,j] + fsm.fsnd[i,j])
      fsm.sathh[i,j] = 10^(0.17 - 0.63*fsm.fcly[i,j] - 1.58*fsm.fsnd[i,j])
      fsm.Vsat[i,j] = 0.505 - 0.037*fsm.fcly[i,j] - 0.142*fsm.fsnd[i,j]
      fsm.Vcrit[i,j] = fsm.Vsat[i,j]*(fsm.sathh[i,j]/3.364)^(1/fsm.b[i,j])
      hcon_min = (hcon_clay^fsm.fcly[i,j]) * (hcon_sand^(1 - fsm.fcly[i,j]))
      fsm.hcon_soil[i,j] = (hcon_air^fsm.Vsat[i,j]) * (hcon_min^(1 - fsm.Vsat[i,j]))
    end
  end

  # Initial soil profiles
  fsat = 0.5*ones(fsm.Nsoil)
  Tprof = 285*ones(fsm.Nsoil)
  for k = 1:fsm.Nsoil
    fsm.theta[k,:,:] .= fsat[k]*fsm.Vsat[:,:]
    fsm.Tsoil[k,:,:] .= Tprof[k]
  end
  fsm.Tsrf[:,:] .= fsm.Tsoil[1,:,:]

  # Derived canopy parameters
  fsm.canh[:,:] .= 12500*fsm.VAI[:,:]
  fsm.scap[:,:] .= fsm.cvai*fsm.VAI[:,:]

  #skyview-factor
  fsm.fsky_terr = landuse["res_skyvf"]
  fsm.fsky_terr[isnan.(fsm.fsky_terr)] .= 1.0 #TODO lus.sbg_skyvf does not exist in .mat-file

  fsm.slopemu = sqrt.(landuse["dhdxdy"] ./ 2)
  fsm.xi = (sqrt(2) .* landuse["sd"]) ./ fsm.slopemu

  fsm.Ld[:,:] .= 1.0

  latlon!(fsm, landuse)

  fsm.dem .= 1234.0 #landuse["dem"] #TODO: what to do with NaNs in DEM?

  fsm.tilefrac .= 1   # TODO read this variable from file...

  # Read states from file
  if !isempty(state_file)    # TODO this only work is Nx and Ny == 1
    tmp = readlines(state_file)

    fsm.albs[1,1] = parse.(Float64,split(tmp[1]))[1]
    fsm.Ds[:,1,1] .= parse.(Float64,split(tmp[2]))
    fsm.Nsnow[1,1] = parse.(Float64,split(tmp[3]))[1]
    fsm.Qcan[1,1] = parse.(Float64,split(tmp[4]))[1]
    fsm.Sice[:,1,1] .= parse.(Float64,split(tmp[5]))
    fsm.Sliq[:,1,1] .= parse.(Float64,split(tmp[6]))
    fsm.Sveg[1,1] = parse.(Float64,split(tmp[7]))[1]
    fsm.Tcan[1,1] = parse.(Float64,split(tmp[8]))[1]
    fsm.theta[:,1,1] .= parse.(Float64,split(tmp[9]))
    fsm.Tsnow[:,1,1] .= parse.(Float64,split(tmp[10]))
    fsm.Tsoil[:,1,1] .= parse.(Float64,split(tmp[11]))
    fsm.Tsrf[1,1] = parse.(Float64,split(tmp[12]))[1]
    fsm.fsnow[1,1] = parse.(Float64,split(tmp[13]))[1]
    fsm.Tveg[1,1] = parse.(Float64,split(tmp[14]))[1]
    fsm.snowdepthmin[1,1] = parse.(Float64,split(tmp[15]))[1]
    fsm.snowdepthmax[1,1] = parse.(Float64,split(tmp[16]))[1]
    fsm.snowdepthhist[:,1,1] .= parse.(Float64,split(tmp[17]))
    fsm.swemin[1,1] = parse.(Float64,split(tmp[18]))[1]
    fsm.swemax[1,1] = parse.(Float64,split(tmp[19]))[1]
    fsm.swehist[:,1,1] .= parse.(Float64,split(tmp[20]))
  end

end

function latlon!(fsm, landuse) #TODO
  fsm.lat .= 46.83
  fsm.lon .= 9.8092

end