using Libdl

function compile()

  mods = ["MODULES.F90" "MODE_WRITE.F90"]

  routines = ["LWRADTOPO.F90" "SWRADTOPO.F90" "CANOPY.F90" "CUMULATE.F90" "DRIVE.F90" "DUMP.F90" "EBALFOR.F90" ;;
              "EBALSRF_SBG.F90" "EBALSRF.F90" "FSM2.F90" "LUDCMP.F90" "OPEN_FILES.F90" "PHYSICS.F90" "QSAT.F90" ;;
              "RADIATION.F90" "SETUP.F90" "SNOW.F90" "SNOWCOVERFRACTION.F90" "SOIL.F90" ;;
              "SFEXCH.F90" "THERMAL.F90" "TRIDIAG.F90" "OUTPUT.F90"]
  
  if Sys.iswindows()
    flags = ["-cpp", "-ffpe-trap=overflow", "-O3"]
  else
    flags = ["-cpp", "-ffpe-trap=overflow", "-O3"]
  end

  run(`gfortran $flags $mods $routines -o FSM2_TXT_64`)

  rm.(filter!(f -> endswith(f, "mod"), readdir()))

  if Sys.iswindows()
    mv("FSM2_TXT_64.exe", "../FSM2_TXT_64.exe", force=true)
  else
    mv("FSM2_TXT_64", "../FSM2_TXT_64", force=true)
  end

end

