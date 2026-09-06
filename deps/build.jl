using Libdl

# Build SnowSlide routines (upstream files kept verbatim; the wrapper
# provides a standalone entry point for the ccall from Julia)
snowslide_routines = [
    "MODULES.F90",
    "SNOWSLIDE.F90",
    "SNOW_ABLATION.F90",
    "SWE_FROM_HS.F90",
    "SNOWSLIDE_WRAPPER.F90",
]

# Build SnowTran3D routines (upstream files kept verbatim; the wrapper
# provides a standalone entry point for the ccall from Julia)
snowtran3d_routines = [
    "MODULES.F90",
    "SNOWTRAN3D.F90",
    "HS_FROM_SWE.F90",
    "SNOW_ABLATION.F90",
    "SWE_FROM_HS.F90",
    "SNOWTRAN3D_WRAPPER.F90",
]

if Sys.iswindows()
    flags = ["-m$(Sys.WORD_SIZE)", "-shared", "-O3"]
else
    flags = ["-m$(Sys.WORD_SIZE)", "-shared", "-O3", "-fPIC"]
end

try
    # Build SnowSlide library
    run(`gfortran $flags $snowslide_routines -o libsnowslide.$(Libdl.dlext)`)
    println("Built libsnowslide.$(Libdl.dlext) successfully")

    # Build SnowTran3D library
    run(`gfortran $flags $snowtran3d_routines -o libsnowtran3d.$(Libdl.dlext)`)
    println("Built libsnowtran3d.$(Libdl.dlext) successfully")
catch e
    @warn "Failed to build Fortran libraries (libsnowslide, libsnowtran3d) - ensure gfortran is installed and on PATH. The snowslide and snowtran3d routines will be unavailable." exception = e
finally
    # Clean up any .mod files produced during compilation
    for f in readdir()
        if endswith(f, ".mod")
            rm(f)
        end
    end
end
