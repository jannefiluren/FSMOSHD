# Global constant for library path
const LIBSNOWTRAN3D = joinpath(@__DIR__, "..", "deps", "libsnowtran3d")

"""
    snowtran3d!(fsm, met, snowdepth0, Sice0, dSWE_salt, dSWE_susp, dSWE_subl)

Snow transport by wind using Liston's SnowTran3D model.

Implementation of Liston and Sturm (1998) and Liston et al. (2007) SnowTran3D model via Fortran ccall.

# Arguments
- `fsm::FSM`: Model state structure  
- `met::MET`: Meteo variable structure  
- `snowdepth0::Matrix`: Snow depth of deposited snow (m) - modified in-place
- `Sice0::Matrix`: Ice content of deposited snow (kg/m²) - modified in-place  
- `dSWE_salt::Matrix`: SWE change due to saltation (kg/m²) - output
- `dSWE_susp::Matrix`: SWE change due to suspension (kg/m²) - output
- `dSWE_subl::Matrix`: SWE change due to sublimation (kg/m²) - output
"""
function snowtran3d!(
        fsm::FSM{Tf, Ti}, met::MET{Tf, Ti}, snowdepth0::Matrix{Tf}, Sice0::Matrix{Tf},
        dSWE_salt::Matrix{Tf}, dSWE_susp::Matrix{Tf},
        dSWE_subl::Matrix{Tf}
    ) where {Tf <: Real, Ti <: Integer}

    @unpack Nx, Ny, Nsmax, dt, zRH, zU = fsm
    @unpack Ds_min = fsm
    @unpack Ua, Udir, Ta, RH = met
    @unpack vegsnowd_xy, z0_snow = fsm
    @unpack Ds, Nsnow, fsnow, Sice, Sliq, Tsnow, histowet = fsm
    @unpack dSWE_tot_subl, dSWE_tot_salt, dSWE_tot_susp = fsm
    @unpack Ld = fsm
    @unpack rhos_min, rhos_max, rho_snow = fsm


    # Call standalone Fortran routine following exact signature order
    ccall(
        (:snowtran3d_, LIBSNOWTRAN3D),
        Cvoid,
        (
            Ref{Ti}, Ref{Ti}, Ref{Ti},                     # Nsmax, Nx, Ny
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf}, Ptr{Tf}, Ptr{Tf},   # snowdepth0, Sice0, dSWE_salt, dSWE_susp, dSWE_subl
            Ref{Tf},                                       # Ds_min
            Ptr{Tf}, Ptr{Tf}, Ref{Tf}, Ptr{Tf}, Ptr{Tf},   # Ua, Udir, dt, Ta, RH
            Ref{Tf}, Ref{Tf},                              # zRH, zU
            Ptr{Tf}, Ptr{Tf},                              # vegsnowd_xy, z0_snow
            Ptr{Tf}, Ptr{Ti}, Ptr{Tf}, Ptr{Tf}, Ptr{Tf},   # Ds, Nsnow, fsnow, Sice, Sliq
            Ptr{Tf}, Ptr{Tf},                              # Tsnow, histowet
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf},                     # dSWE_tot_subl, dSWE_tot_salt, dSWE_tot_susp
            Ptr{Tf}, Ref{Tf}, Ref{Tf}, Ref{Tf},
        ),           # Ld, rhos_min, rhos_max, rho_snow
        Nsmax, Nx, Ny,
        snowdepth0, Sice0, dSWE_salt, dSWE_susp, dSWE_subl,
        Ds_min,
        Ua, Udir, dt, Ta, RH,
        zRH, zU,
        vegsnowd_xy, z0_snow,
        Ds, Nsnow, fsnow, Sice, Sliq,
        Tsnow, histowet,
        dSWE_tot_subl, dSWE_tot_salt, dSWE_tot_susp,
        Ld, rhos_min, rhos_max, rho_snow
    )

    return nothing
end
