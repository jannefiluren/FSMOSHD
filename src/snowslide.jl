# Global constant for library path
const LIBSNOWSLIDE = joinpath(@__DIR__, "..", "deps", "libsnowslide")

"""
    snowslide!(fsm, snowdepth0, Sice0, dSWE_slide)

Lateral redistribution of snow through gravity using the SnowSlide model.

Implementation of Bernhardt and Schulz (2010) SnowSlide model via Fortran ccall.
Reference: Quéno et al. (2024)

# Arguments
- `fsm::FSM`: Model state structure  
- `snowdepth0::Matrix`: Snow depth of deposited snow (m) - modified in-place
- `Sice0::Matrix`: Ice content of deposited snow (kg/m²) - modified in-place  
- `dSWE_slide::Matrix`: SWE change due to snow slides (kg/m²) - output
"""
function snowslide!(
        fsm::FSM{Tf, Ti}, snowdepth0::Matrix{Tf},
        Sice0::Matrix{Tf}, dSWE_slide::Matrix{Tf}
    ) where {Tf <: Real, Ti <: Integer}

    @unpack Nx, Ny, Nsmax = fsm
    @unpack fsnow, Ds, dSWE_tot_slide, index_sorted_dem = fsm
    @unpack dem, slope, Shd = fsm
    @unpack dyn_ratio, rho_deposit, slope_min, Shd_min, Ds_min, rho_snow = fsm
    @unpack Sice, Sliq, Nsnow, Tsnow, histowet, rhos_min, rhos_max = fsm

    # Constants
    Tm = Tf(273.15)     # TODO get from constants function...

    # Call standalone Fortran routine
    ccall(
        (:snowslide_, LIBSNOWSLIDE),
        Cvoid,
        (
            Ref{Ti}, Ref{Ti}, Ref{Ti},                    # Nx, Ny, Nsmax
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf},                    # snowdepth0, Sice0, dSWE_slide
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf}, Ptr{Ti},           # fsnow, Ds, dSWE_tot_slide, index_sorted_dem
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf},                    # dem, slope, Shd
            Ref{Tf}, Ref{Tf}, Ref{Tf}, Ref{Tf},           # dyn_ratio, rho_deposit, slope_min, Shd_min
            Ref{Tf}, Ref{Tf}, Ref{Tf},                    # Ds_min, Tm, rho_snow
            Ptr{Tf}, Ptr{Tf}, Ptr{Ti}, Ptr{Tf}, Ptr{Tf},  # Sice, Sliq, Nsnow, Tsnow, histowet
            Ref{Tf}, Ref{Tf},
        ),                            # rhos_min, rhos_max
        Nx, Ny, Nsmax,
        snowdepth0, Sice0, dSWE_slide,
        fsnow, Ds, dSWE_tot_slide, index_sorted_dem,
        dem, slope, Shd,
        dyn_ratio, rho_deposit, slope_min, Shd_min,
        Ds_min, Tm, rho_snow,
        Sice, Sliq, Nsnow, Tsnow, histowet,
        rhos_min, rhos_max
    )

    return nothing
end
