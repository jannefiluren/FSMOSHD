# Global constant for library path
const LIBSNOWSLIDE = joinpath(@__DIR__, "..", "..", "deps", "libsnowslide")

"""
    snowslide!(fsm, w, snowdepth0, Sice0, dSWE_slide)

Lateral redistribution of snow through gravity using the SnowSlide model.

Implementation of Bernhardt and Schulz (2010) SnowSlide model via Fortran ccall.
The upstream Fortran routines are kept verbatim in deps/ and are entered through
the standalone wrapper in deps/SNOWSLIDE_WRAPPER.F90. Model parameters
(dyn_ratio, trig_ratio, rho_deposit, slope_min, Shd_min, rho_snow) are
hardcoded constants in deps/MODULES.F90.
Reference: Quéno et al. (2024)

# Arguments
- `fsm::FSM`: Model state structure
- `w::SnowTransport`: Transport workspace (static arrays, accumulators, constants)
- `snowdepth0::Matrix`: Snow depth of deposited snow (m) - modified in-place
- `Sice0::Matrix`: Ice content of deposited snow (kg/m²) - modified in-place
- `dSWE_slide::Matrix`: SWE change due to snow slides (kg/m²) - output
"""
function snowslide!(
        fsm::FSM{Tf}, w::SnowTransport{Tf}, snowdepth0::Matrix{Tf},
        Sice0::Matrix{Tf}, dSWE_slide::Matrix{Tf}
    ) where {Tf <: Real}

    (; Nx, Ny, Nsmax) = fsm.grid
    (; Ds_min) = fsm.params
    (; rhos_min, rhos_max, tiled_trans_run) = w
    (; fsnow, Ds, Sice, Sliq, Tsnow, histowet, Nsnow) = fsm.state
    (; dSWE_tot_slide, index_sorted_dem, slope, Shd, forestfrac) = w
    (; dem) = fsm.surface

    # The Fortran wrapper is compiled for 4-byte integers. Nsnow is in-out (staged as an Int32
    # copy and read back); index_sorted_dem is read-only (a throwaway Int32 copy). Integer
    # scalars convert through their `Ref{Int32}` slots automatically.
    Nsnow32 = Int32.(Nsnow)
    index_sorted_dem32 = Int32.(index_sorted_dem)

    # Call the standalone Fortran wrapper
    ccall(
        (:snowslide_wrapper_, LIBSNOWSLIDE),
        Cvoid,
        (
            Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Tf},  # Nx, Ny, Nsmax, Ds_min
            Ref{Tf}, Ref{Tf}, Ref{Int32},                 # rhos_min, rhos_max, tiled_trans_run
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf},                    # snowdepth0, Sice0, dSWE_slide
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf}, Ptr{Tf},           # fsnow, Ds, Sice, Sliq
            Ptr{Tf}, Ptr{Tf}, Ptr{Int32},                 # Tsnow, histowet, Nsnow
            Ptr{Tf}, Ptr{Int32},                          # dSWE_tot_slide, index_sorted_dem
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf}, Ptr{Tf},           # dem, slope, Shd, forestfrac
        ),
        Nx, Ny, Nsmax, Ds_min,
        rhos_min, rhos_max, Int32(tiled_trans_run),
        snowdepth0, Sice0, dSWE_slide,
        fsnow, Ds, Sice, Sliq,
        Tsnow, histowet, Nsnow32,
        dSWE_tot_slide, index_sorted_dem32,
        dem, slope, Shd, forestfrac
    )

    Nsnow .= Nsnow32

    return nothing
end
