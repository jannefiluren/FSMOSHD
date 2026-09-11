!-----------------------------------------------------------------------
! Standalone wrapper for the SnowSlide routines
!
! The upstream routines (SNOWSLIDE.F90, SNOW_ABLATION.F90, SWE_FROM_HS.F90)
! read and write their state through modules defined in MODULES.F90. This
! wrapper keeps the upstream files untouched: it receives all state as
! arguments from Julia, copies it into the module variables, calls the
! upstream routine, and copies the results back.
!
! Model parameters (PARAM_SNOWSLIDE, PARAM_SNOWTRAN3D, CONSTANTS) are
! compile-time constants in MODULES.F90 and are therefore not arguments.
!-----------------------------------------------------------------------
subroutine SNOWSLIDE_WRAPPER(Nx_arg, Ny_arg, Nsmax_arg, Ds_min_arg,      &
                             rhos_min_arg, rhos_max_arg, tiled_arg,      &
                             snowdepth0, Sice0, dSWE_slide,              &
                             fsnow_arg, Ds_arg, Sice_arg, Sliq_arg,      &
                             Tsnow_arg, histowet_arg, Nsnow_arg,         &
                             dSWE_tot_slide_arg, index_sorted_dem_arg,   &
                             dem_arg, slope_arg, Shd_arg, forestfrac_arg)

use MODTILE, only: tiled_trans_run

use GRID, only: Nx, Ny, Nsmax, Ds_min

use PARAMETERS, only: rhos_min, rhos_max

use STATE_VARIABLES, only: &
  fsnow, Ds, Sice, Sliq, Tsnow, histowet, Nsnow, &
  dSWE_tot_slide, index_sorted_dem

use LANDUSE, only: dem, slope, Shd, forestfrac

implicit none

integer, intent(in) :: &
  Nx_arg, Ny_arg,          &! Grid dimensions
  Nsmax_arg,               &! Maximum number of snow layers
  tiled_arg                 ! Tiled trans run flag (0 = false, 1 = true)

real, intent(in) :: &
  Ds_min_arg,              &! Minimum possible snow layer thickness (m)
  rhos_min_arg,            &! Minimum snow density (kg/m^3)
  rhos_max_arg              ! Maximum snow density (kg/m^3)

real, intent(inout) :: &
  snowdepth0(Nx_arg,Ny_arg),     &! Snow depth of deposited snow (m)
  Sice0(Nx_arg,Ny_arg),          &! Ice content of deposited snow (kg/m^2)
  dSWE_slide(Nx_arg,Ny_arg)       ! SWE change due to snow slides (kg/m^2)

real, intent(inout) :: &
  fsnow_arg(Nx_arg,Ny_arg),              &! Snow cover fraction
  Ds_arg(Nsmax_arg,Nx_arg,Ny_arg),       &! Snow layer thicknesses (m)
  Sice_arg(Nsmax_arg,Nx_arg,Ny_arg),     &! Ice content of snow layers (kg/m^2)
  Sliq_arg(Nsmax_arg,Nx_arg,Ny_arg),     &! Liquid content of snow layers (kg/m^2)
  Tsnow_arg(Nsmax_arg,Nx_arg,Ny_arg),    &! Snow layer temperatures (K)
  histowet_arg(Nsmax_arg,Nx_arg,Ny_arg), &! Historical variable for past wetting of a layer (0-1)
  dSWE_tot_slide_arg(Nx_arg,Ny_arg)       ! Cumulated SWE change due to snow slides (kg/m^2)

integer, intent(inout) :: &
  Nsnow_arg(Nx_arg,Ny_arg)                ! Number of snow layers

integer, intent(in) :: &
  index_sorted_dem_arg(Nx_arg*Ny_arg,2)   ! Location (i,j) of sorted grid points

real, intent(in) :: &
  dem_arg(Nx_arg,Ny_arg),        &! Terrain elevation (m)
  slope_arg(Nx_arg,Ny_arg),      &! Slope (deg)
  Shd_arg(Nx_arg,Ny_arg),        &! Snow holding depth (m)
  forestfrac_arg(Nx_arg,Ny_arg)   ! Forest fraction

! Set scalar module variables
Nx = Nx_arg
Ny = Ny_arg
Nsmax = Nsmax_arg
Ds_min = Ds_min_arg
rhos_min = rhos_min_arg
rhos_max = rhos_max_arg
tiled_trans_run = (tiled_arg /= 0)

! Copy state into module arrays (allocation on assignment)
fsnow = fsnow_arg
Ds = Ds_arg
Sice = Sice_arg
Sliq = Sliq_arg
Tsnow = Tsnow_arg
histowet = histowet_arg
Nsnow = Nsnow_arg
dSWE_tot_slide = dSWE_tot_slide_arg
index_sorted_dem = index_sorted_dem_arg
dem = dem_arg
slope = slope_arg
Shd = Shd_arg
forestfrac = forestfrac_arg

! Call the upstream routine
call SNOWSLIDE(snowdepth0, Sice0, dSWE_slide)

! Copy modified state back
fsnow_arg = fsnow
Ds_arg = Ds
Sice_arg = Sice
Sliq_arg = Sliq
Tsnow_arg = Tsnow
histowet_arg = histowet
Nsnow_arg = Nsnow
dSWE_tot_slide_arg = dSWE_tot_slide

end subroutine SNOWSLIDE_WRAPPER
