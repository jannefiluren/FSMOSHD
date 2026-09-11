!-----------------------------------------------------------------------
! Standalone wrapper for the SnowTran3D routines
!
! The upstream routines (SNOWTRAN3D.F90, SNOW_ABLATION.F90, SWE_FROM_HS.F90,
! HS_FROM_SWE.F90) read and write their state through modules defined in
! MODULES.F90. This wrapper keeps the upstream files untouched: it receives
! all state as arguments from Julia, copies it into the module variables,
! calls the upstream routine, and copies the results back.
!
! Model parameters (PARAM_SNOWTRAN3D, CONSTANTS_SNOWTRAN3D, CONSTANTS) are
! compile-time constants in MODULES.F90 and are therefore not arguments.
!-----------------------------------------------------------------------
subroutine SNOWTRAN3D_WRAPPER(Nx_arg, Ny_arg, Nsmax_arg, Ds_min_arg,       &
                              dt_arg, zU_arg, zRH_arg,                     &
                              rhos_min_arg, rhos_max_arg, tiled_arg,       &
                              snowdepth0, Sice0,                           &
                              dSWE_salt, dSWE_susp, dSWE_subl,             &
                              Ua_arg, Udir_arg, Ta_arg, RH_arg,            &
                              veg_shd_arg, z0_snow_arg,                    &
                              fsnow_arg, Ds_arg, Sice_arg, Sliq_arg,       &
                              Tsnow_arg, histowet_arg, Nsnow_arg,          &
                              dSWE_tot_subl_arg, dSWE_tot_salt_arg,        &
                              dSWE_tot_susp_arg,                           &
                              cellsize_arg, dem_arg, forestfrac_arg)

use MODTILE, only: tiled_trans_run

use GRID, only: Nx, Ny, Nsmax, Ds_min

use DRIVING, only: dt, zU, zRH, Ua, Udir, Ta, RH

use PARAMETERS, only: rhos_min, rhos_max

use PARAMMAPS, only: veg_shd, z0_snow

use STATE_VARIABLES, only: &
  fsnow, Ds, Sice, Sliq, Tsnow, histowet, Nsnow, &
  dSWE_tot_subl, dSWE_tot_salt, dSWE_tot_susp

use LANDUSE, only: cellsize, dem, forestfrac

implicit none

integer, intent(in) :: &
  Nx_arg, Ny_arg,          &! Grid dimensions
  Nsmax_arg,               &! Maximum number of snow layers
  tiled_arg                 ! Tiled trans run flag (0 = false, 1 = true)

real, intent(in) :: &
  Ds_min_arg,              &! Minimum possible snow layer thickness (m)
  dt_arg,                  &! Timestep (s)
  zU_arg,                  &! Wind speed measurement height (m)
  zRH_arg,                 &! Relative humidity measurement height (m)
  rhos_min_arg,            &! Minimum snow density (kg/m^3)
  rhos_max_arg              ! Maximum snow density (kg/m^3)

real, intent(inout) :: &
  snowdepth0(Nx_arg,Ny_arg),     &! Snow depth of snowdrift accumulation (m)
  Sice0(Nx_arg,Ny_arg),          &! Ice content of snowdrift accumulation (kg/m^2)
  dSWE_salt(Nx_arg,Ny_arg),      &! SWE change due to saltation (kg/m^2)
  dSWE_susp(Nx_arg,Ny_arg),      &! SWE change due to suspension (kg/m^2)
  dSWE_subl(Nx_arg,Ny_arg)        ! SWE change due to snowdrift sublimation (kg/m^2)

real, intent(in) :: &
  Ua_arg(Nx_arg,Ny_arg),         &! Wind speed (m/s)
  Udir_arg(Nx_arg,Ny_arg),       &! Wind direction (degrees, clockwise from N)
  Ta_arg(Nx_arg,Ny_arg),         &! Air temperature (K)
  RH_arg(Nx_arg,Ny_arg),         &! Relative humidity (%)
  veg_shd_arg(Nx_arg,Ny_arg),    &! Vegetation snow holding depth (m)
  z0_snow_arg(Nx_arg,Ny_arg)      ! Roughness length of snow (m)

real, intent(inout) :: &
  fsnow_arg(Nx_arg,Ny_arg),              &! Snow cover fraction
  Ds_arg(Nsmax_arg,Nx_arg,Ny_arg),       &! Snow layer thicknesses (m)
  Sice_arg(Nsmax_arg,Nx_arg,Ny_arg),     &! Ice content of snow layers (kg/m^2)
  Sliq_arg(Nsmax_arg,Nx_arg,Ny_arg),     &! Liquid content of snow layers (kg/m^2)
  Tsnow_arg(Nsmax_arg,Nx_arg,Ny_arg),    &! Snow layer temperatures (K)
  histowet_arg(Nsmax_arg,Nx_arg,Ny_arg), &! Historical variable for past wetting of a layer (0-1)
  dSWE_tot_subl_arg(Nx_arg,Ny_arg),      &! Cumulated SWE change due to snowdrift sublimation (kg/m^2)
  dSWE_tot_salt_arg(Nx_arg,Ny_arg),      &! Cumulated SWE change due to saltation (kg/m^2)
  dSWE_tot_susp_arg(Nx_arg,Ny_arg)        ! Cumulated SWE change due to suspension (kg/m^2)

integer, intent(inout) :: &
  Nsnow_arg(Nx_arg,Ny_arg)                ! Number of snow layers

real, intent(in) :: &
  cellsize_arg(Nx_arg,Ny_arg),   &! Grid cell size (m)
  dem_arg(Nx_arg,Ny_arg),        &! Terrain elevation (m)
  forestfrac_arg(Nx_arg,Ny_arg)   ! Forest fraction

! Set scalar module variables
Nx = Nx_arg
Ny = Ny_arg
Nsmax = Nsmax_arg
Ds_min = Ds_min_arg
dt = dt_arg
zU = zU_arg
zRH = zRH_arg
rhos_min = rhos_min_arg
rhos_max = rhos_max_arg
tiled_trans_run = (tiled_arg /= 0)

! Copy driving data and parameter maps into module arrays (allocation on assignment)
Ua = Ua_arg
Udir = Udir_arg
Ta = Ta_arg
RH = RH_arg
veg_shd = veg_shd_arg
z0_snow = z0_snow_arg

! Copy state into module arrays
fsnow = fsnow_arg
Ds = Ds_arg
Sice = Sice_arg
Sliq = Sliq_arg
Tsnow = Tsnow_arg
histowet = histowet_arg
Nsnow = Nsnow_arg
dSWE_tot_subl = dSWE_tot_subl_arg
dSWE_tot_salt = dSWE_tot_salt_arg
dSWE_tot_susp = dSWE_tot_susp_arg
cellsize = cellsize_arg
dem = dem_arg
forestfrac = forestfrac_arg

! Call the upstream routine
call SNOWTRAN3D(snowdepth0, Sice0, dSWE_salt, dSWE_susp, dSWE_subl)

! Copy modified state back
fsnow_arg = fsnow
Ds_arg = Ds
Sice_arg = Sice
Sliq_arg = Sliq
Tsnow_arg = Tsnow
histowet_arg = histowet
Nsnow_arg = Nsnow
dSWE_tot_subl_arg = dSWE_tot_subl
dSWE_tot_salt_arg = dSWE_tot_salt
dSWE_tot_susp_arg = dSWE_tot_susp

end subroutine SNOWTRAN3D_WRAPPER
