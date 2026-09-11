!-----------------------------------------------------------------------
! Snow transport by wind - Liston's scheme adapted to FSM 2.0
! References: Liston and Sturm (1998), Liston et al. (2007)
! Author: Glen E. Liston
! Adapted to FSM 2.0 by Louis Quéno
! Reference: Quéno et al. (2024)
! ATTENTION: y is the W->E axis, while x is the S->N
! South of (i,j): (i-1,j)
! North of (i,j): (i+1,j)
! West of (i,j): (i,j-1)
! East of (i,j): (i,j+1)
!-----------------------------------------------------------------------
subroutine SNOWTRAN3D(snowdepth0,Sice0,dSWE_salt,dSWE_susp,dSWE_subl)

use CONSTANTS, only: &
  pi,                   &! pi
  rho_wat,              &! Density of water (kg/m^3)
  Tm                     ! Melting point (K)

use DRIVING, only: &  
  Ua,                   &! Wind speed (m/s)
  Udir                   ! Wind direction (degrees, clockwise from N)

use GRID, only: &
  Nx,Ny                  ! Grid dimensions

use PARAMMAPS, only: &
  veg_shd                ! Vegetation snow holding depth (m)

use STATE_VARIABLES, only: &
  Ds,                   &! Snow layer thicknesses (m)
  Nsnow,                &! Number of snow layers
  fsnow,                &! Snow cover fraction 
  Sice,                 &! Ice content of snow layers (kg/m^2)
  Sliq,                 &! Liquid content of snow layers (kg/m^2)
  Tsnow                  ! Snow layer temperatures (K)

use LANDUSE, only: &
  cellsize,          &! Grid cell size (m)
  dem                 ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  flag_variable_Utau_t, &! Flag for variable threshold friction velocitiy (.TRUE.) or constant (.FALSE.)
  Utau_t_const,         &! Constant threshold friction velocity (m/s) used if flag_variable_Utau_t = .FALSE.
  rho_snow               ! Constant snow density (kg/m^3)

use CONSTANTS_SNOWTRAN3D, only: &
  wind_min               ! Minimum wind speed to compute snow transport (m/s)

implicit none

real, intent(inout) :: &
  snowdepth0(Nx,Ny),    &! Snow depth of snowdrift accumulation, averaged over the grid cell (m)
  Sice0(Nx,Ny),         &! Ice content of snowdrift accumulation(kg/m^2)
  dSWE_salt(Nx,Ny),     &! SWE change due to saltation (kg/m^2)
  dSWE_susp(Nx,Ny),     &! SWE change due to suspension (kg/m^2)
  dSWE_subl(Nx,Ny)       ! SWE change due to snowdrift sublimation (kg/m^2)

integer :: &
  i,j                    ! Point counters

integer :: &
  index_ue(Nx,2*Ny+1),  &! Wind index array E
  index_uw(Nx,2*Ny+1),  &! Wind index array W
  index_vn(Ny,2*Nx+1),  &! Wind index array N
  index_vs(Ny,2*Nx+1)    ! Wind index array S

logical :: &
  flag_blowing_snow      ! Blowing snow flag

real :: &
  delta_WE,             &! Grid cell size in WE direction (m)
  delta_SN,             &! Grid cell size in SN direction (m)
  Ua_max                 ! Maximum wind speed on the domain (m/s)

real :: &
  Qsalt(Nx,Ny),         &! Saltation flux (kg/m/s)
  Qsalt_u(Nx,Ny),       &! x component of saltation flux (kg/m/s)
  Qsalt_v(Nx,Ny),       &! y component of saltation flux (kg/m/s)
  conc_salt(Nx,Ny),     &! Saltation-layer reference-level mass concentration (kg/m^3)
  Qsusp(Nx,Ny),         &! Suspension flux (kg/m/s)
  Qsusp_u(Nx,Ny),       &! x component of suspension flux (kg/m/s)
  Qsusp_v(Nx,Ny),       &! y component of suspension flux (kg/m/s)
  Qsubl(Nx,Ny),         &! Sublimation flux (kg/m^2/s)
  Ds_soft(Nx,Ny),       &! Soft snow thickness (m)
  snowthickness(Nx,Ny), &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  uwind(Nx,Ny),         &! x component of wind speed (m/s)
  vwind(Nx,Ny),         &! y component of wind speed (m/s)
  Utau(Nx,Ny),          &! Friction velocity (m/s)
  Utau_t(Nx,Ny),        &! Threshold friction velocity (m/s)
  h_star(Nx,Ny),        &! Height of the saltation layer (m)
  veg_z0(Nx,Ny),        &! Vegetation roughness length (m)
  z_0(Nx,Ny)             ! Surface roughness length (m)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialization of delta_WE and delta_SN
delta_WE = cellsize(1,1)
delta_SN = cellsize(1,1)

! Initialize Ds_soft and snowdepth0
Ds_soft(:,:) = 0.0
snowdepth0(:,:) = 0.0

! Initialize all fluxes at 0
Qsalt(:,:) = 0.0
Qsalt_u(:,:) = 0.0
Qsalt_v(:,:) = 0.0
dSWE_salt(:,:) = 0.0
conc_salt(:,:) = 0.0
Qsusp(:,:) = 0.0
Qsusp_u(:,:) = 0.0
Qsusp_v(:,:) = 0.0
dSWE_susp(:,:) = 0.0
Qsubl(:,:) = 0.0
dSWE_subl(:,:) = 0.0

! Initialization of uwind and vwind
uwind = Ua * cos(-Udir * pi/180.0 - pi/2.0)
vwind = Ua * sin(-Udir * pi/180.0 - pi/2.0)

! Initialization of maximum wind speed on the domain (m/s)
! Initialization of snowthickness
Ua_max = 0.0
do j = 1, Ny
  do i = 1, Nx

    Ua_max = max(Ua_max,Ua(i,j))
    snowthickness(i,j) = sum(Ds(:,i,j))

  end do
end do

! Initialization of veg_z0 (copied from preprocess_code.f)
! Define the roughness lengths for each of the vegetation types.
! Note that this is pretty arbitrary, because these values are
! really only used when there is no blowing snow, and thus have
! no impact on the simulation except to provide a non-zero value
! for any such parts of the domain.
veg_z0(:,:) = 0.25 * veg_shd(:,:)

! Update the thicknesses of the soft snow layer.
call compute_soft_snow(Ds_soft)

! Update the threshold friction velocity.
if (.NOT. flag_variable_Utau_t) then
  Utau_t(:,:) = Utau_t_const
else
  call surface_snow(Utau_t)
end if

! Set the blowing snow flag to false until it is clear that we will
! have blowing snow.
flag_blowing_snow = .FALSE.

! If the wind speed is lower that some threshold, then don't
! need to to any of the snow transport computations.
if (Ua_max >= wind_min) then

  ! Get the wind direction indexing arrays for this particular
  ! wind event (time step).
  call getdirection(index_ue,index_uw,index_vn,index_vs,uwind,vwind)

  ! Solve for Utau and z_0 if snow is saltating, else solve assuming
  ! z_0 is known from snow depth and/or veg type, and solve for Utau.
  call solveUtau(Utau,z_0,h_star,snowthickness,veg_z0,flag_blowing_snow,Utau_t,Ds_soft)

  ! Update Ds_soft by removing too dense layers compared to U_tau
  call update_soft_snow(Utau,Ds_soft)

  ! If the blowing snow flag indicates wind transported snow
  ! somewhere within the domain (flag_blowing_snow = .TRUE.), run the saltation
  ! and suspension models.
  if (flag_blowing_snow) then

    ! Solve for the saltation flux.
    call saltation(delta_WE,Utau,Utau_t,snowthickness,delta_SN, &
                   index_ue,index_uw,index_vn,index_vs,uwind,vwind,Ds_soft, &
                   Qsalt,Qsalt_u,Qsalt_v)

    ! Solve for the suspension flux.
    call suspension(Utau,z_0,h_star,Utau_t,uwind,vwind, &
                    Qsalt,conc_salt,Qsusp,Qsusp_u,Qsusp_v,Qsubl)

  end if

end if

! Compute the new snow depth due to accumulation from precipitation,
! saltation, and suspension, and the mass loss due to sublimation.
call accum(snowthickness, &
           delta_WE,delta_SN, &
           index_ue,index_uw,index_vn,index_vs, &
           flag_blowing_snow, &
           uwind,vwind,Ds_soft, &
           Qsalt_u,Qsalt_v, &
           dSWE_salt, &
           Qsusp_u,Qsusp_v, &
           dSWE_susp, &
           Qsubl, &
           dSWE_subl, &
           snowdepth0,Sice0)

end subroutine SNOWTRAN3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_soft_snow(Ds_soft)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers
  histowet,          &! Historical variable for past wetting of a layer (0-1)
  Sliq                ! Liquid content of snow layers (kg/m^2)

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

implicit none

real, intent(inout) :: &
  Ds_soft(Nx,Ny)      ! Soft snow thickness (m)

integer :: &
  i,j,               &! Point counters
  k                   ! Layer counter

! New determination of soft snow based on layer wetting history and temperature
! Now histowet is used by default
do j = 1, Ny
  do i = 1, Nx

    k = 1
    do while (k <= Nsnow(i,j) .AND. Sliq(k,i,j) < epsilon(Sliq) .AND. histowet(k,i,j) < 0.5)
      Ds_soft(i,j) = Ds_soft(i,j) + Ds(k,i,j)
      k = k + 1
    end do

  end do
end do

end subroutine compute_soft_snow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_soft_snow(Utau,Ds_soft)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq                ! Liquid content of snow layers (kg/m^2)

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

implicit none

real, intent(in) :: &
  Utau(Nx,Ny)         ! Friction velocity (m/s)

real, intent(inout) :: &
  Ds_soft(Nx,Ny)      ! Soft snow thickness (m)

integer :: &
  i,j,               &! Point counters
  k                   ! Layer counter

logical :: &
  hard_layer_reached  ! A hard layer has been reached

real :: &
  Ds_soft_old,       &! Soft snow thickness before update (m)
  rho_layer,         &! Snow layer density
  Utau_t_layer        ! Threshold friction velocity of each layer (m/s)

do j = 1, Ny
  do i = 1, Nx

    Ds_soft_old = Ds_soft(i,j)
    Ds_soft(i,j) = 0.0
    hard_layer_reached = .FALSE.
    k = 1

    do while ((k <= Nsnow(i,j)) .AND. (Ds_soft_old - Ds_soft(i,j) > epsilon(Ds_soft)) .AND. (.NOT. hard_layer_reached))

      ! Calculate snow layer density
      rho_layer = (Sice(k,i,j) + Sliq(k,i,j)) / Ds(k,i,j) / fsnow(i,j)

      ! Calculate the snow threshold friction velocity.
      if (rho_layer <= 300.0) then
        Utau_t_layer = 0.10 * exp(0.003 * rho_layer)
      else
        Utau_t_layer = 0.005 * exp(0.013 * rho_layer)
      end if

      if (Utau(i,j) - Utau_t_layer > epsilon(Utau)) then
        Ds_soft(i,j) = Ds_soft(i,j) + Ds(k,i,j)
      else
        hard_layer_reached = .TRUE.
      end if

      k = k + 1

    end do

  end do
end do

end subroutine update_soft_snow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine accum(snowthickness, &
                 delta_WE,delta_SN, &
                 index_ue,index_uw,index_vn,index_vs, &
                 flag_blowing_snow, &
                 uwind,vwind,Ds_soft, &
                 Qsalt_u,Qsalt_v, &
                 dSWE_salt, &
                 Qsusp_u,Qsusp_v, &
                 dSWE_susp, &
                 Qsubl, &
                 dSWE_subl, &
                 snowdepth0,Sice0)

use MODTILE, only: &
  tiled_trans_run              ! Tiled trans run flag

use CONSTANTS, only: &
  rho_wat                 ! Density of water (kg/m^3)

use DRIVING, only: &
  dt                      ! Timestep (s)

use GRID, only: &
  Nx,Ny                   ! Grid dimensions

use PARAMMAPS, only: &
  veg_shd                 ! Vegetation snow holding depth (m)

use STATE_VARIABLES, only: &
  Ds,                    &! Snow layer thicknesses (m)
  fsnow,                 &! Snow cover fraction 
  Sice,                  &! Ice content of snow layers (kg/m^2)
  Sliq,                  &! Liquid content of snow layers (kg/m^2)
  dSWE_tot_subl,         &! Cumulated SWE change due to snowdrift sublimation (kg/m^2)
  dSWE_tot_salt,         &! Cumulated SWE change due to saltation (kg/m^2)
  dSWE_tot_susp           ! Cumulated SWE change due to suspension (kg/m^2)

  use LANDUSE, only: &
  dem,                  &! Terrain elevation (m)
  forestfrac             ! Forest fraction

use PARAM_SNOWTRAN3D, only: &
  rho_snow                ! Constant snow density (kg/m^3)

implicit none

real, intent(in) :: &
  uwind(Nx,Ny),          &! WE component of wind speed (m/s)
  vwind(Nx,Ny)            ! SN component of wind speed (m/s)

real, intent(in) :: &
  delta_WE,              &! Grid cell size in WE direction (m)
  delta_SN                ! Grid cell size in SN direction (m)

logical, intent(in) :: &
  flag_blowing_snow       ! Blowing snow flag

real, intent(inout) :: &
  Qsalt_u(Nx,Ny),        &! x component of saltation flux (kg/m/s)
  Qsalt_v(Nx,Ny),        &! y component of saltation flux (kg/m/s)
  dSWE_salt(Nx,Ny),      &! SWE change due to saltation (kg/m^2)
  Qsusp_u(Nx,Ny),        &! x component of suspension flux (kg/m/s)
  Qsusp_v(Nx,Ny),        &! y component of suspension flux (kg/m/s)
  dSWE_susp(Nx,Ny),      &! SWE change due to suspension (kg/m^2)
  Qsubl(Nx,Ny),          &! Sublimation flux (kg/m^2/s)
  dSWE_subl(Nx,Ny),      &! SWE change due to snowdrift sublimation (kg/m^2)
  snowthickness(Nx,Ny),  &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  Ds_soft(Nx,Ny),        &! Soft snow thickness (m)
  snowdepth0(Nx,Ny),     &! Snow depth of snowdrift accumulation, averaged over the grid cell (m)
  Sice0(Nx,Ny)            ! Ice content of transported snow (kg/m^2)

integer, intent(in) :: &
  index_ue(Nx,2*Ny+1),   &! Wind index array E
  index_uw(Nx,2*Ny+1),   &! Wind index array W
  index_vn(Ny,2*Nx+1),   &! Wind index array N
  index_vs(Ny,2*Nx+1)     ! Wind index array S

integer :: &
  i,j                     ! Point counters

real :: &
  dSWE_subl_loss,        &! SWE loss due to sublimation (kg/m^2)
  dh_subl_loss,          &! Snow depth loss due to sublimation (m)
  Ds_hard,               &! Hard snow thickness (m)
  swe_loc,               &! Total snowpack SWE (kg/m^2)
  snowdmin                ! Minimum snow depth to allow transport (m)

if (flag_blowing_snow) then

  ! SALTATION
  call getnewdepth(delta_WE,delta_SN,Qsalt_u, &
                   Qsalt_v,dSWE_salt,index_ue,index_uw, &
                   index_vn,index_vs, &
                   snowthickness,Ds_soft,snowdepth0,Sice0)

  ! SUSPENSION
  call getnewdepth(delta_WE,delta_SN,Qsusp_u, &
                   Qsusp_v,dSWE_susp,index_ue,index_uw, &
                   index_vn,index_vs, &
                   snowthickness,Ds_soft,snowdepth0,Sice0)

  ! SUBLIMATION
  do j = 1, Ny
    do i = 1, Nx

      ! Make adjustments for the case where there is no snow available
      ! on the ground (or captured within the vegetation) to be
      ! eroded.
      Ds_hard = snowthickness(i,j) - Ds_soft(i,j)
      snowdmin = max(veg_shd(i,j),Ds_hard)
      swe_loc = sum(Sice(:,i,j) + Sliq(:,i,j))

      ! Convert Qsubl to sublimated snow depth dh_subl_loss.
      dSWE_subl(i,j) = Qsubl(i,j) * dt * fsnow(i,j)
      dSWE_subl_loss = - dSWE_subl(i,j)
      ! The following min should actually never be necessary, because handled by the snowdmin check later.
      ! Kept for clarity and safety.
      dSWE_subl_loss = min(dSWE_subl_loss,swe_loc)
      dSWE_subl(i,j) = - dSWE_subl_loss
      if (dSWE_subl_loss > epsilon(dSWE_subl_loss)) then
        call HS_FROM_SWE(dSWE_subl_loss,dh_subl_loss,i,j)
      else
        dSWE_subl_loss = 0.0
        dh_subl_loss = 0.0
        dSWE_subl(i,j) = 0.0
      end if

      if (snowthickness(i,j) > snowdmin .and. fsnow(i,j) > epsilon(fsnow)) then
        if (snowthickness(i,j) - dh_subl_loss / fsnow(i,j) <= snowdmin) then
          dh_subl_loss = (snowthickness(i,j) - snowdmin) * fsnow(i,j)
          if (dh_subl_loss > epsilon(dh_subl_loss)) then
            call SWE_FROM_HS(dh_subl_loss,dSWE_subl_loss,i,j)
            dSWE_subl(i,j) = - dSWE_subl_loss
          else
            dh_subl_loss = 0.0
            dSWE_subl_loss = 0.0
            dSWE_subl(i,j) = 0.0
          end if
        end if
      else
        dSWE_subl(i,j) = 0.0
        dSWE_subl_loss = 0.0
        dh_subl_loss = 0.0
      end if

      ! The SWE loss due to wind only occurs on the open part of the pixel.
      ! If it is a tiled run or open run, the transferred SWE is then weighted differently.
      ! If tiled run, in the open tile, SWE loss stays the same since it will be weighted later when combining the tiles.
      if (.not. tiled_trans_run) then
        dSWE_subl(i,j) = dSWE_subl(i,j) * (1.0 - forestfrac(i,j))
        dSWE_subl_loss = dSWE_subl_loss * (1.0 - forestfrac(i,j))
        ! Now, we may need to adjust dh_subl_loss.
        if (dSWE_subl_loss > epsilon(dSWE_subl_loss)) then
          call HS_FROM_SWE(dSWE_subl_loss,dh_subl_loss,i,j)
        else
          dSWE_subl_loss = 0.0
          dSWE_subl(i,j) = 0.0
          dh_subl_loss = 0.0
        end if
      end if

      if (dSWE_subl_loss > epsilon(dSWE_subl_loss) .and. dh_subl_loss > epsilon(dh_subl_loss)) then
        call SNOW_ABLATION(dh_subl_loss,dSWE_subl_loss,i,j)
      end if

      ! Update the snow layer thicknesses
      snowthickness(i,j) = sum(Ds(:,i,j))

    end do
  end do

end if

! Update cumulated sublimation and transport
do j = 1, Ny
  do i = 1, Nx

    ! Fill summing arrays of the sublimation and transport quantities.
    dSWE_tot_subl(i,j) = dSWE_tot_subl(i,j) + dSWE_subl(i,j)
    dSWE_tot_salt(i,j) = dSWE_tot_salt(i,j) + dSWE_salt(i,j) 
    dSWE_tot_susp(i,j) = dSWE_tot_susp(i,j) + dSWE_susp(i,j)

  end do
end do

end subroutine accum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine suspension(Utau,z_0,h_star,Utau_t,uwind,vwind, &
                      Qsalt,conc_salt,Qsusp,Qsusp_u,Qsusp_v,Qsubl)

use CONSTANTS, only: &
  vkman               ! Von Karman constant

use DRIVING, only: &
  Ta,                &! Air temperature (K)
  RH                  ! Relative humidity (%)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

use CONSTANTS_SNOWTRAN3D, only: &
  Up_const,          &! Constant coefficient for calculation of U_p
  dz_susp,           &! dz in suspension layer
  ztop_susp,         &! Height of the top of suspension layer (m)
  fall_vel,          &! Particle-settling velocity, s in L&S 1998 eq. 12 (m/s)
  Ur_const            ! Constant coefficient beta in L&S 1998, eq. 13

implicit none

real, intent(in) :: &
  Qsalt(Nx,Ny),      &! Saltation flux (kg/m/s)
  uwind(Nx,Ny),      &! x component of wind speed (m/s)
  vwind(Nx,Ny),      &! y component of wind speed (m/s)
  Utau_t(Nx,Ny),     &! Threshold friction velocity (m/s)
  Utau(Nx,Ny),       &! Friction velocity (m/s)
  h_star(Nx,Ny),     &! Height of the saltation layer (m)
  z_0(Nx,Ny)          ! Surface roughness length (m)

real, intent(inout) :: &
  conc_salt(Nx,Ny),  &! Saltation-layer reference-level mass concentration (kg/m^3)
  Qsusp(Nx,Ny),      &! Suspension flux (kg/m/s)
  Qsusp_u(Nx,Ny),    &! x component of suspension flux (kg/m/s)
  Qsusp_v(Nx,Ny),    &! y component of suspension flux (kg/m/s)
  Qsubl(Nx,Ny)        ! Sublimation flux (kg/m^2/s)

integer :: &
  i,j,               &! Point counters
  nzsteps,           &! Number of vertical steps
  iz                  ! Vertical steps counter

real :: &
  conc,              &! Concentration of the suspended snow at height z (kg/m^3)
  prd,               &! Auxiliary variable for calculations
  U_p,               &! Horizontal particle velocity within the saltation layer (m/s)
  U_r,               &! Wind speed at height h_star (m/s)
  phistar_Cr,        &! Left member of L&S 1998 eq. 13
  V_susp,            &! Sublimation loss rate coefficient in the suspension layer (s^-1)
  V_salt,            &! Sublimation loss rate coefficient in the saltation layer (s^-1)
  z                   ! Height (m)

! Compute the mass concentration of suspended snow according to
! Kind (1992).

do j = 1, Ny
  do i = 1, Nx

    if (Qsalt(i,j) > epsilon(Qsalt)) then
      ! LQ: the present subroutine (suspension) is only called in the case there is saltation (flag_blowing_snow = .TRUE.)
      ! in which case (see third case of solveUtau): h* = (h_const/C_z) z0 = 13.3 z0
      ! So no need to check h*, the following log is always positive
      U_r = Utau(i,j)/vkman * log(h_star(i,j)/z_0(i,j))
      phistar_Cr = Utau(i,j)/U_r * Ur_const
      prd = phistar_Cr * Utau(i,j) / fall_vel
      U_p = Up_const * Utau_t(i,j)

      ! Compute the concentration in the saltation layer (kg/m^3).
      conc_salt(i,j) = Qsalt(i,j) / (h_star(i,j) * U_p)

      nzsteps = int((ztop_susp - h_star(i,j)) / dz_susp)

      Qsusp(i,j) = 0.0
      Qsubl(i,j) = 0.0

      do iz = 1, nzsteps
        z = h_star(i,j) + 0.5 * dz_susp + real(iz - 1) * dz_susp

        ! Compute the concentration of the suspended snow at height z.
        conc = conc_salt(i,j) * ((prd + 1.0) * &
               (z/h_star(i,j))**((-fall_vel)/(vkman*Utau(i,j))) - &
               prd)
        conc = max(conc,0.0)

        ! Only do the integration if the concentration is non-zero.
        if (conc > epsilon(conc)) then

          ! Compute the sublimation due to suspension.
          call getsublim(z,RH(i,j),Ta(i,j),Utau(i,j), &
                         z_0(i,j),V_susp,V_salt,Utau_t(i,j),.TRUE.)

          ! Perform the quadrature (summation), without the constants.
          if (z == z_0(i,j)) z = 1.2 * z_0(i,j)
          Qsusp(i,j) = Qsusp(i,j) + conc * log(z/z_0(i,j)) * dz_susp
          Qsubl(i,j) = Qsubl(i,j) + conc * V_susp * dz_susp

        end if

      end do

      ! Finish the quadratures.
      ! Include the constants for Qsusp.
      Qsusp(i,j) = Utau(i,j) / vkman * Qsusp(i,j)

      ! Include the sublimation contribution due to saltation.
      z = h_star(i,j) / 2.0
      call getsublim(z,RH(i,j),Ta(i,j),Utau(i,j), &
                     z_0(i,j),V_susp,V_salt,Utau_t(i,j),.FALSE.)

      Qsubl(i,j) = Qsubl(i,j) + V_salt * conc_salt(i,j) * h_star(i,j)

    else
      conc_salt(i,j) = 0.0
      Qsusp(i,j) = 0.0
      Qsubl(i,j) = 0.0
    end if

  end do
end do

! Separate the east-west and the north-south suspended transport
! components; the vector sum should equal Qsusp.
do j = 1, Ny
  do i = 1, Nx

    Qsusp_u(i,j) = Qsusp(i,j) * abs(uwind(i,j)) / &
                   sqrt(uwind(i,j)**2 + vwind(i,j)**2)
    Qsusp_v(i,j) = Qsusp(i,j) * abs(vwind(i,j)) / &
                   sqrt(uwind(i,j)**2 + vwind(i,j)**2)

  end do
end do

end subroutine suspension

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine saltation(delta_WE,Utau,Utau_t,snowthickness,delta_SN, &
                     index_ue,index_uw,index_vn,index_vs,uwind,vwind,Ds_soft, &
                     Qsalt,Qsalt_u,Qsalt_v)

use CONSTANTS, only: &
  grav,                    &! Acceleration due to gravity (m/s^2)
  rho_air                   ! Density of air (kg/m^3)

use GRID, only: &
  Nx,Ny                     ! Grid dimensions

use PARAMMAPS, only: &
  veg_shd                   ! Vegetation snow holding capacity (m)

use LANDUSE, only: &
  dem,                     &! Terrain elevation (m)
  forestfrac                ! Forest fraction

use PARAM_SNOWTRAN3D, only: &
  flag_boundary_condition, &! Boundary condition flag
  blowby                    ! Fraction of the saltation flux transferred downwind

use CONSTANTS_SNOWTRAN3D, only: &
  fetch,                   &! Equilibrium fetch distance (m)
  xmu                       ! Scaling constant for non-equilibrium saltation transport

implicit none

real, intent(in) :: &
  delta_WE,                &! Grid cell size in WE direction (m)
  delta_SN                  ! Grid cell size in SN direction (m)

real, intent(in) :: &
  snowthickness(Nx,Ny),    &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  uwind(Nx,Ny),            &! x component of wind speed (m/s)
  vwind(Nx,Ny),            &! y component of wind speed (m/s)
  Utau_t(Nx,Ny),           &! Threshold friction velocity (m/s)
  Utau(Nx,Ny),             &! Friction velocity (m/s)
  Ds_soft(Nx,Ny)            ! Soft snow thickness (m)

real, intent(inout) :: &
  Qsalt(Nx,Ny),            &! Saltation flux (kg/m/s)
  Qsalt_u(Nx,Ny),          &! x component of saltation flux (kg/m/s)
  Qsalt_v(Nx,Ny)            ! y component of saltation flux (kg/m/s)

integer, intent(in) :: &
  index_ue(Nx,2*Ny+1),     &! Wind index array E
  index_uw(Nx,2*Ny+1),     &! Wind index array W
  index_vn(Ny,2*Nx+1),     &! Wind index array N
  index_vs(Ny,2*Nx+1)       ! Wind index array S

integer :: &
  i,j,k,                   &! Point counters
  istart,iend,             &! Point counters boundaries
  jstart,jend               ! Point counters boundaries

real :: &
  Qsalt_max(Nx,Ny),        &! Maximum possible saltation flux (kg/m/s)
  Qsalt_maxu(Nx,Ny),       &! x component of Qsalt_max (kg/m/s)
  Qsalt_maxv(Nx,Ny)         ! y component of Qsalt_max (kg/m/s)

real :: &
  dUtau,                   &! Utau difference (m/s)
  scale_EW,                &! Scaling coefficient for Eqn. 9 in L&S 1998
  scale_NS                  ! Scaling coefficient for Eqn. 9 in L&S 1998

! Compute the maximum possible saltation flux, assuming that
! an abundance of snow is available at the surface.
do j = 1, Ny
  do i = 1, Nx

    ! For a given wind speed, find Qsalt_max.
    Qsalt_max(i,j) = 0.68 * rho_air / grav * &
                     Utau_t(i,j) / Utau(i,j) * (Utau(i,j)**2 - Utau_t(i,j)**2)
    Qsalt_max(i,j) = max(Qsalt_max(i,j),0.0)

    ! Now weight the max saltation flux for the u and v wind
    ! components, where the vector sum should equal Qsalt_max.
    Qsalt_maxu(i,j) = Qsalt_max(i,j) * abs(uwind(i,j)) / &
                      sqrt(uwind(i,j)**2 + vwind(i,j)**2)
    Qsalt_maxv(i,j) = Qsalt_max(i,j) * abs(vwind(i,j)) / &
                      sqrt(uwind(i,j)**2 + vwind(i,j)**2)

  end do
end do

! Define an upwind boundary condition.  If flag_boundary_condition = .TRUE. then it is
! assumed that the inflow saltation flux has reached steady state.
! If flag_boundary_condition = .FALSE. then the saltation flux is assumed to be zero.
! The boundary condition is implemented by initializing the arrays
! to Qsalt_max, and since upwind boundaries are not called in
! the Qsalt computation, they stay in effect for the future 
! accumulation/erosion computation.
if (.NOT. flag_boundary_condition) then
  do j = 1, Ny
    do i = 1, Nx

      ! Zero incoming flux at the boundaries.
      Qsalt_u(i,j) = 0.0
      Qsalt_v(i,j) = 0.0

    end do
  end do
else
  do j = 1, Ny
    do i = 1, Nx

      ! Steady-state (maximum) incoming flux at the boundaries.
      Qsalt_u(i,j) = Qsalt_maxu(i,j)
      Qsalt_v(i,j) = Qsalt_maxv(i,j)

    end do
  end do
end if

! Define the scaling coefficients for Eqn. 9 in L&S 1998. Don't
! let them be greater than 1.0 or you will make more snow than
! there was before.
scale_EW =  xmu * delta_WE / fetch
scale_EW = min(1.0,scale_EW)
scale_NS =  xmu * delta_SN / fetch
scale_NS = min(1.0,scale_NS)

! Consider WESTERLY winds.
do i = 1, Nx
  do k = 1, index_uw(i,1)
    jstart = index_uw(i,k*2)+1
    jend = index_uw(i,k*2+1)
    do j = jstart, jend
      dUtau = Utau(i,j) - Utau(i,j-1)
      if (dUtau >= epsilon(dUtau)) then
        Qsalt_u(i,j) = Qsalt_u(i,j-1) + scale_EW * &
                       (Qsalt_maxu(i,j) - Qsalt_u(i,j-1))
      else
!        Qsalt_u(i,j) = min(Qsalt_u(i,j-1),Qsalt_maxu(i,j))
        if (Qsalt_u(i,j-1) < Qsalt_maxu(i,j)) then
          Qsalt_u(i,j) = Qsalt_u(i,j-1)
        else
          Qsalt_u(i,j) = max(blowby*Qsalt_u(i,j-1),Qsalt_maxu(i,j))
        end if
      end if
      ! No flux from fully-forested cells.
      ! For safety, because very high snow holding depth should anyway prevent it.
      if (forestfrac(i,j) > 0.9) then
        Qsalt_u(i,j) = 0.0
      end if
      ! Downwind forest cover reduces the flow.
      if (j+1 <= Ny) then
        if (forestfrac(i,j+1) > 0.9) then
          Qsalt_u(i,j) = 0.0
        else
          Qsalt_u(i,j) = Qsalt_u(i,j) * (1.0 - forestfrac(i,j+1))
        end if
      end if
    end do
  end do
end do

! Consider EASTERLY winds.
do i = 1, Nx
  do k = 1, index_ue(i,1)
    jend = index_ue(i,k*2)
    jstart = index_ue(i,k*2+1)-1
    do j = jstart, jend,-1
      dUtau = Utau(i,j) - Utau(i,j+1)
      if (dUtau >= epsilon(dUtau)) then
        Qsalt_u(i,j) = Qsalt_u(i,j+1) + scale_EW * &
                       (Qsalt_maxu(i,j) - Qsalt_u(i,j+1))
      else
!        Qsalt_u(i,j) = min(Qsalt_u(i,j+1),Qsalt_maxu(i,j))
        if (Qsalt_u(i,j+1) < Qsalt_maxu(i,j)) then
          Qsalt_u(i,j) = Qsalt_u(i,j+1)
        else
          Qsalt_u(i,j) = max(blowby*Qsalt_u(i,j+1),Qsalt_maxu(i,j))
        end if
      end if
      ! No flux from fully-forested cells.
      ! For safety, because very high snow holding depth should anyway prevent it.
      if (forestfrac(i,j) > 0.9) then
        Qsalt_u(i,j) = 0.0
      end if
      ! Downwind forest cover reduces the flow.
      if (j-1 >= 1) then
        if (forestfrac(i,j-1) > 0.9) then
          Qsalt_u(i,j) = 0.0
        else        
          Qsalt_u(i,j) = Qsalt_u(i,j) * (1.0 - forestfrac(i,j-1))
        end if
      end if
    end do
  end do
end do

! Consider SOUTHERLY winds.
do j = 1, Ny
  do k = 1, index_vs(j,1)
    istart = index_vs(j,k*2)+1
    iend = index_vs(j,k*2+1)
    do i = istart, iend
      dUtau = Utau(i,j) - Utau(i-1,j)
      if (dUtau >= epsilon(dUtau)) then
        Qsalt_v(i,j) = Qsalt_v(i-1,j) + scale_NS * &
                       (Qsalt_maxv(i,j) - Qsalt_v(i-1,j))
      else
!        Qsalt_v(i,j) = min(Qsalt_v(i-1,j),Qsalt_maxv(i,j))
        if (Qsalt_v(i-1,j) < Qsalt_maxv(i,j)) then
          Qsalt_v(i,j) = Qsalt_v(i-1,j)
        else
          Qsalt_v(i,j) = max(blowby*Qsalt_v(i-1,j),Qsalt_maxv(i,j))
        end if
      end if
      ! No flux from fully-forested cells.
      ! For safety, because very high snow holding depth should anyway prevent it.
      if (forestfrac(i,j) > 0.9) then
        Qsalt_v(i,j) = 0.0
      end if
      ! Downwind forest cover reduces the flow.
      if (i+1 <= Nx) then
        if (forestfrac(i+1,j) > 0.9) then
          Qsalt_v(i,j) = 0.0
        else        
          Qsalt_v(i,j) = Qsalt_v(i,j) * (1.0 - forestfrac(i+1,j))
        end if
      end if
    end do
  end do
end do

! Consider NORTHERLY winds.
do j = 1, Ny
  do k = 1, index_vn(j,1)
    iend = index_vn(j,k*2)
    istart = index_vn(j,k*2+1)-1
    do i = istart, iend,-1
      dUtau = Utau(i,j) - Utau(i+1,j)
      if (dUtau >= epsilon(dUtau)) then
        Qsalt_v(i,j) = Qsalt_v(i+1,j) + scale_NS * &
                       (Qsalt_maxv(i,j) - Qsalt_v(i+1,j))
      else
!        Qsalt_v(i,j) = min(Qsalt_v(i+1,j),Qsalt_maxv(i,j))
        if (Qsalt_v(i+1,j) < Qsalt_maxv(i,j)) then
          Qsalt_v(i,j) = Qsalt_v(i+1,j)
        else
          Qsalt_v(i,j) = max(blowby*Qsalt_v(i+1,j),Qsalt_maxv(i,j))
        end if
      end if
      ! No flux from fully-forested cells.
      ! For safety, because very high snow holding depth should anyway prevent it.
      if (forestfrac(i,j) > 0.9) then
        Qsalt_v(i,j) = 0.0
      end if
      ! Downwind forest cover reduces the flow.
      if (i-1 >= 1) then
        if (forestfrac(i-1,j) > 0.9) then
          Qsalt_v(i,j) = 0.0
        else        
          Qsalt_v(i,j) = Qsalt_v(i,j) * (1.0 - forestfrac(i-1,j))
        end if
      end if
    end do
  end do
end do

! Combine the u and v components to yield the total saltation flux
! at each grid cell.
do j = 1, Ny
  do i = 1, Nx

    Qsalt(i,j) = Qsalt_u(i,j) + Qsalt_v(i,j)

  end do
end do

! Adjust Qsalt to account for the availablity of snow for transport;
! taking into consideration whether there is snow on the ground,
! the holding depth of the vegetation, etc..
do j = 1, Ny
  do i = 1, Nx

    if (snowthickness(i,j) <= veg_shd(i,j)) then
      Qsalt(i,j) = 0.0
      Qsalt_u(i,j) = 0.0
      Qsalt_v(i,j) = 0.0
    end if
    if (Ds_soft(i,j) <= epsilon(Ds_soft)) then
      Qsalt(i,j) = 0.0
      Qsalt_u(i,j) = 0.0
      Qsalt_v(i,j) = 0.0
    end if

  end do
end do

end subroutine saltation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine solveUtau(Utau,z_0,h_star,snowthickness,veg_z0,flag_blowing_snow,Utau_t,Ds_soft)

use CONSTANTS, only: &
  grav,                 &! Acceleration due to gravity (m/s^2)
  vkman                  ! Von Karman constant

use DRIVING, only: &
  Ua,                   &! Wind speed (m/s)
  zU                     ! Wind speed measurement height (m)

use GRID, only: &
  Nx,Ny                  ! Grid dimensions

use PARAMMAPS, only: &
  veg_shd,              &! Vegetation snow holding depth (m)
  z0_snow                ! Roughness length of snow (m)

use LANDUSE, only: &
  dem                    ! Terrain elevation (m)

use CONSTANTS_SNOWTRAN3D, only: &
  C_z,                  &! Coefficient 0.12 in Liston and Sturm (1998) eq. 5 p. 500
  h_const                ! Coefficient 1.6 in Liston and Sturm (1998) eq. 14 p. 501

implicit none

real, intent(in) :: &
  snowthickness(Nx,Ny), &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  Utau_t(Nx,Ny),        &! Threshold friction velocity (m/s)
  veg_z0(Nx,Ny),        &! Vegetation roughness length (m)
  Ds_soft(Nx,Ny)         ! Soft snow thickness (m)

logical, intent(inout) :: &
  flag_blowing_snow           ! Blowing snow flag

real, intent(inout) :: &
  Utau(Nx,Ny)            ! Friction velocity (m/s)

real, intent(out) :: &  
  h_star(Nx,Ny),        &! Height of the saltation layer (m)
  z_0(Nx,Ny)             ! Surface roughness length (m)

integer :: &
  i,j                    ! Point counters

logical :: &
  threshold_flag         ! Flag if u* is above threshold

real :: &
  threshold,            &! Threshold for u* (m/s)
  guess,                &! Initial guess for Utau
  sfrac,                &! Depth-fraction of vegetation covered by snow
  Utautmp,              &! Temporary Utau variable (m/s)
  windtmp,              &! Temporary wind speed variable (m/s)
  wind_max,             &! Maximum wind speed (m/s)
  z_0_tmp                ! Temporary z_0 variable (m)

! We are entering the subroutine with the blowing snow flag set to .FALSE.
! Then, if snow is found to blow in any
! domain grid cell, we set the flag to (flag_blowing_snow = .TRUE.).

! Build the Utau array.
guess = 0.1
do j = 1, Ny
  do i = 1, Nx

    ! Determine whether snow is saltating (this influences how Utau
    ! and z_0 are computed).
    if (snowthickness(i,j) <= veg_shd(i,j)) then

      ! Saltation will not occur.
      sfrac = snowthickness(i,j) / max(veg_shd(i,j),veg_z0(i,j)) ! Eq. 3, LS 1998
      z_0(i,j) = sfrac * z0_snow(i,j) + (1.0 - sfrac) * veg_z0(i,j) ! Eq. 4, LS 1998
      z_0_tmp = min(0.25*zU,z_0(i,j))
      Utau(i,j) = Ua(i,j) * vkman / log(zU/z_0_tmp)
      h_star(i,j) = z_0(i,j) * h_const / C_z

    else if (Ds_soft(i,j) <= epsilon(Ds_soft)) then

      ! Saltation will not occur.
      z_0(i,j) = z0_snow(i,j)
      Utau(i,j) = Ua(i,j) * vkman / log(zU/z_0(i,j))
      h_star(i,j) = z_0(i,j)

    else

      ! Saltation may occur.  Test for that possibility by assuming that
      ! saltation is present, solving for Utau and z_0, and comparing
      ! whether Utau exceeds Utau_t.  If it does not, set z_0 to that
      ! of snow and recompute Utau.

      ! To help insure that the iteration converges, set the minimum
      ! wind speed to be 1.0 m/s, and the maximum wind speed to be
      ! 30 m/s at 10-m height.
      windtmp = max(1.0,Ua(i,j))
      wind_max = 30.0 * log(zU/z0_snow(i,j))/log(10.0/z0_snow(i,j))
      windtmp = min(windtmp,wind_max) 

      ! For u* over 0.6, use the relation z0 = 0.00734 u* - 0.0022,
      ! instead of Equation (5) in Liston and Sturm (1998).  Note that
      ! for windspeeds greater than about 35 m/s this will have to be
      ! modified for the solution algorithm to converge (because the
      ! roughness length will start to be higher than the obs height!).
      threshold = 0.6/vkman * log(zU/0.0022)
      if (windtmp <= threshold) then
        threshold_flag = .FALSE.
      else
        threshold_flag = .TRUE.
      end if

      call solve1(Utautmp,guess,windtmp, &
                  threshold_flag)

      if (Utautmp > Utau_t(i,j)) then

        ! We have saltation.
        Utau(i,j) = Utautmp
        z_0(i,j) = C_z * Utau(i,j)**2 / (2.0 * grav)
        h_star(i,j) = h_const * Utau(i,j)**2 / (2.0 * grav)
        ! Set the blowing snow flag to .TRUE.
        flag_blowing_snow = .TRUE.

      else

        ! We do not have saltation, but the vegetation is covered by snow.
        ! Because we have determined that we do not have saltation, make
        ! sure Utau does not exceed Utau_t.
        z_0(i,j) = z0_snow(i,j)
        Utau(i,j) = Ua(i,j) * vkman / log(zU/z_0(i,j))
        Utau(i,j) = min(Utau(i,j),Utau_t(i,j))
        h_star(i,j) = z_0(i,j) * h_const / C_z

      end if
    end if

  end do
end do

end subroutine solveUtau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine solve1(xnew,guess,windtmp,threshold_flag)

use CONSTANTS, only: &
  grav,              &! Acceleration due to gravity (m/s^2)
  vkman               ! Von Karman constant

use DRIVING, only: &
  zU                  ! Wind speed measurement height (m)

use CONSTANTS_SNOWTRAN3D, only: &
  C_z                 ! Coefficient 0.12 in Liston and Sturm (1998) eq. 5 p. 500

implicit none

real, intent(in) :: &
  guess,             &! Initial guess for Utau (m/s)
  windtmp             ! Wind speed (m/s)

logical, intent(in) :: &
  threshold_flag      ! Flag if u* is above threshold

real, intent(inout) :: &  
  xnew                ! Solution for Utau (m/s)

integer :: &
  i,                 &! Point counter
  maxiter             ! Number of iterations

real :: &
  tol,               &! Convergence threshold
  old,               &! Temporary xnew
  fprime,            &! Temporary function 1
  funct               ! Temporary function 2

tol = 1.0e-3
maxiter = 20
old = guess

if (threshold_flag) then

  old = 0.6
  do i = 1, maxiter
    fprime = - 1.0 + windtmp * vkman * 0.00734 / (0.00734 * old - 0.0022) * &
            (log(zU) - log(0.00734 * old - 0.0022))**(-2)
    funct = - old + windtmp * vkman * &
           (log(zU) - log(0.00734 * old - 0.0022))**(-1)
    xnew = old - funct/fprime
    if (abs(xnew - old) < tol) return
    old = xnew
  end do

else

  do i = 1, maxiter
    fprime = - 1.0 + 2.0 / old * windtmp * vkman * &
            (log(zU) - log(C_z/(2.0*grav)) - 2.0*log(old))**(-2)
    funct = - old + windtmp * vkman * &
           (log(zU) - log(C_z/(2.0*grav)) - 2.0*log(old))**(-1)
    xnew = old - funct/fprime
    if (abs(xnew - old) < tol) return
    old = xnew
  end do

end if

write(*,*) 'max iteration exceeded when solving for Utau, Utau=',old

end subroutine solve1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getsublim(z,RH,Ta,Utau,z_0,V_susp,V_salt,Utau_t,flag_susp)

use CONSTANTS, only: &
  hcon_air,          &! Thermal conductivity of air (W/m/K)
  Ls,                &! Latent heat of sublimation (J/kg)
  pi,                &! pi
  Rair,              &! Gas constant for air (J/K/kg)
  Runi,              &! Universal gas constant (J/kmol/K)
  rho_ice,           &! Density of ice (kg/m^3)
  Tm,                &! Melting point (K)
  vkman,             &! Von Karman constant
  xM,                &! Molecular weight of water (kg/kmol)
  visc_air            ! Kinematic viscosity of air (m^2/s)

use DRIVING, only: &
  zRH                 ! Relative humidity measurement height (m)

implicit none

real, intent(in) :: &
  z,                 &! Height (m)
  RH,                &! Relative humidity (%)
  Ta,                &! Air temperature (K)
  Utau_t,            &! Threshold friction velocity (m/s)
  Utau,              &! Friction velocity (m/s)
  z_0                 ! Surface roughness length (m)

logical, intent(in) :: &
  flag_susp           ! .TRUE.: sublimation due to suspension, .FALSE.: sublimation due to saltation

real, intent(out) :: &
  V_susp,            &! Sublimation loss rate coefficient in the suspension layer (s^-1)
  V_salt              ! Sublimation loss rate coefficient in the saltation layer (s^-1)

real :: &
  D,                 &! Diffusivity of water vapor in the atmosphere (m^2/s)
  rho_sat,           &! Saturation density of water vapour (kg/m^3)
  rh_offset,         &! Relative humidity offset (fraction 0->1)
  sigma,             &! Atmospheric undersaturation of water vapour with respect to ice at height z
  alpha,             &! Coefficient alpha defined in L&S (1998), eq. A-5 p. 515 
  rbar_r,            &! Mean radius of snow particles at height z (m)
  xmbar,             &! Mean particle mass at height z (m)
  rbar,              &! Radius of a snow particle of mean particle mass xmbar at height z (m)
  u_z,               &! Wind speed at height z (m)
  x_r,               &! Fluctuating velocity component (m/s), defined in L&S (1998), eq. A-18 p. 516
  wbar,              &! Mean terminal fall velocity of suspended snow (m/s), L&S (1998), eq. A-17 p. 516
  V_r,               &! Ventilation velocity for turbulent suspension (m/s)
  V_rsalt,           &! Ventilation velocity for saltation (m/s)
  xN_r,              &! Particle Reynolds number Re
  xNu,               &! Nusselt number
  xSh,               &! Sherwood number
  tmp1,              &! Auxiliary variable 1
  tmp2,              &! Auxiliary variable 2
  top,               &! Top of eq. A-6, p. 516, L&S (1998)
  bottom              ! Bottom of eq. A-6, p. 516, L&S (1998)

D = 2.06e-5 * (Ta/273.)**(1.75)
! rho_sat = 0.622 * 10.0**(11.40 - 2353./Ta) / (Rair * Ta)
rho_sat = 0.622 / (Rair * Ta) * 610.78 * exp(21.875 * (Ta - Tm) / (Ta - 7.66))

! Assume that the rh varies according to a modification to 
! Pomeroy's humidity variation with height equation.
rh_offset = 1.0 - 0.027 * log(zRH)
sigma = (0.01 * RH - 1.0) * (rh_offset + 0.027 * log(z))
sigma = min(0.0,sigma)
sigma = max(-1.0,sigma)

alpha = 4.08 + 12.6 * z
rbar_r = 4.6e-5 * z**(-0.258)
xmbar = 4.0/3.0 * pi * rho_ice * rbar_r**3 * (1.0 + 3.0/alpha + 2.0/alpha**2)
rbar = ((3.0 * xmbar) / (4.0 * pi * rho_ice))**(0.33)
u_z = Utau/vkman * log(z/z_0)
x_r = 0.005 * u_z**(1.36)
wbar = 1.1e7 * rbar**(1.8)

if (flag_susp) then

  ! Compute the sublimation loss rate coefficient for the suspension
  ! layer.
  V_r = wbar + 3.0 * x_r * cos(pi/4.0)
  xN_r = 2.0 * rbar * V_r / visc_air
  xNu = 1.79 + 0.606 * xN_r**(0.5)
  xSh = xNu
  tmp1 = (Ls * xM)/(Runi * Ta) - 1.0
  tmp2 = hcon_air * Ta * xNu
  top = 2.0 * pi * rbar * sigma
  bottom = Ls/tmp2 * tmp1 + 1.0/(D * rho_sat * xSh)
  V_susp = (top/bottom)/xmbar
  V_salt = 0.0

else

  ! Compute the sublimation loss rate coefficient for the saltation
  ! layer.
  V_rsalt = 0.68 * Utau + 2.3 * Utau_t
  xN_r = 2.0 * rbar * V_rsalt / visc_air
  xNu = 1.79 + 0.606 * xN_r**(0.5)
  xSh = xNu
  tmp1 = (Ls * xM)/(Runi * Ta) - 1.0
  tmp2 = hcon_air * Ta * xNu
  top = 2.0 * pi * rbar * sigma
  bottom = Ls/tmp2 * tmp1 + 1.0/(D * rho_sat * xSh)
  V_salt = (top/bottom)/xmbar
  V_susp = 0.0

end if

end subroutine getsublim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getdirection(index_ue,index_uw,index_vn,index_vs,uwind,vwind)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

implicit none

integer, intent(inout) :: &
  index_ue(Nx,2*Ny+1), &! Wind index array E
  index_uw(Nx,2*Ny+1), &! Wind index array W
  index_vn(Ny,2*Nx+1), &! Wind index array N
  index_vs(Ny,2*Nx+1)   ! Wind index array S

real, intent(in) :: &
  uwind(Nx,Ny),      &! x component of wind speed (m/s)
  vwind(Nx,Ny)        ! y component of wind speed (m/s)

integer :: &
  i,j,               &! Point counters
  npairs              ! ?

real :: &
  sign1,sign2         ! Sign indicators of wind

! Index whether the winds are blowing east or west.  The first
! column of the index array is the number of pairs of begining
! and ending array index of blocks of wind running in the same
! direction.

! Sweep looking for WESTERLY winds, looking for positive numbers.
do i = 1, Nx

  if (uwind(i,1) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 > 0.0) then
    npairs = 1
    index_uw(i,2) =  1
  else
    npairs = 0
  end if

  do j = 2, Ny

    if (uwind(i,j-1) <= 0.0) then
      sign1 = -1.0
    else
      sign1 = 1.0
    end if
    if (uwind(i,j) <= 0.0) then
      sign2 = -1.0
    else
      sign2 = 1.0
    end if

    if (sign2 /= sign1) then
      ! We have a sign change.
      if (sign2 > 0.0) then
        ! We have gone from negative to positive, indicating the start
        ! of a new positive group.
        npairs = npairs + 1
        index_uw(i,npairs*2) = j
      else
        ! We have gone from positive to negative, indicating the end of
        ! the group.
        index_uw(i,npairs*2+1) = j - 1
      end if
    end if

  end do

  if (uwind(i,Ny) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 > 0.0) then
    index_uw(i,npairs*2+1) = Ny
  end if
  index_uw(i,1) = npairs
end do

! Sweep looking for EASTERLY winds, looking for negative numbers.
do i = 1, Nx

  if (uwind(i,1) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 < 0.0) then
    npairs = 1
    index_ue(i,2) = 1
  else
    npairs = 0
  end if

  do j = 2, Ny

    if (uwind(i,j-1) <= 0.0) then
      sign1 = -1.0
    else
      sign1 = 1.0
    end if
    if (uwind(i,j) <= 0.0) then
      sign2 = -1.0
    else
      sign2 = 1.0
    end if

    if (sign2 /= sign1) then
      ! We have a sign change.
      if (sign2 < 0.0) then
        ! We have gone from positive to negative, indicating the start
        ! of a new negative group.
        npairs = npairs + 1
        index_ue(i,npairs*2) = j
      else
        ! We have gone from negative to positive, indicating the end of
        ! the group.
        index_ue(i,npairs*2+1) = j - 1
      end if
    end if

  end do

  if (uwind(i,Ny) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 < 0.0) then
    index_ue(i,npairs*2+1) = Ny
  end if
  index_ue(i,1) = npairs
end do

! Sweep looking for SOUTHERLY winds, looking for positive numbers.
do j = 1, Ny

  if (vwind(1,j) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 > 0.0) then
    npairs = 1
    index_vs(j,2) = 1
  else
    npairs = 0
  end if

  do i = 2, Nx

    if (vwind(i-1,j) <= 0.0) then
      sign1 = -1.0
    else
      sign1 = 1.0
    end if
    if (vwind(i,j) <= 0.0) then
      sign2 = -1.0
    else
      sign2 = 1.0
    end if

    if (sign2 /= sign1) then
      ! We have a sign change.
      if (sign2 > 0.0) then
        ! We have gone from negative to positive, indicating the start
        ! of a new positive group.
        npairs = npairs + 1
        index_vs(j,npairs*2) = i
      else
        ! We have gone from positive to negative, indicating the end of
        ! the group.
        index_vs(j,npairs*2+1) = i - 1
      end if
    end if

  end do

  if (vwind(Nx,j) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 > 0.0) then
    index_vs(j,npairs*2+1) = Nx
  end if
  index_vs(j,1) = npairs
end do

! Sweep looking for NORTHERLY winds, looking for negative numbers.
do j = 1, Ny

  if (vwind(1,j) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 < 0.0) then
    npairs = 1
    index_vn(j,2) = 1
  else
    npairs = 0
  end if

  do i = 2, Nx

    if (vwind(i-1,j) <= 0.0) then
      sign1 = -1.0
    else
      sign1 = 1.0
    end if
    if (vwind(i,j) <= 0.0) then
      sign2 = -1.0
    else
      sign2 = 1.0
    end if

    if (sign2 /= sign1) then
      ! We have a sign change.
      if (sign2 < 0.0) then
        ! We have gone from positive to negative, indicating the start
        ! of a new negative group.
        npairs = npairs + 1
        index_vn(j,npairs*2) = i
      else
        ! We have gone from negative to positive, indicating the end of
        ! the group.
        index_vn(j,npairs*2+1) = i - 1
      end if
    end if

  end do

  if (vwind(Nx,j) <= 0.0) then
    sign1 = -1.0
  else
    sign1 = 1.0
  end if

  if (sign1 < 0.0) then
    index_vn(j,npairs*2+1) = Nx
  end if
  index_vn(j,1) = npairs
end do

end subroutine getdirection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getnewdepth(delta_WE,delta_SN,Qs_u,Qs_v, &
                       dSWE_s, &
                       index_ue,index_uw,index_vn,index_vs, &
                       snowthickness,Ds_soft,snowdepth0,Sice0)

use MODTILE, only: &
  tiled_trans_run              ! Tiled run flag

use DRIVING, only: &
  dt                     ! Timestep (s)

use GRID, only: &
  Nx,Ny                  ! Grid dimensions

use PARAMMAPS, only: &
  veg_shd                ! Vegetation snow holding depth (m)

use STATE_VARIABLES, only: &
  Ds,                   &! Snow layer thicknesses (m)
  fsnow,                &! Snow cover fraction 
  Sice,                 &! Ice content of snow layers (kg/m^2)
  Sliq                   ! Liquid content of snow layers (kg/m^2)

use LANDUSE, only: &
  dem,                  &! Terrain elevation (m)
  forestfrac             ! Forest fraction

use PARAM_SNOWTRAN3D, only: &
  rho_snow               ! Constant snow density (kg/m^3)

implicit none

real, intent(in) :: &
  delta_WE,             &! Grid cell size in WE direction (m)
  delta_SN               ! Grid cell size in SN direction (m)

integer, intent(in) :: &
  index_ue(Nx,2*Ny+1),  &! Wind index array E
  index_uw(Nx,2*Ny+1),  &! Wind index array W
  index_vn(Ny,2*Nx+1),  &! Wind index array N
  index_vs(Ny,2*Nx+1)    ! Wind index array S

real, intent(inout) :: &
  snowthickness(Nx,Ny), &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  Qs_u(Nx,Ny),          &! x component of either Qsalt or Qsusp (kg/m/s)
  Qs_v(Nx,Ny),          &! y component of either Qsalt or Qsusp (kg/m/s)
  dSWE_s(Nx,Ny),        &! SWE change due to either saltation or suspension (kg/m^2)
  Ds_soft(Nx,Ny),       &! Soft snow thickness (m)
  snowdepth0(Nx,Ny),    &! Snow depth of snowdrift accumulation, averaged over the grid cell (m)
  Sice0(Nx,Ny)           ! Surface layer ice content for transported snow (kg/m^2)

real :: &
  dh_s_u(Nx,Ny),        &! x component of snow depth change due to either saltation or suspension (m)
  dh_s_v(Nx,Ny),        &! y component of snow depth change due to either saltation or suspension (m)
  dSWE_s_u(Nx,Ny),      &! x component of SWE change due to either saltation or suspension (kg/m^2)
  dSWE_s_v(Nx,Ny),      &! y component of SWE change due to either saltation or suspension (kg/m^2)
  dSWE_s_u_loss(Nx,Ny), &! x component of SWE loss due to either saltation or suspension (kg/m^2)
  dSWE_s_v_loss(Nx,Ny), &! y component of SWE loss due to either saltation or suspension (kg/m^2)
  dh_s_u_loss(Nx,Ny),   &! x component of snow depth loss due to either saltation or suspension (m)
  dh_s_v_loss(Nx,Ny),   &! y component of snow depth loss due to either saltation or suspension (m)
  dSWE_s_u_gain(Nx,Ny), &! x component of SWE gain due to either saltation or suspension (kg/m^2)
  dSWE_s_v_gain(Nx,Ny)   ! y component of SWE gain due to either saltation or suspension (kg/m^2)

real :: &
  dh_s,                 &! Snow depth change due to either saltation or suspension (m)
  dh_s_u_gain,          &! x component of snow depth gain due to either saltation or suspension (m)
  dh_s_v_gain,          &! y component of snow depth gain due to either saltation or suspension (m)
  dh_s_gain,            &! Snow depth gain due to either saltation or suspension (m)
  dh_s_loss,            &! Snow depth loss due to either saltation or suspension (m)
  dSWE_s_loss,          &! SWE loss due to either saltation or suspension (kg/m^2)
  dSWE_s_gain,          &! SWE gain due to either saltation or suspension (kg/m^2)
  dSWE_s_u_loss_tmp,    &! Temporary variable for x component of SWE loss due to either saltation or suspension (kg/m^2)
  dSWE_s_v_loss_tmp,    &! Temporary variable for y component of SWE loss due to either saltation or suspension (kg/m^2)
  dh_s_u_loss_tmp,      &! Temporary variable for x component of snow depth loss due to either saltation or suspension (m)
  dh_s_v_loss_tmp,      &! Temporary variable for y component of snow depth loss due to either saltation or suspension (m)
  Ds_hard,              &! Hard snow thickness (m)
  eps,                  &! Epsilon
  weight_u,             &! x weight for dSWE_s
  weight_v,             &! y weight for dSWE_s
  snowdmin,             &! Minimum snow depth to allow transport (m)
  swe_loc                ! SWE in the snow covered part of the pixel (kg/m^2)

integer :: &
  i,j,k,                &! Point counters
  istart,iend,          &! Point counters boundaries
  jstart,jend            ! Point counters boundaries  

! Define an upwind boundary condition for saltation (here I have
! assumed that the transport is in equilibrium).
dh_s_u(:,:) = 0.0
dh_s_v(:,:) = 0.0
dSWE_s_u(:,:) = 0.0
dSWE_s_v(:,:) = 0.0
dSWE_s_u_gain(:,:) = 0.0
dSWE_s_v_gain(:,:) = 0.0
dSWE_s_u_loss(:,:) = 0.0
dSWE_s_v_loss(:,:) = 0.0
dh_s_u_loss(:,:) = 0.0
dh_s_v_loss(:,:) = 0.0

! Consider WESTERLY winds.
do i = 1, Nx
  do k = 1, index_uw(i,1)
    jstart = index_uw(i,k*2)+1
    jend = index_uw(i,k*2+1)
    do j = jstart, jend
      swe_loc = sum(Sice(:,i,j) + Sliq(:,i,j))
      dSWE_s_u_loss(i,j) = dt * Qs_u(i,j) * fsnow(i,j) / delta_WE
      dSWE_s_u_gain(i,j) = dt * Qs_u(i,j-1) * fsnow(i,j-1) / delta_WE
      ! The SWE loss due to wind only occurs on the open part of the pixel.
      ! The SWE gained by neighbouring is then weighted by the open fraction of the pixel.
      ! The SWE gain due to wind only occurs on the open part of the pixel.
      ! If it is a tiled run or open run, the transferred SWE is then weighted differently.
      if (tiled_trans_run) then
        ! In the open tile, SWE loss stays the same since it will be weighted later when combining the tiles.
        ! SWE gain from neighbouring pixel is first weighted by the open fraction of origin pixel, 
        ! then divided by the open fraction of the target pixel in the open tile, 
        ! since it will be counter-weighted later when combining the tiles.
        if (forestfrac(i,j) <= 0.9) then
          dSWE_s_u_gain(i,j) = dSWE_s_u_gain(i,j) * (1.0 - forestfrac(i,j-1)) / (1.0 - forestfrac(i,j))
          ! else: the flux is already zero
        end if
      else
        dSWE_s_u_loss(i,j) = dSWE_s_u_loss(i,j) * (1.0 - forestfrac(i,j))
        dSWE_s_u_gain(i,j) = dSWE_s_u_gain(i,j) * (1.0 - forestfrac(i,j-1))
      end if
      dSWE_s_u_loss(i,j) = min(dSWE_s_u_loss(i,j),swe_loc) ! No need to adjust Qs_u here because if thresholded, it will be done in the next if loop anyway
      if (dSWE_s_u_loss(i,j) > epsilon(dSWE_s_u_loss)) then
        dSWE_s_u_loss_tmp = dSWE_s_u_loss(i,j)
        call HS_FROM_SWE(dSWE_s_u_loss_tmp,dh_s_u_loss_tmp,i,j)
        dh_s_u_loss(i,j) = dh_s_u_loss_tmp
      else
        dSWE_s_u_loss(i,j) = 0.0
        dh_s_u_loss(i,j) = 0.0
      end if
      dh_s_u_gain = dSWE_s_u_gain(i,j) / rho_snow
      dh_s_u(i,j) = dh_s_u_gain - dh_s_u_loss(i,j)
      dSWE_s_u(i,j) = dSWE_s_u_gain(i,j) - dSWE_s_u_loss(i,j)

      ! Make adjustments for the case where there is no snow available
      ! on the ground (or captured within the vegetation) to be
      ! eroded.
      Ds_hard = snowthickness(i,j) - Ds_soft(i,j)
      snowdmin = max(veg_shd(i,j),Ds_hard)
      if (snowthickness(i,j) > snowdmin .and. fsnow(i,j) > epsilon(fsnow)) then
        if (snowthickness(i,j) - dh_s_u_loss(i,j) / fsnow(i,j) <= snowdmin) then
          dh_s_u_loss(i,j) = (snowthickness(i,j) - snowdmin) * fsnow(i,j)
          if (dh_s_u_loss(i,j) > epsilon(dh_s_u_loss)) then
            dh_s_u_loss_tmp = dh_s_u_loss(i,j)
            call SWE_FROM_HS(dh_s_u_loss_tmp,dSWE_s_u_loss_tmp,i,j)
            dSWE_s_u_loss(i,j) = dSWE_s_u_loss_tmp
            ! Same open tile weighting as before.
            if (.not. tiled_trans_run) then
              dSWE_s_u_loss(i,j) = dSWE_s_u_loss(i,j) * (1.0 - forestfrac(i,j))
            end if
          else
            dh_s_u_loss(i,j) = 0.0
            dSWE_s_u_loss(i,j) = 0.0
          end if
          dh_s_u(i,j) = dh_s_u_gain - dh_s_u_loss(i,j)
          dSWE_s_u(i,j) = dSWE_s_u_gain(i,j) - dSWE_s_u_loss(i,j)
          Qs_u(i,j) = Qs_u(i,j-1) - dSWE_s_u(i,j) * delta_WE / dt / fsnow(i,j)
        end if
      else
        Qs_u(i,j) = 0.0
        dh_s_u(i,j) = dh_s_u_gain
        dSWE_s_u_loss(i,j) = 0.0
        dh_s_u_loss(i,j) = 0.0
        dSWE_s_u(i,j) = dSWE_s_u_gain(i,j)
      end if
    end do
  end do
end do

! Consider EASTERLY winds.
do i = 1, Nx
  do k = 1, index_ue(i,1)
    jend = index_ue(i,k*2)
    jstart = index_ue(i,k*2+1)-1
    do j = jstart, jend,-1
      swe_loc = sum(Sice(:,i,j) + Sliq(:,i,j))
      dSWE_s_u_loss(i,j) = dt * Qs_u(i,j) * fsnow(i,j) / delta_WE
      dSWE_s_u_gain(i,j) = dt * Qs_u(i,j+1) * fsnow(i,j+1) / delta_WE
      ! The SWE loss due to wind only occurs on the open part of the pixel.
      ! The SWE gained by neighbouring is then weighted by the open fraction of the pixel.
      ! The SWE gain due to wind only occurs on the open part of the pixel.
      ! If it is a tiled run or open run, the transferred SWE is then weighted differently.
      if (tiled_trans_run) then
        ! In the open tile, SWE loss stays the same since it will be weighted later when combining the tiles.
        ! SWE gain from neighbouring pixel is first weighted by the open fraction of origin pixel, 
        ! then divided by the open fraction of the target pixel in the open tile, 
        ! since it will be counter-weighted later when combining the tiles.
        if (forestfrac(i,j) <= 0.9) then
          dSWE_s_u_gain(i,j) = dSWE_s_u_gain(i,j) * (1.0 - forestfrac(i,j+1)) / (1.0 - forestfrac(i,j))
          ! else: the flux is already zero
        end if
      else
        dSWE_s_u_loss(i,j) = dSWE_s_u_loss(i,j) * (1.0 - forestfrac(i,j))
        dSWE_s_u_gain(i,j) = dSWE_s_u_gain(i,j) * (1.0 - forestfrac(i,j+1))
      end if
      dSWE_s_u_loss(i,j) = min(dSWE_s_u_loss(i,j),swe_loc) ! No need to adjust Qs_u here because if thresholded, it will be done in the next if loop anyway
      if (dSWE_s_u_loss(i,j) > epsilon(dSWE_s_u_loss)) then
        dSWE_s_u_loss_tmp = dSWE_s_u_loss(i,j)
        call HS_FROM_SWE(dSWE_s_u_loss_tmp,dh_s_u_loss_tmp,i,j)
        dh_s_u_loss(i,j) = dh_s_u_loss_tmp
      else
        dSWE_s_u_loss(i,j) = 0.0
        dh_s_u_loss(i,j) = 0.0
      end if
      dh_s_u_gain = dSWE_s_u_gain(i,j) / rho_snow
      dh_s_u(i,j) = dh_s_u_gain - dh_s_u_loss(i,j)
      dSWE_s_u(i,j) = dSWE_s_u_gain(i,j) - dSWE_s_u_loss(i,j)

      ! Make adjustments for the case where there is no snow available
      ! on the ground (or captured within the vegetation) to be
      ! eroded.
      Ds_hard = snowthickness(i,j) - Ds_soft(i,j)
      snowdmin = max(veg_shd(i,j),Ds_hard)
      if (snowthickness(i,j) > snowdmin .and. fsnow(i,j) > epsilon(fsnow)) then
        if (snowthickness(i,j) - dh_s_u_loss(i,j) / fsnow(i,j) <= snowdmin) then
          dh_s_u_loss(i,j) = (snowthickness(i,j) - snowdmin) * fsnow(i,j)
          if (dh_s_u_loss(i,j) > epsilon(dh_s_u_loss)) then
            dh_s_u_loss_tmp = dh_s_u_loss(i,j)
            call SWE_FROM_HS(dh_s_u_loss_tmp,dSWE_s_u_loss_tmp,i,j)
            dSWE_s_u_loss(i,j) = dSWE_s_u_loss_tmp
            ! Same open tile weighting as before.
            if (.not. tiled_trans_run) then
              dSWE_s_u_loss(i,j) = dSWE_s_u_loss(i,j) * (1.0 - forestfrac(i,j))
            end if
          else
            dh_s_u_loss(i,j) = 0.0
            dSWE_s_u_loss(i,j) = 0.0
          end if
          dh_s_u(i,j) = dh_s_u_gain - dh_s_u_loss(i,j)
          dSWE_s_u(i,j) = dSWE_s_u_gain(i,j) - dSWE_s_u_loss(i,j)
          Qs_u(i,j) = Qs_u(i,j+1) - dSWE_s_u(i,j) * delta_WE / dt / fsnow(i,j)
        end if
      else
        Qs_u(i,j) = 0.0
        dh_s_u(i,j) = dh_s_u_gain
        dSWE_s_u_loss(i,j) = 0.0
        dh_s_u_loss(i,j) = 0.0
        dSWE_s_u(i,j) = dSWE_s_u_gain(i,j)
      end if
    end do
  end do
end do

! Consider SOUTHERLY winds.
do j = 1, Ny
  do k = 1, index_vs(j,1)
    istart = index_vs(j,k*2)+1
    iend = index_vs(j,k*2+1)
    do i = istart, iend
      swe_loc = sum(Sice(:,i,j) + Sliq(:,i,j))
      dSWE_s_v_loss(i,j) = dt * Qs_v(i,j) * fsnow(i,j) / delta_SN
      dSWE_s_v_gain(i,j) = dt * Qs_v(i-1,j) * fsnow(i-1,j) / delta_SN
      ! The SWE loss due to wind only occurs on the open part of the pixel.
      ! The SWE gained by neighbouring is then weighted by the open fraction of the pixel.
      ! The SWE gain due to wind only occurs on the open part of the pixel.
      ! If it is a tiled run or open run, the transferred SWE is then weighted differently.
      if (tiled_trans_run) then
        ! In the open tile, SWE loss stays the same since it will be weighted later when combining the tiles.
        ! SWE gain from neighbouring pixel is first weighted by the open fraction of origin pixel, 
        ! then divided by the open fraction of the target pixel in the open tile, 
        ! since it will be counter-weighted later when combining the tiles.
        if (forestfrac(i,j) <= 0.9) then
          dSWE_s_v_gain(i,j) = dSWE_s_v_gain(i,j) * (1.0 - forestfrac(i-1,j)) / (1.0 - forestfrac(i,j))
          ! else: the flux is already zero
        end if
      else
        dSWE_s_v_loss(i,j) = dSWE_s_v_loss(i,j) * (1.0 - forestfrac(i,j))
        dSWE_s_v_gain(i,j) = dSWE_s_v_gain(i,j) * (1.0 - forestfrac(i-1,j))
      end if
      dSWE_s_v_loss(i,j) = min(dSWE_s_v_loss(i,j),swe_loc) ! No need to adjust Qs_v here because if thresholded, it will be done in the next if loop anyway
      if (dSWE_s_v_loss(i,j) > epsilon(dSWE_s_v_loss)) then
        dSWE_s_v_loss_tmp = dSWE_s_v_loss(i,j)
        call HS_FROM_SWE(dSWE_s_v_loss_tmp,dh_s_v_loss_tmp,i,j)
        dh_s_v_loss(i,j) = dh_s_v_loss_tmp
      else
        dSWE_s_v_loss(i,j) = 0.0
        dh_s_v_loss(i,j) = 0.0
      end if
      dh_s_v_gain = dSWE_s_v_gain(i,j) / rho_snow
      dh_s_v(i,j) = dh_s_v_gain - dh_s_v_loss(i,j)
      dSWE_s_v(i,j) = dSWE_s_v_gain(i,j) - dSWE_s_v_loss(i,j)

      ! Make adjustments for the case where there is no snow available
      ! on the ground (or captured within the vegetation) to be
      ! eroded.
      Ds_hard = snowthickness(i,j) - Ds_soft(i,j)
      snowdmin = max(veg_shd(i,j),Ds_hard)
      if (snowthickness(i,j) > snowdmin .and. fsnow(i,j) > epsilon(fsnow)) then
        if (snowthickness(i,j) - dh_s_v_loss(i,j) / fsnow(i,j) <= snowdmin) then
          dh_s_v_loss(i,j) = (snowthickness(i,j) - snowdmin) * fsnow(i,j)
          if (dh_s_v_loss(i,j) > epsilon(dh_s_v_loss)) then
            dh_s_v_loss_tmp = dh_s_v_loss(i,j)
            call SWE_FROM_HS(dh_s_v_loss_tmp,dSWE_s_v_loss_tmp,i,j)
            dSWE_s_v_loss(i,j) = dSWE_s_v_loss_tmp
            ! Same open tile weighting as before.
            if (.not. tiled_trans_run) then
              dSWE_s_v_loss(i,j) = dSWE_s_v_loss(i,j) * (1.0 - forestfrac(i,j))
            end if
          else
            dh_s_v_loss(i,j) = 0.0
            dSWE_s_v_loss(i,j) = 0.0
          end if
          dh_s_v(i,j) = dh_s_v_gain - dh_s_v_loss(i,j)
          dSWE_s_v(i,j) = dSWE_s_v_gain(i,j) - dSWE_s_v_loss(i,j)
          Qs_v(i,j) = Qs_v(i-1,j) - dSWE_s_v(i,j) * delta_SN / dt / fsnow(i,j)
        end if
      else
        Qs_v(i,j) = 0.0
        dh_s_v(i,j) = dh_s_v_gain
        dSWE_s_v_loss(i,j) = 0.0
        dh_s_v_loss(i,j) = 0.0
        dSWE_s_v(i,j) = dSWE_s_v_gain(i,j)
      end if
    end do
  end do
end do

! Consider NORTHERLY winds.
do j = 1, Ny
  do k = 1, index_vn(j,1)
    iend = index_vn(j,k*2)
    istart = index_vn(j,k*2+1)-1
    do i = istart, iend,-1
      swe_loc = sum(Sice(:,i,j) + Sliq(:,i,j))
      dSWE_s_v_loss(i,j) = dt * Qs_v(i,j) * fsnow(i,j) / delta_SN
      dSWE_s_v_gain(i,j) = dt * Qs_v(i+1,j) * fsnow(i+1,j) / delta_SN
      ! The SWE loss due to wind only occurs on the open part of the pixel.
      ! The SWE gained by neighbouring is then weighted by the open fraction of the pixel.
      ! The SWE gain due to wind only occurs on the open part of the pixel.
      ! If it is a tiled run or open run, the transferred SWE is then weighted differently.
      if (tiled_trans_run) then
        ! In the open tile, SWE loss stays the same since it will be weighted later when combining the tiles.
        ! SWE gain from neighbouring pixel is first weighted by the open fraction of origin pixel, 
        ! then divided by the open fraction of the target pixel in the open tile, 
        ! since it will be counter-weighted later when combining the tiles.
        if (forestfrac(i,j) <= 0.9) then
          dSWE_s_v_gain(i,j) = dSWE_s_v_gain(i,j) * (1.0 - forestfrac(i+1,j)) / (1.0 - forestfrac(i,j))
          ! else: the flux is already zero
        end if
      else
        dSWE_s_v_loss(i,j) = dSWE_s_v_loss(i,j) * (1.0 - forestfrac(i,j))
        dSWE_s_v_gain(i,j) = dSWE_s_v_gain(i,j) * (1.0 - forestfrac(i+1,j))
      end if
      dSWE_s_v_loss(i,j) = min(dSWE_s_v_loss(i,j),swe_loc) ! No need to adjust Qs_v here because if thresholded, it will be done in the next if loop anyway
      if (dSWE_s_v_loss(i,j) > epsilon(dSWE_s_v_loss)) then
        dSWE_s_v_loss_tmp = dSWE_s_v_loss(i,j)
        call HS_FROM_SWE(dSWE_s_v_loss_tmp,dh_s_v_loss_tmp,i,j)
        dh_s_v_loss(i,j) = dh_s_v_loss_tmp
      else
        dSWE_s_v_loss(i,j) = 0.0
        dh_s_v_loss(i,j) = 0.0
      end if
      dh_s_v_gain = dSWE_s_v_gain(i,j) / rho_snow
      dh_s_v(i,j) = dh_s_v_gain - dh_s_v_loss(i,j)
      dSWE_s_v(i,j) = dSWE_s_v_gain(i,j) - dSWE_s_v_loss(i,j)

      ! Make adjustments for the case where there is no snow available
      ! on the ground (or captured within the vegetation) to be
      ! eroded.
      Ds_hard = snowthickness(i,j) - Ds_soft(i,j)
      snowdmin = max(veg_shd(i,j),Ds_hard)
      if (snowthickness(i,j) > snowdmin .and. fsnow(i,j) > epsilon(fsnow)) then
        if (snowthickness(i,j) - dh_s_v_loss(i,j) / fsnow(i,j) <= snowdmin) then
          dh_s_v_loss(i,j) = (snowthickness(i,j) - snowdmin) * fsnow(i,j)
          if (dh_s_v_loss(i,j) > epsilon(dh_s_v_loss)) then
            dh_s_v_loss_tmp = dh_s_v_loss(i,j)
            call SWE_FROM_HS(dh_s_v_loss_tmp,dSWE_s_v_loss_tmp,i,j)
            dSWE_s_v_loss(i,j) = dSWE_s_v_loss_tmp
            ! Same open tile weighting as before.
            if (.not. tiled_trans_run) then
              dSWE_s_v_loss(i,j) = dSWE_s_v_loss(i,j) * (1.0 - forestfrac(i,j))
            end if
          else
            dh_s_v_loss(i,j) = 0.0
            dSWE_s_v_loss(i,j) = 0.0
          end if
          dh_s_v(i,j) = dh_s_v_gain - dh_s_v_loss(i,j)
          dSWE_s_v(i,j) = dSWE_s_v_gain(i,j) - dSWE_s_v_loss(i,j)
          Qs_v(i,j) = Qs_v(i+1,j) - dSWE_s_v(i,j) * delta_SN / dt / fsnow(i,j)
        end if
      else
        Qs_v(i,j) = 0.0
        dh_s_v(i,j) = dh_s_v_gain
        dSWE_s_v_loss(i,j) = 0.0
        dh_s_v_loss(i,j) = 0.0
        dSWE_s_v(i,j) = dSWE_s_v_gain(i,j)
      end if
    end do
  end do
end do

! Update the snow depth changes due to saltation or suspension transport from the
! the east and west, and north and south.
eps = 1e-6
do j = 1, Ny
  do i = 1, Nx

    ! LQ: I can't understand the meaning of this weighting, I keep it as original SNOWTRAN3D though
    weight_u = abs(dh_s_u(i,j)) / &
               (abs(dh_s_u(i,j)) + abs(dh_s_v(i,j)) + eps)

    weight_v = abs(dh_s_v(i,j)) / &
               (abs(dh_s_u(i,j)) + abs(dh_s_v(i,j)) + eps)

    dSWE_s_u(i,j) = weight_u * dSWE_s_u(i,j)
    dSWE_s_v(i,j) = weight_v * dSWE_s_v(i,j)
    dSWE_s_u_gain(i,j) = weight_u * dSWE_s_u_gain(i,j)
    dSWE_s_v_gain(i,j) = weight_v * dSWE_s_v_gain(i,j)
    dSWE_s_u_loss(i,j) = weight_u * dSWE_s_u_loss(i,j)
    dSWE_s_v_loss(i,j) = weight_v * dSWE_s_v_loss(i,j)

    dSWE_s(i,j) = dSWE_s_u(i,j) + dSWE_s_v(i,j)
    dSWE_s_gain = dSWE_s_u_gain(i,j) + dSWE_s_v_gain(i,j)
    dSWE_s_loss = dSWE_s_u_loss(i,j) + dSWE_s_v_loss(i,j)
    dh_s_gain = dSWE_s_gain / rho_snow
    if (dSWE_s_loss > epsilon(dSWE_s_loss)) then
      call HS_FROM_SWE(dSWE_s_loss,dh_s_loss,i,j)
    else
      dSWE_s_loss = 0.0
      dh_s_loss = 0.0
    end if
    dh_s = dh_s_gain - dh_s_loss

    Ds_soft(i,j) = max(Ds_soft(i,j) - dh_s_loss,0.0)

    if (dSWE_s_loss > epsilon(dSWE_s_loss) .and. dh_s_loss > epsilon(dh_s_loss)) then
      call SNOW_ABLATION(dh_s_loss,dSWE_s_loss,i,j)
    end if

    ! Net mass gain for this grid cell at this time step.
    if (dSWE_s_gain > epsilon(dSWE_s_gain) .and. dh_s_gain > epsilon(dh_s_gain)) then

      ! Add to the existing top layer.
      Sice0(i,j) = Sice0(i,j) + dSWE_s_gain
      snowdepth0(i,j) = snowdepth0(i,j) + dh_s_gain

    end if

    ! Update the snow layer thicknesses
    snowthickness(i,j) = sum(Ds(:,i,j))

  end do
end do

end subroutine getnewdepth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine surface_snow(Utau_t)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use LANDUSE, only: &
  dem                 ! Terrain elevation (m)

use PARAM_SNOWTRAN3D, only: &
  rho_snow            ! Constant snow density (kg/m^3)

use STATE_VARIABLES, only: &
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Ds                  ! Snow layer thicknesses (m)

implicit none

real, intent(inout) :: &
  Utau_t(Nx,Ny)       ! Threshold friction velocity (m/s)

integer :: &
  i,j                 ! Point counters

real :: &
  rho_surf_snow       ! Density of top snow layer  (kg/m^3)


do j = 1, Ny
  do i = 1, Nx

    rho_surf_snow = rho_snow
    if (Nsnow(i,j) > 0 .and. Ds(1,i,j) > epsilon(Ds)) then
      rho_surf_snow = (Sice(1,i,j) + Sliq(1,i,j)) / Ds(1,i,j) / fsnow(i,j)
    end if

    ! Calculate the snow threshold friction velocity.
    if (rho_surf_snow <= 300.0) then
      Utau_t(i,j) = 0.10 * exp(0.003 * rho_surf_snow)
    else
      Utau_t(i,j) = 0.005 * exp(0.013 * rho_surf_snow)
    end if

  end do
end do

end subroutine surface_snow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!