!-----------------------------------------------------------------------
! Surface exchange coefficients
!-----------------------------------------------------------------------
subroutine SFEXCH(gs1,KH,KHa,KHg,KHv,KWg,KWv,Usc)

use, intrinsic :: iso_fortran_env, only: dp=>real64

use MODCONF, only: CANMOD, ZOFFST, EXCHNG, OSHDTN, SNFRAC

use MODPERT, only: Z0PERT

use MODTILE, only: TILE, tthresh

use CONSTANTS, only: &
  grav,              &! Acceleration due to gravity (m/s^2)
  vkman               ! Von Karman constant

use DRIVING, only: &
  Ta,                &! Air temperature (K)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Ua,                &! Wind speed (m/s)
  zT,                &! Temperature and humidity measurement height (m)
  zU,                &! Wind measurement height (m)
  z0P                 ! perturbed value for z0
use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only : &
  bstb,              &! Atmospheric stability parameter
  cden,              &! Dense canopy turbulent transfer coefficient
  cveg,              &! Vegetation turbulent transfer coefficient
  gsnf,              &! Snow-free vegetation moisture conductance (m/s)
  rchd,              &! Ratio of displacement height to canopy height
  rchz,              &! Ratio of roughness length to canopy height
  z0sn,              &! Snow roughness length (m)
  wcan,              &! Canopy wind decay coefficient
  zsub,              &! Sub-canopy reference height (m)
  zgf,               &! z0g canopy dependence factors and ranges 
  zgr,               &! z0g canopy dependence range
  khcf                ! diffusivity correction factor 

use PARAMMAPS, only: &
  VAI,               &! Vegetation area index
  z0sf                ! Snow-free surface roughness length (m)

use STATE_VARIABLES, only : &

  Qcan,              &! Canopy air space humidity
  fsnow,             &! Snow cover fraction 
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sveg,              &! Snow mass on vegetation (kg/m^2)
  Tcan,              &! Canopy air space temperature (K)
  Tsrf,              &! Surface skin temperature (K)
  Tveg,              &! Vegetation temperature (K)
  Ds

use LANDUSE, only: &
  dem,               &! Terrain elevation (m) 
  fveg,              &! Canopy cover fraction
  fves,              &! Stand-scale vegetation fraction
  hcan,              &! Canopy height (m)
  tilefrac            ! Grid cell tile fraction

implicit none

real*8, intent(in) :: &
  gs1(Nx,Ny)          ! Surface moisture conductance (m/s)

real*8, intent(out) :: &
  KH(Nx,Ny),         &! Eddy diffusivity for heat to the atmosphere (m/s)
  KHa(Nx,Ny),        &! Eddy diffusivity from the canopy air space (m/s)
  KHg(Nx,Ny),        &! Eddy diffusivity for heat from the ground (m/s)
  KHv(Nx,Ny),        &! Eddy diffusivity for heat from vegetation (m/s)
  KWg(Nx,Ny),        &! Eddy diffusivity for water from the ground (m/s)
  KWv(Nx,Ny),        &! Eddy diffusivity for water from vegetation (m/s)
  Usc(Nx,Ny)           ! Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)

integer :: &
  i,j                 ! Point counters

real*8 :: &
  CD,                &! Drag coefficient
  dh,                &! Displacement height (m)
  fh,                &! Stability factor
  Qs,                &! Saturation humidity
  RiB,               &! Bulk Richardson number
  Ric,               &! Sub-canopy Richardson number
  Tint,              &! Interpolated canopy - ground temperature (K)
  ustar,             &! Friction velocity (m/s)
  zT1,               &! Temperature measurement height with offset (m)
  zU1,               &! Wind measurement height with offset (m)
  z0,                &! Roughness length for momentum (m)
  z0g,               &! Ground surface roughness length (m)
  z0h,               &! Roughness length for heat (m)
  z0v,               &! Vegetation roughness length (m)
  Uc,                &! Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)
  Uso,               &! Windspeed at sub-canopy reference height - open (m/s)
  Usf,               &! Windspeed at sub-canopy reference height - forest (m/s)
  Usub,              &! Wind speed at sub-canopy reference height (m/s)
  Uh,                &! Wind velocity at canopy top (m/s)
  KHh,               &! Eddy diffusivity at canopy top (m/s)
  rad,               &! Aerodynamic resistance between canopy air space and atmosphere in dense canopy, at reference height (s/m)
  rgd,               &! Aerodynamic resistance between canopy air space and atmosphere in dense canopy, at reference height (s/m)
  rgo,               &! Theoretical ground aerodynamic resistance for an open site (s/m) 
  z0loc               ! Elevation-dependent roughness length of snow (m), constant or tuned depending on OSHDTN

do j = 1, Ny
do i = 1, Nx

  if (tilefrac(i,j) < tthresh) goto 1 ! exclude points outside tile of interest

  if (ZOFFST == 0) then
    ! Heights specified above ground
    zU1 = zU
    zT1 = zT
  else ! ZOFFST == 1
    ! Heights specified above canopy top
    zU1 = zU + hcan(i,j)
    zT1 = zT + hcan(i,j)
  endif
  
  if (Z0PERT .eqv. .TRUE.) then
    z0g = z0P(i,j)
  else
    ! Ground roughness length
    if (TILE == 'glacier') then
      z0loc = 0.0009_dp
    else
      if (OSHDTN == 0 .OR. tile == 'forest') then
        z0loc = z0sn
      else ! OSHDTN == 1
        if (dem(i,j) >= 2300) then
          z0loc = 0.003_dp
        else if (dem(i,j) >= 1500_dp) then
          z0loc = 0.03_dp + (dem(i,j) - 1500_dp) / (2300_dp - 1500_dp) * (0.003_dp - 0.03_dp)
        else if (dem(i,j) >= 1200_dp) then  ! simple linear b/w two above values
          z0loc = 0.2_dp + (dem(i,j) - 1200_dp) / (1500_dp - 1200_dp) * (0.03_dp - 0.2_dp)
        else
          z0loc = 0.2_dp
        end if
      end if
    end if
    z0g = z0loc
  end if

  ! BC, stabilize the tuning point runs by using a Ds threshold instead of fsnow.
  ! TODO: test the impact for the grid points.
  if (SNFRAC==3) then
    if (sum(Ds(:,i,j)) <= 0.05_dp) then
      z0g = z0sf(i,j)
    endif
  else
    if (fsnow(i,j) <= epsilon(fsnow)) z0g = z0sf(i,j)
  endif

  ! Additional roughness lengths and friction velocity
  if (EXCHNG == 2) then ! Forest - specific adjustment *GM
    ! Open
    z0h = 0.1_dp * z0g 
    ustar = vkman*Ua(i,j)/log(zU/z0g) 
    rgo = log(zT/z0h)/(vkman*ustar) 

    ! Forest
    if (fveg(i,j) > epsilon(fveg(i,j))) then
      z0g = (zgf + zgr*fveg(i,j)) * z0g    
      z0h = 0.1_dp * z0g
      dh = rchd * hcan(i,j) 
      z0v = rchz * hcan(i,j) 
      ustar = vkman*Ua(i,j)/log((zU1 - dh)/z0v)
      Uh = (ustar/vkman)*log((hcan(i,j) - dh)/z0v)
      KHh = vkman*ustar*(hcan(i,j) - dh)
      Usf = exp(wcan*(zsub/hcan(i,j) - 1_dp))*Uh
    end if
  else
    z0v = rchz*hcan(i,j)
    z0  = (z0v**fveg(i,j)) * (z0g**(1_dp - fveg(i,j)))
    z0h = 0.1_dp*z0
    dh = fveg(i,j)*rchd*hcan(i,j)
    CD = (vkman / log((zU1 - dh)/z0))**2_dp
    ustar = sqrt(CD)*Ua(i,j)
  endif
  Uso = Ua(i,j)*log(zsub/z0g)/log(zU/z0g)

  if (EXCHNG == 0) then
  ! No stability adjustment
    fh = 1_dp
    Ric = 0_dp
  endif
  if (EXCHNG == 1) then
  ! Stability adjustment (Louis et al. 1982, quoted by Beljaars 1992)
    Tint = fveg(i,j)*Tveg(i,j) + (1_dp - fveg(i,j))*Tsrf(i,j)
    RiB = grav*(Ta(i,j) - Tint)*(zU1 - dh)**2_dp / ((zT1 - dh)*Ta(i,j)*Ua(i,j)**2_dp)
    if (RiB > 0.2_dp) RiB = 0.2_dp ! New maximum threshold for RiB
    if (RiB > 0_dp) then 
      fh = 1_dp/(1_dp + 3_dp*bstb*RiB*sqrt(1_dp + bstb*RiB))
    else
      fh = 1_dp - 3_dp*bstb*RiB / (1_dp + 3_dp*bstb**2_dp*CD*sqrt(-RiB*zU1/z0))
    endif
    Ric = grav*(Tcan(i,j) - Tsrf(i,j))*hcan(i,j) / (Tcan(i,j)*ustar**2_dp)
    Ric = max(min(Ric,10.0_dp),0.0_dp)
  endif
  ! Note that currently, fh and Ric are not used in EXCHNG == 2, i.e. no stability correction

  ! Eddy diffusivities
  if (fveg(i,j) == 0_dp) then
    KH(i,j) = fh*vkman*ustar / log(zT1/z0h)
    call QSAT(Ps(i,j),Tsrf(i,j),Qs)
    if (Sice(1,i,j) > epsilon(Sice(1,i,j))  .or. Qa(i,j) > Qs) then
      KWg(i,j) = KH(i,j)
    else
      KWg(i,j) = gs1(i,j)*KH(i,j) / (gs1(i,j) + KH(i,j))
    endif
    Usc(i,j) = Uso
  else
    if (EXCHNG == 2) then
      rad = (log((zT1 - dh)/(hcan(i,j)- dh))/(vkman*ustar) +  & 
        hcan(i,j)*(exp(wcan*(1_dp -(z0v + dh)/hcan(i,j))) - 1_dp)/(wcan*KHh))/khcf
      KHa(i,j) = sqrt(fves(i,j))/ rad 
      Usub = sqrt(fves(i,j))*Usf + (1_dp - sqrt(fves(i,j)))*Uso
      Usub = max(Usub, 0.1_dp)
      rgd = 1_dp / (vkman**2_dp * Usub) * log(zsub/z0h) * log(zsub/z0g); 
      KHg(i,j)  = 1_dp/rgd
      Uc = exp(wcan*((z0v + dh)/hcan(i,j) - 1_dp))*Uh
      KHv(i,j) = VAI(i,j)*sqrt(Uc)/cveg
      Usc(i,j) = Usub
    else
      KHa(i,j) = fh*vkman*ustar / log((zT1 - dh)/z0)
      KHg(i,j) = vkman*ustar*((1_dp - fveg(i,j))*fh/log(z0/z0h) + fveg(i,j)*cden/(1_dp + 0.5_dp*Ric))
      KHv(i,j) = sqrt(ustar)*VAI(i,j)/cveg
    endif
    call QSAT(Ps(i,j),Tsrf(i,j),Qs)
    if (Qcan(i,j) > Qs) then
      KWg(i,j) = KHg(i,j)
    else
      KWg(i,j) = gs1(i,j)*KHg(i,j) / (gs1(i,j) + KHg(i,j))
    endif
    call QSAT(Ps(i,j),Tveg(i,j),Qs)
    if (Sveg(i,j) > epsilon(Sveg(i,j)) .or. Qcan(i,j) > Qs) then
      KWv(i,j) = KHv(i,j)
    else
      KWv(i,j) = gsnf*KHv(i,j) / (gsnf + KHv(i,j))
    endif
  
    if (CANMOD == 0) then
    ! Combined resistances for 0-layer canopy model
      KH(i,j) = KHg(i,j)*(KHa(i,j) + KHv(i,j)) / (KHa(i,j) + KHg(i,j) + KHv(i,j))
      KWg(i,j) = KWg(i,j)*(KHa(i,j) + KWv(i,j)) / (KHa(i,j) + KWg(i,j) + KWv(i,j))
    end if
  end if

  1 continue
  
end do
end do

end subroutine SFEXCH
