!-----------------------------------------------------------------------
! Surface energy balance in open areas or zero-layer forest canopy model
!-----------------------------------------------------------------------
subroutine EBALSRF(Ds1,KH,KHa,KHv,KWg,KWv,ks1,SWsrf,SWveg,Ts1, &
                   Esrf,Eveg,G,H,Hsrf,LE,LEsrf,Melt,Rnet,Rsrf, &
                   LWsci,LWveg)

use, intrinsic :: iso_fortran_env, only: dp=>real64

use MODCONF, only: CANMOD

use MODTILE, only: TILE, tthresh 

use CONSTANTS, only: &
  cp,                &! Specific heat capacity of air (J/K/kg)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  Lv,                &! Latent heat of vapourisation (J/kg)
  Rair,              &! Gas constant for air (J/K/kg)
  Rwat,              &! Gas constant for water vapour (J/K/kg)
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm                  ! Melting point (K)

use DRIVING, only: &
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Ta                  ! Air temperature (K)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMMAPS, only: &
  trcn                ! Canopy transmissivity

use STATE_VARIABLES, only: &
  Sice,              &! Ice content of snow layers (kg/m^2)
  Tcan,              &! Canopy air space temperature (K)
  Tsrf,              &! Surface temperature (K)
  Tveg                ! Vegetation temperature (K)

use LANDUSE, only : &
  fveg,              &! Canopy cover fraction
  tilefrac            ! Grid cell tile fraction

use TEST

implicit none

real*8, intent(in) :: &
  Ds1(Nx,Ny),        &! Surface layer thickness (m)
  KH(Nx,Ny),         &! Eddy diffusivity for heat to the atmosphere (m/s)
  KHa(Nx,Ny),        &! Eddy diffusivity from the canopy air space (m/s)
  KHv(Nx,Ny),        &! Eddy diffusivity for heat from vegetation (m/s)
  KWg(Nx,Ny),        &! Eddy diffusivity for water from the ground (m/s)
  KWv(Nx,Ny),        &! Eddy diffusivity for water from vegetation (m/s)
  ks1(Nx,Ny),        &! Surface layer thermal conductivity (W/m/K)
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny),      &! Net SW radiation absorbed by vegetation (W/m^2)
  Ts1(Nx,Ny)          ! Surface layer temperature (K)

real*8, intent(out) :: &
  Esrf(Nx,Ny),       &! Moisture flux from the surface (kg/m^2/s)
  Eveg(Nx,Ny),       &! Moisture flux from vegetation (kg/m^2/s)
  G(Nx,Ny),          &! Heat flux into surface (W/m^2)
  H(Nx,Ny),          &! Sensible heat flux to the atmosphere (W/m^2)
  Hsrf(Nx,Ny),       &! Sensible heat flux from the surface (W/m^2)
  LE(Nx,Ny),         &! Latent heat flux to the atmosphere (W/m^2)
  LEsrf(Nx,Ny),      &! Latent heat flux from the surface (W/m^2)
  LWsci(Nx,Ny),      &! Subcanopy incoming LWR (W/m^2)
  LWveg(Nx,Ny),      &! Subcanopy incoming LWR (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny),       &! Net radiation (W/m^2)
  Rsrf(Nx,Ny)         ! Net radiation absorbed by the surface (W/m^2)

integer :: & 
  i,j                 ! Point counters

real*8 :: &
  D,                 &! dQsat/dT (1/K)
  dE,                &! Change in surface moisture flux (kg/m^2/s)
  dG,                &! Change in surface heat flux (W/m^2)
  dH,                &! Change in sensible heat flux (W/m^2)
  dR,                &! Change in net radiation (W/m^2)
  dTs,               &! Change in surface skin temperatures (K)
  Lh,                &! Latent heat (J/kg)
  Qs,                &! Saturation humidity
  rho,               &! Air density (kg/m^3)
  Ssub                ! Mass of snow available for sublimation (kg/m^2)


! Initialize arrays for reproducability...

Esrf(:,:) = 0_dp
Eveg(:,:) = 0_dp
G(:,:) = 0_dp
H(:,:) = 0_dp
Hsrf(:,:) = 0_dp
LE(:,:) = 0_dp
LEsrf(:,:) = 0_dp
LWsci(:,:) = 0_dp
LWveg(:,:) = 0_dp
Melt(:,:) = 0_dp
Rnet(:,:) = 0_dp
Rsrf(:,:) = 0_dp


do j = 1, Ny
do i = 1, Nx

  if (tilefrac(i,j) < tthresh) goto 1 ! exclude points outside tile of interest

  if ((CANMOD == 1 .and. fveg(i,j) == 0) .or. CANMOD == 0) then
    Tveg(i,j) = Ta(i,j) 
    Tcan(i,j) = Ta(i,j) 

    ! Saturation humidity and density of air
    call QSAT(Ps(i,j),Tsrf(i,j),Qs)
    Lh = Lv
    if (Tsrf(i,j) < Tm .or. Sice(1,i,j) > epsilon(Sice(1,i,j))) Lh = Ls
    D = Lh*Qs/(Rwat*Tsrf(i,j)**2_dp)
    rho = Ps(i,j) / (Rair*Ta(i,j))

    ! Explicit fluxes
    Esrf(i,j) = rho*KWg(i,j)*(Qs - Qa(i,j))
    G(i,j) = 2_dp*ks1(i,j)*(Tsrf(i,j) - Ts1(i,j))/Ds1(i,j)
    H(i,j) = cp*rho*KH(i,j)*(Tsrf(i,j) - Ta(i,j))
    LE(i,j) = Lh*Esrf(i,j)
    Melt(i,j) = 0_dp
    Rnet(i,j) = SWsrf(i,j) + trcn(i,j)*LW(i,j) - sb*Tsrf(i,j)**4_dp  &
                            + (1_dp - trcn(i,j))*sb*Tveg(i,j)**4_dp

    ! Surface energy balance increments without melt
    dTs = (Rnet(i,j) - G(i,j) - H(i,j) - LE(i,j)) /  &
          (4_dp*sb*Tsrf(i,j)**3_dp + 2_dp*ks1(i,j)/Ds1(i,j) + rho*(cp*KH(i,j) + Lh*D*KWg(i,j)))
    dE = rho*KWg(i,j)*D*dTs
    dG = 2_dp*ks1(i,j)*dTs/Ds1(i,j) 
    dH = cp*rho*KH(i,j)*dTs
    dR = -4_dp*sb*Tsrf(i,j)**3_dp*dTs

    ! Surface melting
    if (Tsrf(i,j) + dTs > Tm .and. Sice(1,i,j) > epsilon(Sice(1,i,j))) then
      Melt(i,j) = sum(Sice(:,i,j))/dt
      dTs = (Rnet(i,j) - G(i,j) - H(i,j) - LE(i,j) - Lf*Melt(i,j)) /  &
            (4_dp*sb*Tsrf(i,j)**3_dp + 2_dp*ks1(i,j)/Ds1(i,j) + rho*(cp*KH(i,j) + Ls*D*KWg(i,j)))
      dE = rho*KWg(i,j)*D*dTs
      dG = 2_dp*ks1(i,j)*dTs/Ds1(i,j)
      dH = cp*rho*KH(i,j)*dTs
      dR = -4_dp*sb*Tsrf(i,j)**3_dp*dTs
      if (Tsrf(i,j) + dTs < Tm) then
        call QSAT(Ps(i,j),Tm,Qs)
        Esrf(i,j) = rho*KWg(i,j)*(Qs - Qa(i,j))  
        G(i,j) = 2_dp*ks1(i,j)*(Tm - Ts1(i,j))/Ds1(i,j)
        H(i,j) = cp*rho*KH(i,j)*(Tm - Ta(i,j))
        LE(i,j) = Ls*Esrf(i,j)
        Rnet(i,j) = SWsrf(i,j) + trcn(i,j)*LW(i,j) - sb*Tm**4_dp  &
                                + (1_dp - trcn(i,j))*sb*Tveg(i,j)**4_dp
        Melt(i,j) = (Rnet(i,j) - H(i,j) - LE(i,j) - G(i,j)) / Lf
        Melt(i,j) = max(Melt(i,j), 0.0_dp)
        dE = 0_dp
        dG = 0_dp
        dH = 0_dp
        dR = 0_dp
        dTs = Tm - Tsrf(i,j)
      end if
    end if
    
    ! In case of glacier without snow, cap Tsrf to 0°C
    ! This adjustment:
    !     - assumes the glacier is an infinite heat reservoir.
    !     - does not conserve energy.
    ! The excess energy would correspond to glacier melting, which we don't track.
    if (TILE == 'glacier') then
      if (Tsrf(i,j) + dTs > Tm .and. Sice(1,i,j) <= epsilon(Sice(1,i,j))) then
        call QSAT(Ps(i,j),Tm,Qs)
        Esrf(i,j) = rho*KWg(i,j)*(Qs - Qa(i,j))  
        G(i,j) = 2_dp*ks1(i,j)*(Tm - Ts1(i,j))/Ds1(i,j)
        H(i,j) = cp*rho*KH(i,j)*(Tm - Ta(i,j))
        LE(i,j) = Ls*Esrf(i,j)
        Rnet(i,j) = SWsrf(i,j) + trcn(i,j)*LW(i,j) - sb*Tm**4_dp  &
                                + (1_dp - trcn(i,j))*sb*Tveg(i,j)**4_dp
        dE = 0_dp
        dG = 0_dp
        dH = 0_dp
        dR = 0_dp
        dTs = Tm - Tsrf(i,j)
      end if
    end if

    ! Update surface temperature and fluxes
    Esrf(i,j) = Esrf(i,j) + dE
    G(i,j) = G(i,j) + dG
    H(i,j) = H(i,j) + dH
    LE(i,j) = Lh*Esrf(i,j)
    Rnet(i,j) = Rnet(i,j) + dR
    Tsrf(i,j) = Tsrf(i,j) + dTs

    ! Sublimation limited by amount of snow after melt
    Ssub = sum(Sice(:,i,j)) - Melt(i,j)*dt
    if (Ssub > epsilon(Ssub) .and. Esrf(i,j)*dt > Ssub) then
      Esrf(i,j) = Ssub / dt
      LE(i,j) = Ls*Esrf(i,j)
      H(i,j) = Rnet(i,j) - G(i,j) - LE(i,j) - Lf*Melt(i,j)
    end if
    Hsrf(i,j) = H(i,j)
    LEsrf(i,j) = LE(i,j)
    Rsrf(i,j) = Rnet(i,j)

    ! Ensure LWsci and LWveg exist as variable even in open runs
    LWsci(i,j) = LW(i,j)
    LWveg(i,j) = 0_dp

    if (CANMOD == 0) then
      ! Add fluxes from canopy in zero-layer model
      Eveg(i,j) = 0_dp
      if (fveg(i,j) > epsilon(fveg(i,j))) then
        Eveg(i,j) = - KWv(i,j)*Esrf(i,j) / (KHa(i,j) + KWv(i,j))
        H(i,j) = KHa(i,j)*H(i,j) / (KHa(i,j) + KHv(i,j))
        Lh = Ls
        if (Tveg(i,j) > Tm) Lh = Lv
        LE(i,j) = LE(i,j) + Lh*Eveg(i,j)
        Rnet(i,j) = Rnet(i,j) + SWveg(i,j) +  &
                    (1_dp - trcn(i,j))*(LW(i,j) + sb*Tsrf(i,j)**4_dp - 2_dp*sb*Tveg(i,j)**4_dp)
      end if
    end if
  end if

  1 continue
  
end do
end do

if (dump_data == 1) then
  open(1, file="temp/test_ebalsrf.txt")
  write(1,*) Esrf
  write(1,*) Eveg
  write(1,*) G
  write(1,*) H
  write(1,*) Hsrf
  write(1,*) LE
  write(1,*) LEsrf
  write(1,*) LWsci
  write(1,*) LWveg
  write(1,*) Melt
  write(1,*) Rnet
  write(1,*) Rsrf
  close(1) 
end if

end subroutine EBALSRF
