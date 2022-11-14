!-----------------------------------------------------------------------
! Saturation specific humidity
!-----------------------------------------------------------------------
subroutine QSAT(P,T,Qs)

use CONSTANTS, only: &
  eps,               &! Ratio of molecular weights of water and dry air
  e0,                &! Saturation vapour pressure at Tm (Pa)
  Tm                  ! Melting point (K)

implicit none

real*8, intent(in) :: &
  P,                 &! Air pressure (Pa)
  T                   ! Temperature (K)

real*8, intent(out) :: &
  Qs                  ! Saturation specific humidity

real*8 :: &
  Tc,                &! Temperature (C)
  es                  ! Saturation vapour pressure (Pa)

Tc = T - Tm
if (Tc > 0) then
  es = e0*exp(17.5043*Tc / (241.3 + Tc))
else
  es = e0*exp(22.4422*Tc / (272.186 + Tc))
end if
Qs = eps*es / P

end subroutine QSAT
