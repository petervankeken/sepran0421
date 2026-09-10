real(kind=8) function get_prespi(z,T,iph)
use coeff
implicit none
real(kind=8),intent(in) :: z,T
integer,intent(in) :: iph
! CY85 definition appears to say:
! get_prespi=phz0(iph)-z-gamma(iph)*(T-pht0(iph))
! but in a subtle twist they define z as the vertical coordinate, not depth, so:
get_prespi=z-phz0(iph)-gamma(iph)*(T-pht0(iph))
end function get_prespi

real(kind=8) function get_prespibar(z,iph)
use coeff
implicit none
real(kind=8),intent(in) :: z
integer,intent(in) :: iph
! CY85 definition appears to say:
! get_prespi=phz0(iph)-z-gamma(iph)*(T-pht0(iph))
! but in a subtle twist they define z as the vertical coordinate, not depth, so:
get_prespibar=z-phz0(iph)
end function get_prespibar
