! multiply viscous dissipation by viscosity
! icheld = 10   phi contains sqrt(edot::edot)
!          11   phi contains edot::edot
subroutine pemultphivisc(icheld)
use sepmodulevecs
use sepmodulekmesh
use sepran_arrays
implicit none
integer :: icheld
interface
   subroutine pemultiphivisc01(icheld,phi,eta,npoint)
   real(kind=8) :: phi(npoint),eta(npoint)
   integer :: npoint,icheld
   end subroutine pemultiphivisc01
end interface

call pemultiphivisc01(icheld,ks(iphi)%sol,ks(ivisc)%sol,npoint)

end subroutine pemultphivisc

subroutine pemultiphivisc01(icheld,phi,eta,npoint)
implicit none
real(kind=8) :: phi(npoint),eta(npoint)
integer :: npoint
integer :: icheld


if (icheld==10) then
   phi(1:npoint) = phi(1:npoint)*phi(1:npoint)*eta(1:npoint)
else if (icheld==11) then
   phi(1:npoint) = phi(1:npoint)*eta(1:npoint)
endif


end subroutine pemultiphivisc01
