! PECOPHI
! copy viscous dissipation (in iphi) and radiogenic heating (in iheat) to user
! PvK 12 August 2020 (rewritten from 1997)
subroutine pecophi(user,iphi,iheat)
use sepmodulekmesh
use sepmodulekprob
use coeff
implicit none
integer :: iphi,iheat
real(kind=8) :: user(*)

call pecophi01(user,ks(iphi)%sol,ks(iheat)%sol,npoint)

end subroutine pecophi

subroutine pecophi01(user,phi,internal_heating,npoint)
use coeff
implicit none
integer :: npoint
real(kind=8) :: user(*),phi(*),internal_heating(*)

user(1:npoint)=0.0_8

if (iqtype>0) then
   ! multiplication with rho should have been done in findheatgen()
   user(1:npoint)=internal_heating(1:npoint)
endif
if (DiRa>0) then
   user(1:npoint)=user(1:npoint)+DiRa*phi(1:npoint)
endif

end subroutine pecophi01
