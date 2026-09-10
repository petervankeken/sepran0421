! PvK July 2024
!   PEVRMS
!
!   Calculates u*u + v*v per nodal point from the
!   solution vector and places the result in array user.
!
!   For element 900/903 with fixed number of degrees of 
!   freedom (e.g., 900 with ISHAPE=6, 903 with ISHAPE=3)
!   or for variable (e.g., 903 with ISHAPE=4)
!
!   PvK, 17-4-89/930107/970105/
!   PvK 190700
!   Updated for 0619 PvK 070119
subroutine pevrms(vrms)
!!use sepmodulemain
use sepran_arrays
use sepmodulekmesh
use coeff
use geometry
implicit none
real(kind=8),intent(inout) :: vrms(2)
real(kind=8) :: volint,outval(3),vrms2
integer :: ihelp,iinvol(5)

call sepactsolbf1(isol(1))
call pevrms00(isol(1),iuser_here,user_here)
! standard volume integral
iinvol=0
call integr(iinvol,outval,kmesh,kprob,isol,iuser_here,user_here)
vrms(1) = sqrt(abs(outval(1))/volume)

! alternative test for comparison with Scott November 2021
vrms(2) = sqrt(sum(user_here(6:6+npoint-1))/npoint)

end subroutine pevrms

subroutine pevrms00(isol,iuser_here,user_here)
!use sepmodulemain
use sepmodulekmesh
use sepmodulekprob
use sepmodulevecs
use coeff
use geometry
implicit none
integer,intent(in) :: isol
integer,intent(inout) :: iuser_here(*)
real(kind=8),intent(inout) :: user_here(*)
integer :: ihelp



call pevrms01(npoint,ks(isol)%sol,nunkpi,nphysi,indprfi,kprobfi,indprpi,kprobpi,user_here(6))

! correct length has already been set
!iuser_here(1) = 100
!iuser_here(2) = 1
iuser_here(3) = 0
iuser_here(4) = 0
iuser_here(5) = 0
iuser_here(6) = 7
iuser_here(7) = intrule900
iuser_here(8) = icoor900
iuser_here(9) = 0
iuser_here(10) = 2001
iuser_here(11) = 6


end subroutine pevrms00

subroutine pevrms01(npoint,usol,nunkp,nphys,indprf,kprobf,indprp,kprobp,user_here) 
implicit none
integer :: npoint,nunkp,indprf,kprobf(*)
integer :: indprp,nphys,kprobp(npoint,nphys)
real(kind=8) :: usol(*),user_here(*)
integer :: j1,j2,i

do i=1,npoint
   if (indprf == 0 .and. indprp == 0) then
      j1 = (i-1)*nunkp+1
      j2 = (i-1)*nunkp+2
   else if (indprp /= 0) then
      j1 = kprobp(i,1)
      j2 = kprobp(i,2)
   else
      j1 = kprobf(i) + 1
      j2 = j1 + 1
   endif
   user_here(i) = usol(j1)*usol(j1)+usol(j2)*usol(j2)
enddo
 
end subroutine pevrms01
