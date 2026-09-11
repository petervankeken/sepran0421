subroutine subtract_average_rotation(cally)
use sepmodulecomio
use sepran_arrays
use geometry
use control
use mtime
implicit none
integer :: iinder(7),ihelp
real(kind=8) :: rotation,volint,rotation2
real(kind=4) :: t11,t12
character(len=80) :: namedof(3),cally

iinder(1)=7
iinder(2)=1  ! ichoice: output is vector of special structure
iinder(3)=0  ! real vector
iinder(4)=5  ! icheld: compute curl v
iinder(5)=1  ! ix (default)
iinder(6)=1  ! jdegfd (default)
iinder(7)=12 ! ivec for curl v
call deriv(iinder,icurl,kmesh1,kprob1,isol1,iuser_here,user_here)

! Put curl v in user for integration
call prep_volint_rotation(icurl,iuser_here,user_here)
! actual integration
rotation=volint(0,1,1,kmesh1,kprob1,isol1,iuser_here,user_here,ihelp)/volume

! point everything to isol1
call sepactsolbf1(isol1)
! subtract rotation from isol1
call subtract_ar01(isol1,rotation)

! recompute rotation
call deriv(iinder,icurl,kmesh1,kprob1,isol1,iuser_here,user_here)
call prep_volint_rotation(icurl,iuser_here,user_here)
! actual integration
rotation2=volint(0,1,1,kmesh1,kprob1,isol1,iuser_here,user_here,ihelp)/volume

if (petest.and.print_node) write(irefwr,'(''PINFO(subtract_average_rotation) is called from '',a80)') cally
!if (print_node) write(irefwr,'(''PINFO(subtract_average_rotation): rotation was reduced from '',e15.7,'' to '',e15.7)') rotation,rotation2
if (print_node) then
   write(LU_ROTATION,'(3e15.7)') time_now,rotation,rotation2
   call flush(LU_ROTATION)
endif
!if (print_node) then
!  ourplotname='veloc.0002'
!  call plotvc(1,2,isol1,isol1,kmesh1,kprob1,15d0,1d0,0.1d0)
!  namedof(1)='velocity'
!  write(ourplotname,'(''UV'',i4.4)') 2
!  call sol2vtu(2,isol1,ourplotname,0,namedof)
!  if (petest) call instop
!endif

end subroutine subtract_average_rotation

subroutine subtract_ar01(isol1,rotation)
use sepmodulekmesh
use sepmodulekprob
use sepmodulevecs
use coeff
use geometry
implicit none
integer :: isol1
real(kind=8) :: rotation

call subtract_ar02(ks(isol1)%sol,coor,indprfi,kprobfi,indprpi,kprobpi,npoint,nphys,nunkp,rotation)

end subroutine subtract_ar01

! subtract rigid body rotation around z-axis from usol.
! rotation contains the average curl over the domain
! PvK 082814
subroutine subtract_ar02(usol,coor,indprf,kprobf,indprp,kprobp,npoint,nphys,nunkp,rotation)
implicit none
real(kind=8) :: usol(*),coor(2,*),rotation
integer :: npoint,indprf,kprobf(*),indprp
integer :: nphys,nunkp,kprobp(npoint,nphys)

integer :: i,j,j1,j2
real(kind=8) :: u,v,vr,vth,r,theta,x,y,costh,sinth,fac
real(kind=8),parameter ::  pi=3.1415926535898_8

fac = rotation/2d0
if (indprf == 0 .and. indprp == 0) then
   do i=1,npoint
      u = usol(2*i-1)
      v = usol(2*i)
      x = coor(1,i)
      y = coor(2,i)
      r = sqrt(x*x+y*y)
      if (x>0) then
         theta = acos(y/r)
      else
         theta = 2*pi-acos(y/r)
      endif
      u = u + r*fac*cos(theta)
      v = v - r*fac*sin(theta)
      usol(2*i-1) = u
      usol(2*i) = v
   enddo
else if (indprp /= 0) then
   do i=1,npoint
      j1 = kprobp(i,1)
      j2 = kprobp(i,2)
      u = usol(j1)
      v = usol(j2)
      x = coor(1,i)
      y = coor(2,i)
      r = sqrt(x*x+y*y)
      if (x>0) then
         theta = acos(y/r)
      else
         theta = 2*pi-acos(y/r)
      endif

      u = u + r*fac*cos(theta)
      v = v - r*fac*sin(theta)
      usol(j1) = u
      usol(j2) = v
!     write(irefwr,'(i5,8f12.3)') i,x,y,u,v,r,theta,vr,vth
    enddo
else if (indprf /= 0) then
    do i=1,npoint
       j = kprobf(i)
       u = usol(2*j-1)
       v = usol(2*j)
       x = coor(1,i)
       y = coor(2,i)
       r = sqrt(x*x+y*y)
       if (x>0) then
          theta = acos(y/r)
       else
          theta = 2*pi-acos(y/r)
       endif
       u = u + r*fac*cos(theta)
       v = v - r*fac*sin(theta)
       usol(2*j-1) = u
       usol(2*j) = v
    enddo
endif

end subroutine subtract_ar02

subroutine prep_volint_rotation(icurl,iuser_here,user_here)
use sepmodulekmesh
use sepmodulekprob
use sepmodulevecs
use coeff
use geometry
implicit none
integer :: icurl,iuser_here(*)
real(kind=8) :: user_here(*)

call sepactsolbf1(icurl)
call prep_vr02(ks(icurl)%sol,user_here(6))
iuser_here(1) = 100
iuser_here(2) = 1
iuser_here(3) = 0
iuser_here(4) = 0
iuser_here(5) = 0
iuser_here(6) = 7
iuser_here(7) = intrule900
iuser_here(8) = icoor900
iuser_here(9) = 0
iuser_here(10) = 2001
iuser_here(11) = 6

end subroutine prep_volint_rotation

subroutine prep_vr02(curl,user_here)
use sepmodulekmesh
use sepmodulekprob
use sepmodulevecs
use coeff
use geometry
implicit none
real(kind=8),intent(in) :: curl(*)
real(kind=8),intent(out) :: user_here(*)

user_here(1:npoint)=curl(1:npoint)

end subroutine prep_vr02
