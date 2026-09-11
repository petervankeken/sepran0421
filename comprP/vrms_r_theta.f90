! compute rms of total velocity and of radial and polar components
! For cylindrical benchmark (Davies et al.)
! PvK May 2019
subroutine vrms_r_theta(isol,vrmsa)
use sepmodulecomio
use sepmodulekmesh
use sepmodulekprob
use sepran_arrays
use coeff
use geometry
use control
implicit none
integer :: isol
real(kind=8) :: vrmsa(*)
real(kind=8),dimension(:),allocatable :: userl
integer :: iuserl(100)

integer :: allocate_status,i,ihelp(5)
real(kind=8) :: volint


if (.not.allocated(userl)) allocate(userl(6+3*npoint),stat=allocate_status)
if (allocate_status /= 0) then
   if (print_node) write(irefwr,*) 'PERROR(vrms_r_theta): allocate_status = ',allocate_status
   call instop
endif
userl(1)=6+3*npoint

call vrms_rt01(npoint,ks(isol)%sol,nunkpi,nphysi,indprfi,kprobfi, &
     &                indprpi,kprobpi,coor,userl(6))

iuserl(1) = 100
iuserl(2) = 1
iuserl(3) = 0
iuserl(4) = 0
iuserl(5) = 0
iuserl(6) = 7
iuserl(7) = intrule900
iuserl(8) = icoor900
iuserl(9) = 0
iuserl(10) = 2001
iuserl(11) = 6

vrmsa(2)=sqrt(volint(0,1,1,kmesh1,kprob1,isol,iuserl,userl,ihelp)/volume)
iuserl(11) = 6+npoint
vrmsa(3)=sqrt(volint(0,1,1,kmesh1,kprob1,isol,iuserl,userl,ihelp)/volume)
iuserl(11) = 6+2*npoint
vrmsa(4)=volint(0,1,1,kmesh1,kprob1,isol,iuserl,userl,ihelp)/volume
do i=1,npoint
   userl(5+i)=sqrt(userl(5+i)*userl(5+i)+userl(5+i+npoint)*userl(5+i+npoint))
enddo
iuserl(11)=6
vrmsa(1)=sqrt(vrmsa(2)*vrmsa(2)+vrmsa(3)*vrmsa(3))

if (allocated(userl)) deallocate(userl)

end subroutine vrms_r_theta

subroutine vrms_rt01(npoint,usol,nunkp,nphys,indprf,kprobf,&
     &                      indprp,kprobp,coor,user)
use sepmodulecomio
use control
implicit none
integer :: npoint,nunkp,indprf,kprobf(*)
integer :: indprp,nphys,kprobp(npoint,nphys)
real(kind=8) :: usol(*),user(*),coor(2,*)
real(kind=8) :: x,y,r,theta,u,v,vr,vth
integer :: j1,j2,i
real(kind=8),parameter :: pi=3.14159265358979323846264338_8

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
    u = usol(j1)
    v = usol(j2)
    x=coor(1,i)
    y=coor(2,i)
    r=sqrt(x*x+y*y)
    theta=acos(y/r)
    if (x<=0) theta=2*pi-theta
    vth=cos(theta)*u - sin(theta)*v
    vr =sin(theta)*u + cos(theta)*v
    vth=vth
    user(i) = vth*vth
    user(i+npoint) = vr*vr
    user(i+2*npoint) = vth
!   write(irefwr,'(8f12.3)') x,y,u,v,r,theta,vr,vth
enddo

end subroutine vrms_rt01
