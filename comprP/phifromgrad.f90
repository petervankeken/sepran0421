! PHIFROMGRAD
!
! Compute viscous dissipation from gradients of velocity
! The viscous dissipation needs to be multiplied by viscosity
! outside this subroutine
! 
! Replacement of phifromsepran which may be inaccurate in case
! of compressible convection
!
! PvK 082708
subroutine phifromgrad(kmesh1,kprob1,isol1,iphi)
use control
use coeff
implicit none
integer :: kmesh1(*),kprob1(*),isol1(*),iphi(*)
integer :: igradv(5),idivv(5),istress(5),iuser(100),iinder(200)
real(kind=8) :: anorm,volint,penalty_parameter,pressure,user(100)
real(kind=8) :: gradv11,gradv12,gradv22,div11,stress11,stress12,stress22
integer :: ihelp,nparm,iwork(4,100),i
real(kind=8) :: work(100)

iuser(1)=100
 user(1)=100d0
iinder(1)=8
iinder(2)=1
iinder(3)=0
iinder(5)=0
iinder(6)=0
iinder(7)=0
iinder(8)=2
iuser(6)=7
!     ** itime, modelv, intrule, icoor, mcontv
iuser(2)=1
iuser(7)=0
iuser(8)=1
iuser(9)=intrule900
iuser(10)=icoor900
iuser(11)=mcontv
! eps, rho
iuser(12)=-6
if (compress) then
   iuser(13)=3
else
   iuser(13)=-7
endif
!      *** omega,f1,f2,f3,eta
iuser(14)=0
iuser(15)=0
iuser(16)=0
iuser(17)=0
iuser(18)=-7
if (isolmethod9==0) then
   user(6)=1d-6
else
   user(6)=0d0
endif
user(7)=1d0
! icheld = 2 for grad v
iinder(4)=2
call deriv(iinder,igradv,kmesh1,kprob1,isol1,iuser,user)

do i=1,100
   iwork(1,i)=0
   iwork(2,i)=0
   iwork(3,i)=0
   iwork(4,i)=0
    work(i) = 0d0
enddo
iwork(1,1) = 0
iwork(1,2) = 1
iwork(1,3) = intrule900
iwork(1,4) = icoor900
iwork(1,5) = mcontv
if (isolmethod9==0) then
  ! penaltyfunction parameter
  penalty_parameter = 1d-6
else
  penalty_parameter = 0
endif
! 6: epsilon
iwork(1,6) = 0
 work(6) = penalty_parameter
! 7: rho
if (compress) then
   iwork(1,7) = 3
else
   iwork(1,7) = 0
    work(7)   = 1d0
endif
! *** 12: eta assumed constant 1 here
iwork(1,12) = 0
 work(12)   = 1d0

do i=6,iuser(1)
   iuser(i)=0
enddo
do i=6,nint(user(1))
   user(i)=0d0
enddo
nparm=12
call fil103(1,1,iuser,user,kprob1,nparm,iwork,work,kmesh1)

call phifromgrad001(kmesh1,kprob1,isol1,igradv,iphi)

end subroutine phifromgrad

subroutine phifromgrad001(kmesh,kprob,isol,igradv,iphi)
use sepmodulekmesh
use sepmodulekprob
implicit none
integer :: kmesh,kprob,isol,igradv,iphi
integer :: iu1(2)
real(kind=8) :: u1(2)

integer :: ipcoor,ipdiv,ipgrad,ipphi,i
integer ndef_grad

call sepactsolbf1(isol)

if (iphi==0) then
   iu1(1)=0 ! problem #2
    u1(1)=0d0
   call creavc(0,1002,1,iphi,kmesh,kprob,iu1,u1,iu1,u1)
endif
ndef_grad=4 ! by default for 2D
call phifromgrad002(coor,npoint,ks(igradv)%sol,ndef_grad,ks(iphi)%sol,ks(isol)%sol,nunkp,kprobf,indprf)
end subroutine phifromgrad001

subroutine phifromgrad002(coor,npoint,grad,ngrad,phi,uvp,nunkp,kprobf,indprf)
use control
use coeff
implicit none
integer :: npoint,ndiv,ngrad,nphi,nunkp,indprf,kprobf(*)
real(kind=8) :: coor(2,npoint),grad(ngrad,npoint)
real(kind=8) :: phi(*)
real(kind=8) :: uvp(*)
integer :: i,j1
real(kind=8) :: stress11,stress12,stress21,stress22,divu

if (indprf.ne.0) then
  ! number of dofs per point is not constant
   do i=1,npoint
      ! location of upward velocity component
      j1 = kprobf(i)+2
      if (divufromeos) then
        divu = -DiG*uvp(j1)
        stress11 = 2*grad(1,i) - 2d0/3d0*divu
        stress22 = 2*grad(4,i) - 2d0/3d0*divu
        stress12 = grad(2,i)+grad(3,i)
      else
        stress11 = 4d0/3d0*grad(1,i) - 2d0/3d0*grad(4,i)
        stress22 = 4d0/3d0*grad(4,i) - 2d0/3d0*grad(1,i)
        stress12 = grad(2,i) + grad(3,i)
      endif
      phi(i) = stress11*grad(1,i)+stress22*grad(4,i) + stress12*stress12
  enddo
else
  ! number of dofs per point is constant 
  do i=1,npoint
      ! location of upward velocity component
      j1 = (i-1)*nunkp+2
      if (divufromeos) then
        divu = DiG*uvp(j1)
        stress11 = 2*grad(1,i) - 2d0/3d0*divu
        stress22 = 2*grad(4,i) - 2d0/3d0*divu   
        stress12 = grad(2,i)+grad(3,i)
      else
        stress11 = 4d0/3d0*grad(1,i) - 2d0/3d0*grad(4,i)
        stress22 = 4d0/3d0*grad(4,i) - 2d0/3d0*grad(1,i)
        stress12 = grad(2,i) + grad(3,i)
      endif
      phi(i) = stress11*grad(1,i)+stress22*grad(4,i) + stress12*stress12
  enddo
endif

end subroutine phifromgrad002





