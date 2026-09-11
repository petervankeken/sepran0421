subroutine find_C_int(ichoice,iuser,user,C_int)
use sepmodulekmesh
use sepmodulevecs
use sepmodulecomio
use sepran_arrays
use control
use coeff
use convparam
use geometry
implicit none
real(kind=8) :: user(*),C_int
integer :: iuser(*)
real(kind=8) :: volint,rinvec(2),x,y,funccf
integer :: ihelp,i,ichoice,iinvec(2),iprho
!     include 'pecof800.inc'
!     include 'ccc.inc'
!     include 'dimensional.inc'
!     include 'compression.inc'

iuser(2:iuser(1))=0
user(7:int(user(1)))=0.0_8
iuser(2)=1
iuser(6)=7
iuser(7)=intrule800+100*interpol800
iuser(8)=icoorsystem
iuser(9)=0
! specify f=rhobar
! (7) rho: store per point in user (even if constant)
iuser(10) = 2001
iprho=10
if (iprho+npoint-1>nint(user(1))) then
   if (print_node) then
      write(irefwr,*) 'PERROR(find_C_int): array user is too small at: ',nint(user(1))
      write(irefwr,*) '  we need it here to be at least: ',10+npoint-1
   endif
   call instop
endif
iuser(11) = iprho
if (compress.and.ratio_method) then
   ! define vector for reference density
   do i=1,npoint 
      x=coor(1,i)
      y=coor(2,i)
      user(iprho+i-1)=funccf(3,x,y,y)
   enddo
else
   ! if .not.ratiomethod density of particles defines rhobar here
   user(iprho:iprho+npoint-1)=1.0_8 
endif
if (print_node) write(irefwr,*) 'user(iprho: )',user(iprho)

C_int = volint(0,2,1,kmesh1,kprob1,idens,iuser,user,ihelp)


end subroutine find_C_int

subroutine fill_rho_C_int(npoint,coor,user)
implicit none
integer :: npoint,i
real(kind=8) :: coor(2,*),user(*)
real(kind=8) :: funccf,x,y


end subroutine fill_rho_C_int
