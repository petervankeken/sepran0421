! PvK July 2024
! Start from scratch on a bm2d version with an eye out on comprP 
subroutine fillcoef800()
use sepmodulekmesh
use sepran_arrays
use control
implicit none
!  call coef800_print(npoint,coor,user(ihorvel),user(ivervel),user(ibeta),user(iright),user(irhocp),user(itemp),user(iuser_cond))
interface
  subroutine coef800_print(npoint,coor,horvel,vervel,beta,right,rhocp,temp,cond)
     integer,intent(in) :: npoint
     real(kind=8),intent(in),dimension(*) :: horvel,vervel,beta,right,rhocp,temp,cond
     real(kind=8),intent(in),dimension(2,*) :: coor
  end subroutine coef800_print
end interface
integer :: iuser_cond,ihorvel,ivervel,ibeta,iright,irhocp,itemp,ip,i
logical :: first=.true.
integer :: intrule800,interpol800,icoor800
save first

! debug - works to give a conductive solution
!iuser_here(2)=1
!iuser_here(6)=7
!iuser_here(7)=0
!iuser_here(8)=0 ! metupw
!iuser_here(9)=0 ! intrule
!iuser_here(10)=0 ! icoorsystem
!iuser_here(11)=0 ! not yet used
!iuser_here(12)=-6 ! k11
!iuser_here(15)=-6 ! k22
!user_here(6)=1.0_8
!return

! pointers to information in user_here
iuser_cond=10
ihorvel=10+npoint
ivervel=10+2*npoint
ibeta=10+3*npoint
iright=10+4*npoint
irhocp=10+5*npoint
itemp=10+6*npoint

! zero out arrays
iuser_here(2:iuser_here(1))=0
 user_here(2:nint(user_here(1)))=0.0_8

iuser_here(2)=1 ! make sure sepran knows there is only one element group
iuser_here(6)=15 ! start information about element group at iuser_here position 15

! define integer information for coefficients
ip=14
! (1) not used
iuser_here(ip+1) = 0
! (2) type of upwinding ; avoid with P2 for now
iuser_here(ip+2) = metupw
! (3) intrule
iuser_here(ip+3) = intrule800+100*interpol800  ! let element routines define integration 
! (4) icoor: type of coordinate system
iuser_here(ip+4) = icoor800    ! Cartesian
! (5) not yet used
iuser_here(ip+5) = 0

! (6) k11
iuser_here(ip+6) = 2001
iuser_here(ip+7) = iuser_cond
ip=ip+1 ! account for increment since type 2001 takes two integer spaces
! (7) k12
iuser_here(ip+7) = 0
! (8) k13
iuser_here(ip+8) = 0
! (9) k22
iuser_here(ip+9) = 2001
iuser_here(ip+10) = iuser_cond
ip=ip+1
! (10) k23
iuser_here(ip+10) = 0
! (11) k33
iuser_here(ip+11) = 0

conductive=.false.
! velocity components
if (conductive) then
  iuser_here(ip+12) = 0  ! u=0
  iuser_here(ip+13) = 0  ! v=0
else
  ! (12) u
  iuser_here(ip+12) = 2001
  iuser_here(ip+13) = ihorvel
  ip=ip+1 ! account for increment since type 2001 takes two integer spaces
  ! (13) v
  iuser_here(ip+13) = 2001
  iuser_here(ip+14) = ivervel
  ip=ip+1 ! account for increment since type 2001 takes two integer spaces
endif

! (14) w==0
iuser_here(ip+14) = 0

! (15) beta
iuser_here(ip+15) = 2001
iuser_here(ip+16) = ibeta
ip=ip+1 ! account for increment since type 2001 takes two integer spaces

! (16) q    Note we use iright because we may add viscous dissipation to form
!           the right-hand side
iuser_here(ip+16)=2001
iuser_here(ip+17)=iright
ip=ip+1 ! account for increment since type 2001 takes two integer spaces

! (17) rho*cp
iuser_here(ip+17) = 2001
iuser_here(ip+18) = irhocp
ip=ip+1


! start filling point wise values in user

if (.not.conductive) then
   call pecopy(2,user_here(ihorvel:ihorvel+npoint-1),isol)
   call pecopy(3,user_here(ivervel:ivervel+npoint-1),isol)
else
   user_here(ihorvel:ihorvel+npoint-1)=0.0_8 
   user_here(ivervel:ivervel+npoint-1)=0.0_8 
endif

! assume conductivity and rho*cp are constant 
user_here(iuser_cond:iuser_cond+npoint-1)=1.0_8
user_here(irhocp:irhocp+npoint-1)=1.0_8

if (first.and.print_node) then
   if (.not.petest) first=.false.
   do i=1,100
      write(irefwr,'(i20,i20,e15.7)') i,iuser_here(i),user_here(i)
   enddo
   !write(6,*) 'ihorvel etc.: ',ihorvel,ivervel,ibeta,iright,irhocp,itemp,iuser_cond,ialpha,ivisdip
   call coef800_print(npoint,coor,user_here(ihorvel),user_here(ivervel),user_here(ibeta),user_here(iright), & 
        & user_here(irhocp),user_here(itemp),user_here(iuser_cond))
!        & ks(ialpha)%sol,ks(ivisdip)%sol)
   first=.false.
endif


end subroutine fillcoef800

subroutine coef800_print(npoint,coor,u,v,beta,right,rhocp,temp,cond) !,alpha,phi)
use sepmodulecomio
implicit none
integer,intent(in) :: npoint
real(kind=8),dimension(*),intent(in) :: u,v,beta,right,rhocp,temp,cond ! ,alpha,phi
real(kind=8),dimension(2,*),intent(in) :: coor
integer :: i

write(irefwr,'(''coef800: '',11a12 )') 'x','y','u','v','beta','right','rhocp','temp','cond' ! ,'alpha','phi'
do i=1,npoint,npoint/10
   write(irefwr,'(''         '',11f12.3)') coor(1,i),coor(2,i),u(i),v(i),beta(i),right(i),rhocp(i),temp(i),cond(i)!,alpha(i),phi(i)
enddo

end subroutine coef800_print

