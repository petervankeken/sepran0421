! PvK July 2024
! Start from scratch on a bm2d version with an eye out on comprP 
subroutine fillcoef900(ichoice)
use sepmodulekmesh
use sepran_arrays
use control
use coeff
use geometry
implicit none
integer,intent(in) :: ichoice
!  call coef900_print(npoint,coor,user(ihorvel),user(ivervel),user(ibeta),user(iright),user(irhocp),user(itemp),user(iuser_cond))
interface
  subroutine coef900_print(npoint,coor,rho,f1,f2,eta)
     integer,intent(in) :: npoint
     real(kind=8),intent(in),dimension(*) :: rho,f1,f2,eta
     real(kind=8),intent(in),dimension(2,*) :: coor
  end subroutine coef900_print
  subroutine pefilbuoy(ic,fbuoy,isol)
     integer,intent(in) :: ic,isol
     real(kind=8),intent(inout) :: fbuoy(:)
  end subroutine pefilbuoy
end interface
integer :: iprho,ipbuoy1,ipbuoy2,ipeta,ip,ipuser,i,ic
logical :: first=.true.
save first



! zero out arrays
iuser_here(2:iuser_here(1))=0
 user_here(2:nint(user_here(1)))=0.0_8

iuser_here(2)=1 ! make sure sepran knows there is only one element group
iuser_here(6)=8 ! start information about element group at iuser_here position 8
ipuser=10

! define integer information for coefficients
ip=7
iuser_here(ip+1) = 0 ! boring old stokes
iuser_here(ip+2) = 1 ! type of constitutive equation
iuser_here(ip+3) = intrule900+100*interpol900  ! let element routines define integration 
iuser_here(ip+4) = icoor900    ! Cartesian
iuser_here(ip+5) = mcontv ! 0=incompressible
!write(6,*) 'mcontv: ',mcontv

! (6) eps: penalty function parameter
iuser_here(ip+6) = -7
 user_here(7)    = 0.0_8
if (itype_stokes==900) user_here(7)=penalty_parameter
!if (pedebug) write(irefwr,*) 'penalty function parameter: ',user_here(7)

! (7) rho: store per point in user (even if constant)
iuser_here(ip+7) = 2001
iuser_here(ip+8) = ipuser
iprho = ipuser
if (compress) then
   ! define vector for reference density
   ! call filluser900(npoint,ks(irho)%sol,user_here(iprho))
   if (print_node) write(irefwr,*) 'PERROR(fillcoef900): not yet suited for compress'
   call instop
else
   user_here(iprho:iprho+npoint-1)=1.0_8 ! constant rho for now
endif
ipuser=ipuser+npoint
ip=ip+1 ! to account for offset due to use of 2001

! (8) omega: constant 0
iuser_here(ip+8) = 0

ipbuoy1=1
ipbuoy2=1
! (10) f1: look up per nodal point in user
iuser_here(ip+9) = 2001
iuser_here(ip+10) = ipuser
ipbuoy1=ipuser
if (cyl.and.ichoice==1) then
   ! take into account thermal buoyancy
   call sepactsolbf1(isol(2))
   call sepgetmeshinfo(ndim,npoint)
   ic=1
   call pefilbuoy(ic,user_here(ipbuoy1:ipbuoy1+npoint-1),isol(2)) ! ,iadia,ipress,idens,ialpha,irho)
else
   user_here(ipbuoy1:ipbuoy1+npoint-1)=0.0_8
endif
ipuser=ipuser+npoint
ip=ip+1 ! account for increment since type 2001 takes two integer spaces

! (10) f2: look up per nodal point in user
iuser_here(ip+10) = 2001
iuser_here(ip+11) = ipuser
ipbuoy2=ipuser
if (ichoice==1) then ! .and..not.T_buoyancy_through_particles) then
   call sepactsolbf1(isol(2))
   ic=2
   call pefilbuoy(ic,user_here(ipbuoy2:ipbuoy2+npoint-1),isol(2)) ! ,iadia,ipress,idens,ialpha,irho)
else
   ! no buoyancy for Gable plates or when temperature buoyancy is carried by particles
   user_here(ipbuoy2:ipbuoy2+npoint-1)=0.0_8
endif

ipuser=ipuser+npoint
ip=ip+1 ! account for increment since type 2001 takes two integer spaces

! (11) f3: constant 0
iuser_here(ip+11) = 0

! (12) eta. Always store in user when variable & mostly continuous...
iuser_here(ip+12) = 2001
iuser_here(ip+13) = ipuser
ipeta=ipuser
call coefvis()
if (itypv==0) then
   user_here(ipeta:ipeta+npoint-1)=1.0_8
else 
   call sepactsolbf1(ivisc)
   call pecopy(0,user_here(ipeta:ipeta+npoint-1),ivisc)
endif
ipuser = ipuser+npoint
ip=ip+1

if (first.and.print_node) then
   if (.not.petest) first=.false.
   write(irefwr,*) 'iuser/user: '
   write(irefwr,*) 'ichoice,iprho,ipbuoy1,ipbuoy2,ipeta=',ichoice,iprho,ipbuoy1,ipbuoy2,ipeta
   do i=1,100
      write(irefwr,'(i20,i20,e15.7)') i,iuser_here(i),user_here(i)
   enddo
   !write(6,*) 'ihorvel etc.: ',ihorvel,ivervel,ibeta,iright,irhocp,itemp,iuser_cond,ialpha,ivisdip
   call coef900_print(npoint,coor,user_here(iprho),user_here(ipbuoy1),user_here(ipbuoy2),user_here(ipeta))
   first=.false.
endif


end subroutine fillcoef900

subroutine coef900_print(npoint,coor,rho,f1,f2,eta)
use sepmodulemesh
use sepmodulecomio
use control
implicit none
integer,intent(in) :: npoint
real(kind=8),dimension(*),intent(in) :: rho,f1,f2,eta
real(kind=8),dimension(2,*),intent(in) :: coor
real(kind=8) :: r
integer :: i

write(irefwr,'(''coef900: '',11a12,: )') 'x','y','r','rho','f1','f2','eta'
do i=1,npoint,npoint/10
   r=sqrt(coor(1,i)*coor(1,i)+coor(2,i)*coor(2,i))
   write(irefwr,'(''         '',11f12.3,:)') coor(1,i),coor(2,i),r,rho(i),f1(i),f2(i),eta(i)
enddo

end subroutine coef900_print
