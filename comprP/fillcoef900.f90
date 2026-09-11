! fillcoef900
!
! Fill coefficients for Stokes equations directly in iuser/user
! 
! Updated February 2020 to get ready for compressible convection with variable rho
!
! Called from pefilcof_cart that prepares viscosity vector
! as well as pressure vector for 903
! For 900 it is more indirect as we need to fill coefficients to be able to use deriv to find pressure
! so fillcoef900 should be called twice.
subroutine fillcoef900(ichoice,iuserh,userh)
use sepmodulecomio
use sepmodulekmesh
use sepmodulevecs
use geometry
use control
use coeff
use sepran_interface
use sepran_arrays
implicit none
!interface
!  subroutine pecopy(ichoice,user,isol)
!     integer,intent(in) :: ichoice,isol
!     real(kind=8),intent(inout) :: user(:)
!  end subroutine pecopy
!end interface
real(kind=8) :: userh(*)
integer :: ichoice,iuserh(*)
integer :: ipuser,ip,ic,i
integer :: ipbuoy1,ipbuoy2,ipeta,iprho
real(kind=8) :: f1,f2,ftot
logical :: first=.true.


!if (ichoice==4) then
!   ! test why plate solutions are weird
    ! answer: something doesn't seem right when assuming  jchois=iinbld(2)=5
!   call pefilcof_cyl_old(ichoice,iuserh,userh)
!   return
!endif

!call sepgetmeshinfo(ndim,npoint) ! now from sepmodulekmesh


iuserh(3:iuserh(1))=0
 userh(2:nint(userh(1)))=0.0_8

! Define coefficients through iuser/user (bit convoluted but robust)
iuserh(6)=8 ! start information about first and only element group at iuser position 8
ipuser=10  ! pointer to nodal point information. 

! first define integer information: itime, modelv, intrule, icoor, mcontv
! (1) itime: type of Navier-Stokes equation
ip=7
iuserh(ip+1) = 0  ! boring old Stokes
! (2) modelv: type of constitutive equation, specify here through coefficient #12
iuserh(ip+2) = 1
! (3) intrule: type of element integration routines
iuserh(ip+3) = intrule900 + 100*interpol900  ! let element routines define integration but linearly interpolate coefficients to Gauss points if necessary
! (4) icoor: type of coordinate system
iuserh(ip+4) = 0    ! Cartesian
! (5) mcontv: compressibility formulation
iuserh(ip+5) = mcontv    ! 0=incompressible
!write(irefwr,*) 'mcontv: ',mcontv

!pedebug=.true.

! (6) eps: penalty function parameter
iuserh(ip+6) = -7
 userh(7)    = 0.0_8   
if (itype_stokes==900) userh(7)=1.0e-6_8
!if (pedebug) write(irefwr,*) 'penalty function parameter: ',userh(7)

! (7) rho: store per point in user (even if constant)
iuserh(ip+7) = 2001
iuserh(ip+8) = ipuser
iprho = ipuser
if (compress) then
   ! define vector for reference density
   call filluser900(npoint,ks(irho)%sol,userh(iprho))
else
   userh(iprho:iprho+npoint-1)=1.0_8 ! constant rho for now
endif
ipuser=ipuser+npoint
ip=ip+1 ! to account for offset due to use of 2001

! (8) omega: constant 0
iuserh(ip+8) = 0

ipbuoy1=1
ipbuoy2=1
! (10) f1: look up per nodal point in user
iuserh(ip+9) = 2001
iuserh(ip+10) = ipuser
ipbuoy1=ipuser
if (cyl.and.ichoice==1) then
   ! take into account thermal buoyancy
   call sepactsolbf1(isol2)
   call sepgetmeshinfo(ndim,npoint)
   ic=1
   call pefilbuoy(ic,userh(ipbuoy1:ipbuoy1+npoint-1),isol2,iadia,ipress,idens,ialpha) 
else
   userh(ipbuoy1:ipbuoy1+npoint-1)=0.0_8
endif
ipuser=ipuser+npoint
ip=ip+1 ! account for increment since type 2001 takes two integer spaces

! (10) f2: look up per nodal point in user
iuserh(ip+10) = 2001
iuserh(ip+11) = ipuser
ipbuoy2=ipuser
if (ichoice==1.and..not.T_buoyancy_through_particles) then
   call sepactsolbf1(isol2)
   ic=2
   call pefilbuoy(ic,userh(ipbuoy2:ipbuoy2+npoint-1),isol2,iadia,ipress,idens,ialpha)
else
   ! no buoyancy for Gable plates or when temperature buoyancy is carried by particles
   userh(ipbuoy2:ipbuoy2+npoint-1)=0.0_8
endif

ipuser=ipuser+npoint
ip=ip+1 ! account for increment since type 2001 takes two integer spaces

! (11) f3: constant 0
iuserh(ip+11) = 0

! 12: eta. Always store in user when variable
iuserh(ip+12) = 2001
iuserh(ip+13) = ipuser
ipeta=ipuser
if (itypv==0) then
  userh(ipeta:ipeta+npoint-1)=1.0_8
else
  call sepactsolbf1(ivisc)
  call pecopy(0,userh(ipeta:ipeta+npoint-1),ivisc)
endif
ipuser = ipuser+npoint
ip=ip+1


if (first.and.print_node) then
   if (.not.petest) first=.false.
   write(irefwr,*) 'iuser/user: '
   write(irefwr,*) 'ichoice,iprho,ipbuoy1,ipbuoy2,ipeta= ',ichoice,iprho,ipbuoy1,ipbuoy2,ipeta
   do i=1,100
      write(irefwr,'(i20,i20,e15.7)') i,iuserh(i),userh(i)
   enddo
   !write(irefwr,*) 'density at top = ',userh(iprho+npoint-1)
   call coef900_print(npoint,coor,userh(iprho),userh(ipbuoy1),userh(ipbuoy2),userh(ipeta))
   !do i=200,6*npoint,100
   !   write(irefwr,'(i20,20x,e15.7)') i,user(i)
   !enddo
endif

pedebug=.false.

end subroutine fillcoef900

subroutine filluser900(npoint,rho,urho)
implicit none
integer :: npoint
real(kind=8),dimension(*) :: rho,urho

urho(1:npoint)=rho(1:npoint)

end subroutine filluser900
