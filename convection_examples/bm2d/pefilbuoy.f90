! PvK July 2024
! pefilbuoy
! Compute buoyancy force and store in user via fbuoy
! ic=1: f1
!    2: f2
! isol contains the temperature, iadia the adiabatic gradient, ipress the pressure
! and idens the chemical density
subroutine pefilbuoy(ic,fbuoy,isol) !,iadia,ipress,idens,ialpha,irho)
use sepmodulekmesh
use sepmodulemesh
use sepmodulemeshinf
use sepmodulevecs
implicit none
integer,intent(in) :: ic,isol ! ,iadia,ipress,idens,ialpha,irho
real(kind=8),intent(inout) :: fbuoy(:)
interface
   subroutine pefilbuoy01(ic,ndim,npoint,coor,fbuoy,usol) ! ,adia,pressure,density,alpha,rho)
     integer, intent(in) :: ic,ndim,npoint
     real(kind=8),intent(in) :: usol(:) ! ,adia(:),coor(:,:),pressure(:),density(:),alpha(:),rho(:)
     real(kind=8),intent(in) :: coor(:,:)
     real(kind=8),intent(inout) :: fbuoy(:)
   end subroutine pefilbuoy01
end interface

call sepactsolbf1(isol)
call pefilbuoy01(ic,ndim,npoint,coor,fbuoy,ks(isol)%sol) ! ,ks(iadia)%sol,ks(ipress)%sol,ks(idens)%sol,ks(ialpha)%sol,ks(irho)%sol)

end subroutine pefilbuoy

subroutine pefilbuoy01(ic,ndim,npoint,coor,fbuoy,usol) ! ,adia,pressure,density,alpha,rho)
use sepmodulecomio
use geometry
use coeff
use control
implicit none
integer, intent(in) :: ic,ndim,npoint
real(kind=8), intent(in) :: usol(:),coor(:,:)
real(kind=8), intent(inout) :: fbuoy(:)
integer :: i
real(kind=8) :: x,y

if (cyl) then
   if (print_node) then
      write(irefwr,*) 'PERROR(pefilbuoy01): not yet suited for cyl'
      call instop
   endif
else
  do i=1,npoint
     x=coor(1,i)
     y=coor(2,i)
     fbuoy(i)=Ra*usol(i)
  enddo
endif

end subroutine pefilbuoy01
