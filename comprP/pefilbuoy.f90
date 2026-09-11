! pefilbuoy
! Compute buoyancy force and store in user via fbuoy
! ic=1: f1
!    2: f2
! isol contains the temperature, iadia the adiabatic gradient, ipress the pressure
! and idens the chemical density
subroutine pefilbuoy(ic,fbuoy,isol,iadia,ipress,idens,ialpha)
use sepmodulekmesh
use sepmodulemesh
use sepmodulemeshinf
use sepmodulevecs
implicit none
integer,intent(in) :: ic,isol,iadia,ipress,idens,ialpha
real(kind=8),intent(inout) :: fbuoy(:)
interface
   subroutine pefilbuoy01(ic,ndim,npoint,coor,fbuoy,usol,adia,pressure,density,alpha)
     integer, intent(in) :: ic,ndim,npoint
     real(kind=8),intent(in) :: usol(:),adia(:),coor(:,:),pressure(:),density(:),alpha(:)
     real(kind=8),intent(inout) :: fbuoy(:)
   end subroutine pefilbuoy01
end interface

call pefilbuoy01(ic,ndim,npoint,coor,fbuoy,ks(isol)%sol,ks(iadia)%sol,ks(ipress)%sol,ks(idens)%sol,ks(ialpha)%sol)

end subroutine pefilbuoy

subroutine pefilbuoy01(ic,ndim,npoint,coor,fbuoy,usol,adia,pressure,density,alpha)
use sepmodulecomio
use geometry
use convparam
use coeff
use control
implicit none
integer, intent(in) :: ic,ndim,npoint
real(kind=8),intent(in) :: usol(:),adia(:),coor(:,:),pressure(:),density(:),alpha(:)
real(kind=8),intent(inout) :: fbuoy(:)
integer :: i,iph
real(kind=8) :: x,y,r,sint,cost,aload,funccf,alpha_here=1.0_8
real(kind=8) :: prespi,bigGamma,get_prespi
real(kind=8) :: temp_norm,prespi_norm,bigGamma_norm

if (cyl) then
   if (nph>0) then
      if (print_node) write(6,*) 'PERROR(pefilbuoy): not yet suited for phase changes and cyl'
      call instop
   endif
   do i=1,npoint
      x=coor(1,i)
      y=coor(2,i)
      r=sqrt(x*x+y*y)
      sint=x/r
      cost=y/r
      ! always fill alpha(:) even if constant 1
      alpha_here=alpha(i)
      if (compress.and. .not.solve_for_Tperturb) then
         aload=Ra*alpha_here*(usol(i)-adia(i))
      else
         aload=Ra*alpha_here*usol(i)
      endif
      !if (pvk_buoy_special) aload=Ra*usol(i)
      if (pvk_buoy_special) then
         if (print_node) write(irefwr,*) 'pvk_buoy_special ?!?'
         call instop
      endif
      if (fieldC) then
         if (stretch_tracers) then
            aload=aload-Rb_local*density(i)/funccf(3,x,y,y)
         else
            aload=aload-Rb_local*density(i)
         endif
      endif
      if (compress.and..not.TALA) then
         aload = aload - DiG*pressure(i)
      endif
      if (ic==1) then
         fbuoy(i)=aload*sint
      else if (ic==2) then
         fbuoy(i)=aload*cost
      endif
   enddo
else
   ! Cartesian
   open(906,file='fbuoy.dat')
   do i=1,npoint
     x=coor(1,i)
     y=coor(2,i)
     ! always use alpha(:)
     alpha_here=alpha(i)
     if (compress .and. .not.solve_for_Tperturb .and. .not.RaT) then
        fbuoy(i) = Ra*alpha_here*(usol(i)-adia(i))
     else 
        fbuoy(i) = Ra*alpha_here*usol(i)
     endif
     if (fieldC) then
        if (stretch_tracers) then
           ! correct for later multiplication by rhobar
           fbuoy(i)=fbuoy(i)-Rb_local*density(i)/funccf(3,x,y,y)
        else
           fbuoy(i)=fbuoy(i)-Rb_local*density(i)
        endif
     endif
     if (compress.and..not.TALA) then
        fbuoy(i) = fbuoy(i) - DiG*pressure(i)
     endif
     ! always add phase change topography buoyancy even for Di=0
     do iph=1,nph
        prespi=get_prespi(1-y,usol(i),iph)
        bigGamma=0.5+0.5*tanh(prespi/phdz(iph))
        ! make it optional to subtract Gamma_bar from buoyancy term for EBA
        ! needs update for ALA
        if (subtract_bigGamma_norm) then
           temp_norm=pht0(iph) 
           !prespi_norm=(1-y)-phz0(iph)-gamma(iph)*(temp_norm-pht0(iph))
           prespi_norm=get_prespi(1-y,temp_norm,iph)
           bigGamma_norm=0.5_8 + 0.5_8*tanh(prespi_norm/phdz(iph))
        else 
           bigGamma_norm=0.0_8
        endif
        fbuoy(i) = fbuoy(i) - Rb(iph)*(bigGamma-bigGamma_norm)
     enddo
     !if (x<1e-5) write(6,'(''900: '',10e12.5)') y,fbuoy(i),usol(i),bigGamma,bigGamma_norm,prespi,phdz(1),gamma(1),pht0(1),alpha_here
   enddo
endif
!    do i=1,npoint,npoint/10
!       write(6,'(''fbuoy: '',3e15.7)') usol(i),adia(i),fbuoy(i)
!    enddo


end subroutine pefilbuoy01

real(kind=8) function get_prespi(z,T,iph)
use coeff
implicit none
real(kind=8),intent(in) :: z,T
integer,intent(in) :: iph
! CY85 definition appears to say:
! get_prespi=phz0(iph)-z-gamma(iph)*(T-pht0(iph))
! but in a subtle twist they define z as the vertical coordinate, not depth, so:
get_prespi=z-phz0(iph)-gamma(iph)*(T-pht0(iph))
end function get_prespi

real(kind=8) function get_prespibar(z,iph)
use coeff
implicit none
real(kind=8),intent(in) :: z
integer,intent(in) :: iph
! CY85 definition appears to say:
! get_prespi=phz0(iph)-z-gamma(iph)*(T-pht0(iph))
! but in a subtle twist they define z as the vertical coordinate, not depth, so:
get_prespibar=z-phz0(iph)
end function get_prespibar
