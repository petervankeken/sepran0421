! update_rho_adia
! Set up an EoS based on a lookup table, but modify density by integrating the density gradient
! the respective gradients. Take into account phase change influence on density
! See ~/sepran/Convection/theory 
!
! The effects of the phase changes introduces a modified Rb: Rb1=Rb*rhor*rhor/(rho2*rho1)
! where rho2,rho1 are the densities below and above the phase change
! We find this iteratively by first setting up density without phase changes and then
! adding the phase changes

! Note that this process is nonlinear - the scaling of gamma depends on total temperature
! contrast deltaT_dim, which itself depends on the adiabatic temperature gradient. Use
! simple iteration to resolve nonlinearity.

! PvK 091211
subroutine update_rho_adia(ichoice)
use eos
use sepmodulecomio
use coeff
use convparam
use mparallel
use control
implicit none
integer,intent(in) :: ichoice

integer,dimension(10) :: lo_index_ph=0,hi_index_ph=0
integer :: iph,i,npha,N,iz,lu,find_eos_iz,nph_here
integer :: nphmax
parameter(N=2000)
real(kind=8), dimension(10) :: rho1,rho2
real(kind=8), dimension(N) :: rho,temp_d
real(kind=8) :: alpha_mod,cp_mod,prespi,bigGammabar,dGdpi,dGdpi_d,get_prespi
real(kind=8) :: drhoph,drhophtot,rl,adiabat,Di_dim,Di_frac,bigG,y,funccf,find_phase_rho1
real(kind=8) :: rhoph,alpha,cp,z,dz,nph_orig,Hl_fac,rho_loc,get_drhodz_dim,drhodz
character(len=80) :: fname


!if (eos_type==0) then
!   if (print_node) then
!      write(irefwr,*) 'PINFO(update_rho): returning because eos_type=',eos_type
!      write(irefwr,*) 'with this eos_type the table look up is ignored'
!   endif
!   return
!endif
if (.not.delta1K.and..not.Tbar0) then
   if (nph>0) then
     if (print_node) then
        write(irefwr,*) 'PERROR(update_rho): '
        write(irefwr,*) 'only suited for delta1K when using phase changes and .not.Tbar0'
     endif
     call instop
   endif
endif

if (nph<=0) then
   if (print_node) then
      write(irefwr,*) 'PINFO(update_rho): returning because nph=',nph
   endif
   return
endif

if (nph>1) then
   if (print_node) then
      write(irefwr,*) 'PINFO(update_rho): update for nph>1'
   endif
   call instop
endif

if (update_phase_Tbar) then
   if (print_node) then
      write(irefwr,*) 'need to update to bring phase change back in Tbar'
   endif
   call instop
endif

! First dimensional reconstruction of reference state for density without phase changes
!!!!eos_rho_d(1)=rho_dim  ! keep in mind that you may want to update this for rho_d(z=0)<>rho_dim
! Do not overwrite eos_rho_d(1)
nph_here=0
bigG=0.0_8
dGdpi=0.0_8
eos_Gamma(iph,1:eos_np)=0.0_8
do iph=1,nph
   ! Find rho1 and rho2 around each phase change for use in latent heating term; used to construct Rb1 for heat equation
   phase_rho(1:2,iph)=eos_rho_d(1)
   phase_rho(1,iph)=find_phase_rho1(iph,eos_rho_d,eos_z_d,phz0_d(iph),eos_np)
   phase_rho(2,iph)=phase_rho(1,iph)+drho_ph_d(iph)
   do i=2,eos_np
      dz=(eos_z_d(i)-eos_z_d(i-1))/height_dim
      y=(height_dim-eos_z_d(i))/height_dim
      eos_alpha_d(i-1)=funccf(4,y,y,y)
      !drhodz_d=get_drhodz_dim(eos_z_d(i-1),eos_rho_d(i-1),eos_alpha_d(i-1),nph_here)
      drhodz=Di*eos_rho_d(i)
      eos_rho_d(i)=eos_rho_d(i-1)+drhodz*dz
      prespi=(eos_z_d(i)-phz0_d(iph))/height_dim
      bigG=bigGammabar(prespi,phdz_d(iph)/height_dim)
      if (update_phase_rho_via_gradient) then
         dGdpi=2.0_8*height_dim/phdz_d(iph)*bigG*(1.0_8-bigG)
         ! remove bigGamma from rhobar
         eos_rho_d(i)=eos_rho_d(i)+drho_ph_d(iph)*dGdpi*dz
         !if (bigG>0.4.and.bigG<0.6) then
         !   write(6,*) 'bigG: ',bigG,dGdpi,drho_ph_d(iph)*dGdpi*dz,eos_rho_d(i)
         !endif
      else
         ! just add delta_rho_d*bigGamma
         eos_rho_d(i)=eos_rho_d(i)+bigG*drho_ph_d(iph)
      endif
      eos_Gamma(iph,i)=bigG
   enddo
   !write(6,'(i5,4e15.7,:)') i,eos_rho_d(i-1),drhodz,dz,eos_rho_d(i)
enddo
!write(6,*) 'eos_rho_d(eos_np) = ',eos_rho_d(eos_np)
eos_rho=eos_rho_d/rho_dim


end subroutine update_rho_adia

real(kind=8) function get_drhodz_dim(z_d,rho_d,alpha_d,nph_here)
use eos
use convparam
use coeff
implicit none
real(kind=8) :: z_d,rho_d,alpha_d,z_nd
integer :: nph_here

z_nd=z_d/height_dim
get_drhodz_dim=Di*exp(Di*z_nd)*rho_dim

end function get_drhodz_dim

real(kind=8) function bigGammabar(prespi,d)
implicit none
real(kind=8) :: prespi,d

bigGammabar=0.5*(1+tanh(prespi/d))

end function bigGammabar

real(kind=8) function find_phase_rho1(iph,eos_rho_d,eos_z_d,phz0,eos_np)
implicit none
integer,intent(in) :: iph,eos_np
real(kind=8),intent(in) :: eos_rho_d(*),eos_z_d(*),phz0
integer :: i
real(kind=8) :: dz,rl

!write(6,*) 'find_phase_rho1: ',eos_np,phz0
find_phase_rho1=-1.0_8
i=1
do
   i=i+1
   if (i>eos_np) exit
   if (eos_z_d(i)>phz0) then
      dz=phz0-eos_z_d(i-1)
      rl=dz/(eos_z_d(i)-eos_z_d(i-1))
      find_phase_rho1=eos_rho_d(i-1)+rl*(eos_rho_d(i)-eos_rho_d(i-1))
      exit
   endif
enddo

end function find_phase_rho1
