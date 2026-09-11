! update_rho_adia
! Set up an EoS based on a lookup table, but modify density and adiabat by integrating
! the respective gradients. Take into account phase change influence on density and adiabat
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
subroutine update_rho_adia_old(ichoice)
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
real(kind=8) :: alpha_mod,cp_mod,prespi,bigGamma,dGdpi,dGdpi_d,get_prespi
real(kind=8) :: drhoph,drhophtot,rl,adiabat,Di_dim,Di_frac
real(kind=8) :: rhoph,alpha,cp,z,dz,nph_orig,Hl_fac,rho_loc
character(len=80) :: fname


if (eos_type==0) then
   if (print_node) then
      write(irefwr,*) 'PINFO(update_adiabat_info): returning because eos_type=',eos_type
      write(irefwr,*) 'with this eos_type the table look up is ignored'
   endif
   return
endif
if (.not.delta1K) then
   if (nph>0) then
     if (print_node) then
        write(irefwr,*) 'PERROR(update_adiabat_info): '
        write(irefwr,*) 'only suited for delta1K when using phase changes '
     endif
     call instop
   endif
endif

if (nph<=0) then
   if (print_node) then
      write(irefwr,*) 'PINFO(update_adiabat_info): returning because nph=',nph
   endif
   return
endif


if (ichoice==-1) then
!   special case: simulate Leng and Zhong, 2010
   if (print_node) then
      write(irefwr,*) 'PERROR(update_rho_adia): Leng&Zhong formulation needs significant updates'
   endif
   call instop
   !call lengzhong()
   return
endif

if (ichoice==1.or.ichoice==0.or.ignore_latent_heat) then
!   set up adiabat by integrating temperature gradient
!   and modify Rb(1..nph) for compressible flow, but do not
!   *** add latent heating effects to adiabat          
   Hl_fac = 0
else if (ichoice==2) then
!   set up adiabat by integrating temperature gradient
!   and modify Rb(1..nph) for compressible flow; 
!   add latent heating effects to adiabat          
   Hl_fac = 1
else
   if (print_node) then
     write(irefwr,*) 'PERROR(update_rho_adia): ichoice = ',ichoice
     write(irefwr,*) 'but valid options are -1, 0, 1, and 2 only'
   endif
   call instop
endif

if (N < eos_np) then
   write(irefwr,*) 'PERROR(find_rho_adia): temporary array rho '
   write(irefwr,*) ' is too small : ',N,eos_np
   call instop
endif

! Modify upper parts of table for rho to remove crust
if (eos_type == 2) then
  do i=1,8
     eos_rho_d(i)=eos_rho_d(9)
  enddo
endif

drho = drho_rel*rho_dim

open(11,file='buoyancy.numbers')
write(11,*) 'Ra   = ',Ra,alpha_dim,cp_dim,deltaT_dim
do iph=1,nph
  write(11,*) 'Rb,drho   = ',Rb(iph),drho(iph)
  write(11,'(i2,4e15.7)') iph,gamma(iph),glRbRa(iph), gl2RbRaDi(iph),Rb(iph)
enddo
close(11)

!      **** Correct for choice of Di given on input
Di_dim = alpha_dim*g_dim*(R2_dim-R1_dim)/cp_dim
Di_frac = Di/Di_dim
!      write(irefwr,*) 'Di_frac: ',Di_frac

!      iph=1
!      write(irefwr,*) 'iph: ',iph,phz0_d(iph),phdz_d(iph)
!      write(irefwr,*) 'eos_z_d: ',eos_z_d(1),eos_z_d(eos_np)

! Find locations of phase transitions in look up tables
do iph=1,nph
  do i=1,eos_np
   if ((eos_z_d(i+1).ge.(phz0_d(iph)-phdz_d(iph))/1e3) .and.  (eos_z_d(i).le.(phz0_d(iph)-phdz_d(iph))/1e3)) lo_index_ph(iph)=i
   if (eos_z_d(i+1)>=(phz0_d(iph)+phdz_d(iph))/1e3 .and.  eos_z_d(i)<=(phz0_d(iph)+phdz_d(iph))/1e3) hi_index_ph(iph)=i
 enddo
 write(irefwr,'(''phase change location: '',3i5,2f15.3,$)')  iph, lo_index_ph(iph),hi_index_ph(iph),phz0_d(iph),phdz_d(iph)
 write(irefwr,'(2f12.3)') eos_z_d(lo_index_ph(iph)),eos_z_d(hi_index_ph(iph))
 write(irefwr,'('' Rb, Ra, gamma: '',3f12.3)') Rb(iph),Ra,gamma(iph)
enddo

do iph=1,nph
  drho(iph)=drho_rel(iph)*rho_dim
enddo

rho2=rho_dim
rho1=rho_dim
Rb1=Rb
nph_orig=nph
nphmax=nph+1
if (ichoice==0) nphmax=0
do npha=0,nphmax
! skip phase change on first go through         
  drhophtot=0
  !!temp_d = eos_T_d(1)
  rho(1) = eos_rho_d(1)
  Rr=rho_dim*rho_dim/(rho2*rho1)
  Rb1=Rb*Rr
  glRbRa=gamma*Rb1/Ra
  gl2RbRaDi = gamma*Di*glRbRa
  write(irefwr,*) 'Rr: ',npha,Rr(1),Rr(2)
  do iz=1,eos_np-1
     z = eos_z_d(iz)*1e3
     dz = (eos_z_d(iz+1)-eos_z_d(iz))*1e3
     alpha = eos_alpha_d(iz)
     cp = eos_cp_d(iz)
     rho_loc = eos_rho_d(iz)/rho_dim

!    add phase change effect 
!    *** non-dimensionalize alpha,cp
     alpha_mod=alpha/alpha_dim
     cp_mod=cp/cp_dim
     drhoph=0d0
     do iph=1,min(npha,nph)
        !prespi = z-phz0_d(iph)
        prespi=get_prespi(z,0.0_8,iph)
        bigGamma = 0.5d0 + 0.5d0*tanh(prespi/phdz_d(iph))
        rhoph = bigGamma*drho(iph)
        dGdpi = 2d0/phdz(iph)*bigGamma*(1-bigGamma)
        dGdpi_d = 2d0/phdz_d(iph)*bigGamma*(1-bigGamma)
        drhoph = drhoph+drho(iph)*dGdpi_d
        alpha_mod = alpha_mod + Hl_fac*rho_loc*glRbRa(iph)*dGdpi ! rho added 091811
        cp_mod = cp_mod + Hl_fac*gl2RbRaDi(iph)*dGdpi*temp_d(iz)/deltaT_dim
        eos_Gamma(iph,iz)=bigGamma
     enddo ! iph
!    redimensionalize cp,alpha
     cp_mod=cp_mod*cp_dim
     alpha_mod=alpha_mod*alpha_dim
!    integrate density
     rho(iz+1)=rho(iz) +dz*(Di_frac*rho(iz)*g_dim*alpha/(cp*Grueneisen) +drhoph)
!    integrate temperature
     !!temp_d(iz+1) = temp_d(iz) + dz*Di_frac*alpha_mod*g_dim*temp_d(iz)/cp_mod
  enddo ! iz
  eos_Gamma(:,eos_np)=1d0
  if (npha<=nph) then
    rho1=rho(lo_index_ph)
    rho2=rho(hi_index_ph)
    write(irefwr,*) 'rho1: ',rho1(1:nph)
    write(irefwr,*) 'rho2: ',rho2(1:nph)
  endif
enddo !npha

Rb1=Rb*Rr
write(irefwr,*) 'Rr = ',Rr(1:nph)

do i=1,eos_np
  eos_T_d(i)=0.0 !! temp_d(i)
  eos_rho_d(i)=rho(i)
  eos_T(i) = (eos_T_d(i)-T0_dim)
enddo

if (print_node) then
  write(irefwr,*) 'update_rho_adia: eos_T_d  : ',eos_T_d(1), eos_T_d(eos_np)
  write(irefwr,*) 'update_rho_adia: eos_rho_d: ',eos_rho_d(1), eos_rho_d(eos_np)
endif



end subroutine update_rho_adia_old

