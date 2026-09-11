subroutine fillcoef800(iuser,user)
use sepmodulecomio
use sepmodulekmesh
use sepmodulekprob
use control
use coeff
use eos
use sepran_interface
use sepran_arrays
implicit none
!interface
!  subroutine pecopy(ichoice,user,isol)
!     integer,intent(in) :: ichoice,isol
!     real(kind=8),intent(inout) :: user(:)
!  end subroutine pecopy
!end interface
integer :: iuser(*)
real(kind=8) :: user(*)
!integer, parameter :: nph_max=10
!real(kind=8),dimension(nph_max) :: glRbRa,gl2RbRaDi,sign_gamma
real(kind=8) :: prespi,get_prespi,bigGamma,bigGamma_norm,prespi_norm
real(kind=8) :: alpha_mod,cp_mod,alphamax,cpsmin,cpsmax
integer :: ihorvel,ivervel,ibeta,irhocp,iright,itemp,iuser_cond
integer :: ip,i,iph
logical :: first=.true.
save first

if (nph>NPHASE_MAX) then
   if (print_node) then
      write(irefwr,*) 'PERROR(fillcoef800): nph>NPHASE_MAX'
      write(irefwr,*) 'nph, NPHASE_MAX: ',nph,NPHASE_MAX
   endif
   call instop
endif

do iph=1,nph
  if (Ra.eq.0) then
     glRbRa(iph)=0
  else
     glRbRa(iph) = gamma(iph)*Rb1(iph)/Ra ! note use of Rb1 here. 
  endif
  gl2RbRaDi(iph) = gamma(iph)*Di*glRbRa(iph)
! write(6,'(''Rb1 : '',3e15.7)') Rb1(iph),glRbRa(iph),gl2RbRaDi(iph)

  if (phdz(iph).le.0.or.phdz(iph).ge.1) then
     write(6,*) 'PERROR(coef8_2): phdz <= 0 or phdz >= 1'
     write(6,'(''phdz('',i2,'') = '',f12.3)') iph,phdz(iph)
     call instop
  endif
enddo

do iph=1,nph
   if (abs(gamma(iph)) > 1e-7) then
      sign_gamma(iph) = gamma(iph)/abs(gamma(iph))
   else
      sign_gamma(iph) = 0d0
   endif
enddo



if (pedebug.and.print_node) write(irefwr,*) 'fillcoef800: npoint=',npoint
iuser(2:iuser(1))=0
user(2:nint(user(1)))=0.0_8
iuser(2)=1
! Define coefficients through iuser/user (bit convoluted but robust)
iuser(6)=15 ! start information about first and only element group at iuser position 15

if (pedebug.and.print_node) write(irefwr,*) 'fillcoef800: conductive=',conductive

! We will adopt the comprS approach to specify all parameters pointwise even if constant
! pointers to real information on velocity components, beta, right hand side, rho_cp, temperature, and conductivity
iuser_cond=10
ihorvel=10+npoint
ivervel=10+2*npoint
ibeta=10+3*npoint
iright=10+4*npoint
irhocp=10+5*npoint
itemp=10+6*npoint

! first define integer information for coefficients
ip=14
! (1) not used
iuser(ip+1) = 0  
! (2) type of upwinding ; avoid with P2 for now
iuser(ip+2) = metupw
! (3) intrule
iuser(ip+3) = intrule800+100*interpol800  ! let element routines define integration 
! (4) icoor: type of coordinate system
iuser(ip+4) = icoor800    ! Cartesian
! (5) not yet used
iuser(ip+5) = 0    

! (6) k11
iuser(ip+6) = 2001
iuser(ip+7) = iuser_cond
ip=ip+1 ! account for increment since type 2001 takes two integer spaces
! (7) k12
iuser(ip+7) = 0
! (8) k13
iuser(ip+8) = 0
! (9) k22
iuser(ip+9) = 2001
iuser(ip+10) = iuser_cond
ip=ip+1
! (10) k23
iuser(ip+10) = 0
! (11) k33
iuser(ip+11) = 0

! velocity components
if (conductive) then
  iuser(ip+12) = 0  ! u=0
  iuser(ip+13) = 0  ! v=0
else
  ! (12) u
  iuser(ip+12) = 2001
  iuser(ip+13) = ihorvel
  ip=ip+1 ! account for increment since type 2001 takes two integer spaces
  ! (13) v
  iuser(ip+13) = 2001
  iuser(ip+14) = ivervel
  ip=ip+1 ! account for increment since type 2001 takes two integer spaces
endif

! (14) w==0
iuser(ip+14) = 0

! (15) beta
iuser(ip+15) = 2001
iuser(ip+16) = ibeta
ip=ip+1 ! account for increment since type 2001 takes two integer spaces

! (16) q    Note we use iright because we may add viscous dissipation to form
!           the right-hand side
iuser(ip+16)=2001
iuser(ip+17)=iright
ip=ip+1 ! account for increment since type 2001 takes two integer spaces

! (17) rho*cp
iuser(ip+17) = 2001
iuser(ip+18) = irhocp
ip=ip+1
! make sure rest of array is 0
!iuser(ip+18:iuser(1))=0

! constant conductivity for now; replace with copy of conductivity vector
!!!!!user(iuser_cond:iuser_cond+npoint-1)=1.0_8
call pecopy(0,user(iuser_cond:iuser_cond+npoint-1),icond)
!write(irefwr,*) 'conductivity: '
!write(irefwr,*) user(iuser_cond:iuser_cond+npoint-1:npoint/10)
if (.not.conductive) then
   call pecopy(2,user(ihorvel:ihorvel+npoint-1),isol1) 
   call pecopy(3,user(ivervel:ivervel+npoint-1),isol1) 
endif
!write(irefwr,*) 'maxval(ivervel): ',maxval(user(ivervel:ivervel+npoint-1))

! copy temperature into user
call pecopy(0,user(itemp:itemp+npoint-1),isol2)

if (iqtype /= 0) then
   ! find heat production
   call findheatgen(idens,iheat)
endif

if (Di>0) then
   ! find viscous dissipation
   call getvisdip
endif

! copy phi + heat production into coefficients for right hand side vector
call pecophi(user(iright:iright+npoint-1),iphi,iheat)
!write(irefwr,*) 'user(iright): ',user(iright)

! assemble pointwise coefficients 
call fill800_beta_phase(npoint,user(ihorvel),user(ivervel),user(ibeta),user(iright),user(irhocp),user(itemp), &
     & user(iuser_cond),coor,ks(iadia)%sol,ks(ialpha)%sol,ks(icp)%sol,ks(irho)%sol,ks(isol2)%sol)
!write(irefwr,*) 'user(iright): ',user(iright)


! keep track of whether there is real information in user
iuser(7)=iuser_cond
iuser(8)=ihorvel
iuser(9)=ivervel
iuser(10)=ibeta
iuser(11)=iright
iuser(12)=irhocp
if (pedebug.and.print_node) write(irefwr,*) 'fillcoef800: ihorvel/y=',ihorvel,ivervel

if (first .and. .not. conductive .and. print_node) then
   if (.not.petest) first=.false.
   do i=1,100
      write(irefwr,'(i20,i20,e15.7)') i,iuser(i),user(i)
   enddo
   call coef800_print(npoint,coor,user(ihorvel),user(ivervel),user(ibeta),user(iright),user(irhocp),user(itemp),user(iuser_cond))
   !write(6,*) 'stop after coef800_print'
   !call instop
endif


end subroutine fillcoef800

subroutine fill800_beta_phase(N,horvel,vervel,beta,right,rhocp,temp,cond,coor,adia,alpha,cp,rhoa,temp2)
use sepmodulecomio
use coeff
use mdebugTbars
use geometry
use control
use eos
use convparam
implicit none
integer :: N
real(kind=8),dimension(*) :: horvel,vervel,beta,right,rhocp,temp,cond,adia,alpha,cp,rhoa,temp2
real(kind=8) :: coor(2,N)
integer :: i,iph
real(kind=8) :: radial_velocity,x,y,r,cost,sint,th,alphas,rho,cps,rhoalphas,rhocps,rhoprime
real(kind=8) :: betamin,betamax,rightmin,rightmax,vervelmin,vervelmax,alphamin,alphamax,cpsmin,cpsmax
real(kind=8) :: funccf,presp,bigGamma,dGdpi,get_prespibar
real(kind=8) :: rhomminh,rhomaxh,alpha_mod,cp_mod,tempmin,tempmax,z,prespi,get_prespi,temp_norm,prespi_norm
real(kind=8) :: bigGamma_norm,tratio,rightW

!if (cyl) then
!   write(irefwr,*) 'PERROR(coef8_2): not yet suited for cyl'
!   call instop
!endif
!!if (nph>0.and.Di>1e-7) then
!   if (print_node) write(irefwr,*) 'PERROR(coef8_2): not yet suited for phase changes and Di>0'
!   call instop
!endif

!write(irefwr,*) 'density in fill800_beta_phase: ',minval(rhoa(1:N)),maxval(rhoa(1:N)),N
do iph=1,nph
   if (Ra==0) then
       glRbRa(iph)=0.0_8
   else 
       glRbRa(iph)=gamma(iph)*Rb1(iph)/Ra
   endif
  gl2RbRaDi(iph) = gamma(iph)*Di*glRbRa(iph)
! write(6,'(''Rb1 : '',3e15.7)') Rb1(iph),glRbRa(iph),gl2RbRaDi(iph)

  if (phdz(iph).le.0.or.phdz(iph).ge.1) then
     write(6,*) 'PERROR(coef8_2): phdz <= 0 or phdz >= 1'
     write(6,'(''phdz('',i2,'') = '',f12.3)') iph,phdz(iph)
     call instop
  endif
enddo
do iph=1,nph
   if (abs(gamma(iph)) > 1e-7) then
      sign_gamma(iph) = gamma(iph)/abs(gamma(iph))
   else
      sign_gamma(iph) = 0d0
   endif
enddo

alphamin=1e9
alphamax=-1e9
cpsmin=1e9
cpsmax=-1e9
rhomminh=1e9
rhomaxh=-1e9
do i=1,N
   alpha_mod=0.0_8
   cp_mod=0.0_8
   tempmin=1e5
   tempmax=-1e5
   x=coor(1,i)
   y=coor(2,i)
   if (x<1e-6.and.y>0.49999.and.y<0.500001) then
      visdiph=right(i)
   endif
   if (cyl) then
      r=sqrt(x*x+y*y)
      cost=y/r
      sint=x/r
      if (x>=0) then
         th=acos(cost)
      else
         th = 2*pi - acos(cost)
      endif
      radial_velocity = horvel(i)*sin(th) + vervel(i)*cos(th)
   else
      radial_velocity = vervel(i)
   endif

   dGdpi=0.0_8

   ! find local specific heat
   cps=1.0_8
   if (compress) then
      if (eos_type /= 0) then
         cps=cp(i)
      endif
   endif

   ! find local density
   if (compress) then
      rho = rhoa(i)
   else
      rho = 1.0_8
   endif

   !if (ialphatype>=1) then
   ! always fill alpha(:) even if constant 1
   alphas=alpha(i)
   !else
   !   alphas=1.0_8
   !endif

   rhoalphas=rho*alphas
   rhocps=rho*cps
   rhoprime=rho

   if (.not.ignore_latent_heat) then
!    Modify rho*cp and rho*alpha for latent heat of phase changes
     do iph=1,nph
!      *** make sure r2-r1=1 (for cylindrical/axisymmetric case)
       if (cyl) then
          !prespi = (r2-r) - phz0(iph) - gamma(iph)*(temp(i)-pht0(iph))
          prespi = get_prespi(r2-r,temp(i),iph)
          z=r2-r
       else
          !prespi = (1-y) - phz0(iph) - gamma(iph)*(temp(i)-pht0(iph))
          prespi = get_prespi(1-y,temp(i),iph)
          z=1-y
       endif
!      if (eos_type==0) then
         ! for compatibility sake, I guess (adia==0)
         ! bigGamma_norm is ignored...
         !temp_norm = adia(i)
         !prespi_norm = z - phz0(iph) - gamma(iph)*(temp_norm-pht0(iph))
         !prespi_norm=get_prespi(z,temp_norm,iph)
         !bigGamma_norm = 0.5d0 + 0.5d0*tanh(prespi_norm/phdz(iph))
!      else 
!        if (print_node) write(irefwr,*) 'fillcoef800 needs update for phase change with eos_type<>0'
!        call instop
!      endif
       ! First add Ts; it is included in the b.c. for Kequivalent
       if ((compress.or.EBA).and..not.Kequivalent) then
          tratio=temp(i)+Ts_nondimK+dCY85_nondim
       else if (EBA.and.Kequivalent) then
          ! catch for CY85
          tratio=temp(i)+dCY85_nondim
       endif
       if (use_effective_alpha_Hr) then
          ! modify only alpha using Gammabar
          prespi=get_prespibar(z,iph)
          bigGamma = 0.5 + 0.5*tanh(prespi/phdz(iph))
          dGdpi = 2d0/phdz(iph)*bigGamma*(1-bigGamma)
          alpha_mod = alpha_mod + glRbRa(iph)*dGdpi  ! CW notes 
       else
          ! use full formulation with modified alpha and cp using phase boundary topography
          ! Note: prespi depends on depth and temperature
          bigGamma = 0.5 + 0.5*tanh(prespi/phdz(iph))
          dGdpi = 2d0/phdz(iph)*bigGamma*(1-bigGamma)
          if (compress) then
            !alpha_mod = alpha_mod + rho*glRbRa(iph)*dGdpi  ! added rho 091811
            alpha_mod = alpha_mod + glRbRa(iph)*dGdpi  ! CW notes 
          else
            alpha_mod = alpha_mod + glRbRa(iph)*dGdpi
          endif
          cp_mod = cp_mod + gl2RbRaDi(iph)*dGdpi*tratio
       endif
       alphamin=min(alpha_mod,alphamin)
       alphamax=max(alpha_mod,alphamax)
       cpsmin=min(cp_mod,cpsmin)
       cpsmax=max(cp_mod,cpsmax)
     enddo
     if (do_rhoprime) then
        ! Note rhoprime=rho
        rhocps = rhoprime*(cps + cp_mod)
        rhoalphas = rhoprime*(alphas + alpha_mod)
     else
        rhocps = rho*(cps + cp_mod)
        rhoalphas = rho*(alphas + alpha_mod)
     endif
   endif




   beta(i)=rhoalphas*Di*radial_velocity
   !if (i==1) write(irefwr,*) 'right(1)=',right(1)
   ! add T0 from work term 
   if (EBA) then 
      ! add negative effect of T0 in the rhsd
      right(i)=right(i)-rhoalphas*Di*radial_velocity*(Ts_nondimK+dCY85_nondim)
   else if (compress) then
      ! add effect of adiabat with Ta/=Ts when solving for full T
      if (Kequivalent) then
         ! don't do anything because T0 is in the solution
         !rightW=-rhoalphas*Di*radial_velocity*Tso_nondimK
      else if (solve_for_Tperturb) then
         ! add effect of T0 and shifted adiabat (if necessary)
         right(i) = right(i) + rhoalphas*Di*radial_velocity*(Tbars_nondimK-Ts_nondimK-dCY85_nondim)
      else
         ! Add just effect of T0
         right(i) = right(i) + rhoalphas*Di*radial_velocity*(-Ts_nondimK-dCY85_nondim)
      endif
   endif
   if (compress.and.solve_for_Tperturb) then
      ! add diffusion up the adiabat
      ! this is only correct if k=1
      right(i)=right(i)+d2adiabatdz2(y)
      !if (x<1e-6.and.y>0.49999.and.y<0.500001) d2Tbardz2h=d2adiabatdz2(y)
   endif
   !if (x<1e-6.and.y>0.49999.and.y<0.500001) then
   !   betah=beta(i)
   !   rhoalphah=rhoalphas
   !   righth=right(i)
   !   adiah=adia(i)
   !   temph=temp2(i)
   !   upwh = radial_velocity
   !   work1h=betah*temph-rightW
   !endif
   !if (i==1) write(irefwr,*) 'right(1)=',right(1)
   rhocp(i)=rhocps
   !if(x<1e-6.and.y>0.4.and.y<0.6) write(irefwr,'(9e12.5)') y,rhoalphas,rho*alphas,rhocps,rho*cps,beta(i),right(i),radial_velocity,rho


enddo
!write(6,*) 'Ts_nondimK: ',Ts_nondimK
!write(6,*) 'Tbars_nondimK: ',Tbars_nondimK
!write(6,*) 'dCY85_nondim: ',dCY85_nondim
!write(irefwr,*) 'rhoalphas, Di: ',rhoalphas,Di

!betamin=minval(beta(1:N))
!betamax=maxval(beta(1:N))
!write(irefwr,*) 'betamin/max: ',betamin,betamax
!rightmin=minval(right(1:N))
!rightmax=maxval(right(1:N))
!write(irefwr,*) 'rightmin/max: ',rightmin,rightmax
!vervelmax=maxval(vervel(1:N))
!vervelmin=minval(vervel(1:N))
!write(irefwr,*) 'vervelmin/max: ',vervelmin,vervelmax
!if (niter==1) call instop

end subroutine fill800_beta_phase

