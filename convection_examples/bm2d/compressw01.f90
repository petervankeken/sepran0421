subroutine compressw01(npoint,usol1,usol2,work,rhoa,alpha,adia,pprime,coor,indprf,kprobf,indprp,kprobp,nunkp,nphys)
use coeff
use convparam
use geometry
use control
implicit none
integer,intent(in) :: npoint,indprf,kprobf(*),nunkp,indprp,nphys,kprobp(npoint,nphys)
real(kind=8),intent(in) :: usol1(*),usol2(*),coor(2,*),rhoa(*),alpha(*),pprime(*),adia(*)
real(kind=8),intent(inout) :: work(*)
integer :: i,j1,j2,iph
real(kind=8) :: x,y,z,funccf,w,temp,rho,alpha_here,r,th,cost
real(kind=8) :: u,v,dGdpi,bigGamma,get_prespi,prespi,prespi2
real(kind=8) :: prespi_norm,bigGamma_norm,temp_norm,fulltemp,dTdx,dTdy
real(kind=8) :: delta_rho

!write(6,*) 'cp  : ',cp(1)
!write(6,*) 'kprobf: ',kprobf(1),kprobp(1,1)

z=0.0_8
do i=1,npoint
   work(i)=0.0_8
   if (indprf == 0 .and. indprp == 0) then
      j1=(i-1)*nunkp+1
      j2=j1+1
   else if (indprp/=0) then
      j1=kprobp(i,1)
      j2=kprobp(i,2)
   else
      j1=kprobf(i)+1
      j2=j1+1
   endif
   x=coor(1,i)
   y=coor(2,i)
   if (cyl) then
      r=sqrt(x*x+y*y)
      cost=y/r
      if (x>=0) then
         th=acos(cost)
      else
         th=2*pi-acos(cost)
      endif
      u=usol1(j1)
      v=usol1(j2)
      w=u*sin(th)+v*cos(th)
      z=r2-r
   else
      u = usol1(j1)
      w = usol1(j2)
      z=1-y
   endif
   
   if (compress) then
      rho=rhoa(i)
   else
      rho=1.0_8
   endif
   alpha_here=alpha(i)

   temp=usol2(i) 
   !do not do this in computation of phase function!
   !!!!if (compute_work_Tprime) temp=usol2(i)-adia(i)
   bigGamma=0.0_8
   do iph=1,nph
     ! non-dimensional density change
     delta_rho=drho_ph_d(iph)/rho_dim
     if (cyl) then
         !prespi = (r2-r) - phz0(iph) - gamma(iph)*(temp-pht0(iph))
         z=r2-r
      else
         z=1.0_8-y
         !prespi = (1-y) - phz0(iph) - gamma(iph)*(temp-pht0(iph))
      endif
      prespi=get_prespi(z,temp,iph)
      !prespi = z-phz0(iph)-gamma(iph)*(temp-pht0(iph))
      bigGamma = 0.5 + 0.5*tanh(prespi/phdz(iph))
      dGdpi = 2d0/phdz(iph)*bigGamma*(1d0-bigGamma)
      ! make it optional to subtract Gamma_bar from work term
      ! needs update for ALA
      if (subtract_bigGamma_norm) then
         temp_norm=pht0(iph)
         prespi_norm=get_prespi(z,temp_norm,iph)
         bigGamma_norm=0.5_8+0.5_8*tanh(prespi_norm/phdz(iph))
      else
         bigGamma_norm=0.0_8
      endif
 
!     DO NOT do the modifications for alpha & cp below
!     (Christensen & Yuen, JGR, 1985, eqs 20+21) since those
!     are for implementing latent heat into heat equation and
!     actually part of the mechanical work terms
!     DO NOT USE alpha = alpha + glRbRa(iph)*dGdpi
!     DO NOT USE c_p = c_p + gl2RbRaDi(iph)*dGdpi*tratio

      ! PvK 050321 - This is CY85 Equation 22
      ! Note that work is multiplied by Di/Ra outside of this subroutine so this line really reads:
      ! work(i) = work(i) - rho*w*Di*Rb(iph)/Ra*bigGamma
      ! Note you should use Rb not Rb1 here (that is just for the latent heating effect)
      if (rhobar_times_bigGamma) then
         work(i) = work(i) - rho*w*Rb(iph)*(bigGamma-bigGamma_norm)
      else
         work(i) = work(i) - w*Rb(iph)*(bigGamma-bigGamma_norm)
      endif

      if (compress.and.add_dynamic_pressure_term_to_work) then
         ! Now add the phase change term that depends on dynamic pressure. 
         ! note sign: 
         work(i) = work(i) + 1.0_8/arbitrary_p_scale*delta_rho/rho*w*dGdpi*pprime(i) ! *Ra/Di
      endif

      
   enddo

!  vertical advective transport of total temperature. With Kequivalent T0_dim is accounted for in the boundary conditions.
   fulltemp=usol2(i)+dCY85_nondim
   if (EBA) then
     if (delta1K.and..not.Kequivalent) then
        fulltemp = usol2(i) + Ts_nondimK + dCY85_nondim
     else if (.not.Kequivalent) then
        fulltemp = usol2(i) + Ts_nondimK + dCY85_nondim
     endif
   else if (compress) then
     if (Kequivalent) then
        !write(6,*) 'not yet suited for compress and Kequivalent'
        !call instop
        fulltemp=usol2(i)
        ! since Ts_nondimK is in the boundary condition you don't need to do anything here
     else if (solve_for_Tperturb) then
           ! T' notation. Work = T'-Tbar_s+Ts
           fulltemp = usol2(i) - Tbars_nondimK + Ts_nondimK
     else  
          ! full T but not Kequivalent
          fulltemp = usol2(i) + Ts_nondimK
     endif
   endif
   if (compute_work_Tprime) fulltemp=fulltemp-adia(i)
   ! note multiplication by Ra because of future division by it (efficient)
   work(i) = work(i) + rho*alpha_here*w*Ra*fulltemp

   !if (y>0.45.and.y<0.55.and.x<1e-6) then
   !write(6,'(''rho etc.: '',11e13.6)') y,temp,pht0(1),gamma(1),w,phz0(1),phdz(1),prespi,prespi2,bigGamma,work(i)
   !endif


enddo


end subroutine compressw01

