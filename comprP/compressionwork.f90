! COMPRESSIONWORK
! Compute the amount of work done against hydrostatic pressure
! Make sure to define the variability of alpha,rho,c_p in exactly
! the same way as specified in coef800().
!
! PvK 091805
subroutine compressionwork()
use sepmodulekmesh
use sepmodulekprob
use sepran_arrays
implicit none

call sepactsolbf1(isol1)
!write(6,*) 'indprfi: ',indprfi,indprpi,kprobpi,nunkp,nphys,coor(1,1)

call compressw01(npoint,ks(isol1)%sol,ks(isol2)%sol,ks(icompwork)%sol,ks(irho)%sol,ks(ialpha)%sol, &
       & coor,indprfi,kprobfi,indprpi,kprobpi,nunkp,nphys)

end subroutine compressionwork

subroutine compressw01(npoint,usol1,usol2,work,rhoa,alpha,coor,indprf,kprobf,indprp,kprobp,nunkp,nphys)
use coeff
use convparam
use geometry
use control
implicit none
integer :: npoint,indprf,kprobf(*),nunkp,indprp,nphys,kprobp(npoint,nphys)
real(kind=8) :: usol1(*),usol2(*),work(*),coor(2,*),rhoa(*),alpha(*)
integer :: i,j1,j2,iph
real(kind=8) :: x,y,z,funccf,w,temp,rho,c_p,alpha_here,r,th,cost
real(kind=8) :: u,v,dGdpi,bigGamma,get_prespi,prespi,prespi2
real(kind=8) :: prespi_norm,bigGamma_norm,temp_norm

!write(6,*) 'cp  : ',cp(1)
!write(6,*) 'kprobf: ',kprobf(1),kprobp(1,1)

c_p=1.0_8
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
      w = usol1(j2)
      z=1-y
   endif
   
   ! specific heat: unimportant here...
!  if (compress .and. eos_type/=0) then
!     c_p=cp(i)
!  endif
 
   if (compress) then
      rho=rhoa(i)
   else
      rho=1.0_8
   endif
   alpha_here=alpha(i)

   temp=usol2(i) 
   bigGamma=0.0_8
   do iph=1,nph
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
      if (EBA.and.subtract_bigGamma_norm) then
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
      work(i) = work(i) - rho*w*Rb(iph)*bigGamma
   enddo
!  total vertical advective transport. With Kequivalent T0_dim is accounted for in the boundary conditions.
   temp=usol2(i)+dCY85_nondim
   if (EBA) then
     if (delta1K.and..not.Kequivalent) then
        temp = usol2(i) + Ts_nondimK + dCY85_nondim
     else if (.not.Kequivalent) then
        temp = usol2(i) + Ts_nondimK + dCY85_nondim
     endif
   else if (compress) then
     if (Kequivalent) then
        !write(6,*) 'not yet suited for compress and Kequivalent'
        !call instop
        temp=usol2(i)
        ! since Ts_nondimK is in the boundary condition you don't need to do anything here
     else if (solve_for_Tperturb) then
           ! T' notation. Work = T'-Tbar_s+Ts
           temp = usol2(i) - Tbars_nondimK + Ts_nondimK
     else  
          ! full T but not Kequivalent
          temp = usol2(i) + Ts_nondimK
     endif
   endif
   work(i) = work(i) + rho*alpha_here*c_p*w*Ra*temp

   !if (y>0.45.and.y<0.55.and.x<1e-6) then
   !write(6,'(''rho etc.: '',11e13.6)') y,temp,pht0(1),gamma(1),w,phz0(1),phdz(1),prespi,prespi2,bigGamma,work(i)
   !endif


enddo


end subroutine compressw01
