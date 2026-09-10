! COMPRESSIONWORK
! Compute the amount of work done against hydrostatic pressure
! Make sure to define the variability of alpha,rho,c_p in exactly
! the same way as specified in coef800().
!
! Summer 2021. Cian's Compositional ALA analysis
! From just the Stokes equation we can formulate the equality
!     PHI = WORK
! where WORK is -w* [ W1 + W2 + W3 + W4 + W5 ] with
!     W1= (Di*cpr*rhobar*chiTbar*p')/(gamma_r*cvr) 
!     W2 = + Ra*alphabar*rhobar*T' 
!     W3 = - Rb*Gamma'
!     W4 = + delta_rho/rhobar*dGdpi
!     W5 = - (Di*alphabar*w*p')/(gamma_r*cpbar*gammabar)
! where T' may be T-Tbar or just T and Gamma' may be Gamma or Gamma-Gammabar
!
! In the case of chiTbar=1/rhobar and cpbar=gammabar=1 we have W5+W1=0 and:
! 
! PHI = int_omega Di/Ra tau : grad v
! WORK   = w* int_omega [ V1 + V2 + V3 ] where
!     V1 = Di alphabar rhobar T
!     V2 = - Di/Ra Delta_rho_i / rhobar 2/d_i Gammabar ( 1 - Gammabar ) p'
!     V3 = - Di Rb_i / Ra Gamma'
!
!
! PvK 091805
subroutine compressionwork()
use sepmodulekmesh
use sepmodulekprob
use control
use sepran_arrays
implicit none
interface 
   subroutine compressw01(npoint,usol1,usol2,work,rhoa,alpha,adia,pprime,coor,indprf,kprobf,indprp,kprobp,nunkp,nphys)
      integer,intent(in) :: npoint,indprf,kprobf(*),nunkp,indprp,nphys,kprobp(npoint,nphys)
      real(kind=8),intent(in) :: usol1(*),usol2(*),coor(2,*),rhoa(*),alpha(*),pprime(*),adia(*)
      real(kind=8),intent(inout) :: work(*)  
   end subroutine compressw01
   subroutine compressw02(npoint,usol1,usol2,work,rhoa,alpha,adia,press,gradT,coor,indprf,kprobf,indprp,kprobp,nunkp,nphys)
      integer,intent(in) :: npoint,indprf,kprobf(*),nunkp,indprp,nphys,kprobp(npoint,nphys)
      real(kind=8),intent(in) :: usol1(*),usol2(*),coor(2,*),rhoa(*),alpha(*),gradT(*),adia(*),press(*)
      real(kind=8),intent(inout) :: work(*)
   end subroutine compressw02

end interface

call sepactsolbf1(isol(1))
!write(6,*) 'indprfi: ',indprfi,indprpi,kprobpi,nunkp,nphys,coor(1,1)
if (cian_ALA) then
!pvk    ! use  Cian's write up
!pvk    write(6,*) 'isol: ',isol(1),isol(2),icompwork,irho,ialpha
!pvk    write(6,*) 'V1-5: ',(work_components(i),i=1,5)
!pvk    write(6,*) 'adia: ',iadia,ipress,igradT
!pvk    call compressw_cian(npoint,ks(isol(1))%sol,ks(isol(2))%sol,ks(icompwork)%sol,ks(irho)%sol,ks(ialpha)%sol, &
!pvk         & ks(work_components(1))%sol,ks(work_components(2))%sol,ks(work_components(3))%sol, &
!pvk         & ks(work_components(4))%sol,ks(work_components(5))%sol, &
!pvk        & ks(iadia)%sol,ks(ipress)%sol,ks(igradT)%sol,coor,indprfi,kprobfi,indprpi,kprobpi,nunkp,nphys)

else if (CY85_Eq22) then
   call compressw01(npoint,ks(isol(1))%sol,ks(isol(2))%sol,ks(icompwork)%sol,ks(irho)%sol,ks(ialpha)%sol, &
        & ks(iadia)%sol,ks(ipress)%sol,coor,indprfi,kprobfi,indprpi,kprobpi,nunkp,nphys)
else
  !write(6,*) 'compressw02: ',isol(1),isol(2),icompwork,irho,ialpha,iadia,ipress,igradT
  call compressw02(npoint,ks(isol(1))%sol,ks(isol(2))%sol,ks(icompwork)%sol,ks(irho)%sol,ks(ialpha)%sol, &
       & ks(iadia)%sol,ks(ipress)%sol,ks(igradT)%sol,coor,indprfi,kprobfi,indprpi,kprobpi,nunkp,nphys)
endif

end subroutine compressionwork

subroutine compressw_cian(npoint,usol1,usol2,work,rhoa,alpha,V1,V2,V3,V4,V5, &
    & adia,press,gradT,coor,indprf,kprobf,indprp,kprobp,nunkp,nphys)
use coeff
use convparam
use geometry
use control
implicit none
integer :: npoint,indprf,kprobf(*),nunkp,indprp,nphys,kprobp(npoint,nphys)
real(kind=8) :: usol1(*),usol2(*),work(*),coor(2,*),rhoa(*),alpha(*),gradT(*),adia(*),press(*)
real(kind=8) :: V1(*),V2(*),V3(*),V4(*),V5(*)
integer :: i,j1,j2,iph
real(kind=8) :: x,y,z,funccf,w,temp,rho,alpha_here,r,th,cost
real(kind=8) :: u,v,dGdpi,bigGamma,get_prespi,prespi,prespi2
real(kind=8) :: prespi_norm,bigGamma_norm,temp_norm,fulltemp,dTdx,dTdy

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
      u=usol1(j1)
      v=usol1(j2)
      w=v
      z=1-y
   endif
   
   if (compress) then
      rho=rhoa(i)
   else
      rho=1.0_8
   endif
   alpha_here=alpha(i)

!  First vertical advective transport of total temperature. With Kequivalent T0_dim is accounted for in the boundary conditions.
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
     if (compute_work_Tprime) fulltemp=fulltemp-adia(i)
   endif
   ! note multiplication by Ra because of future division by it (efficient...)
   V1(i) = rho*alpha_here*w*Ra*fulltemp

   ! Now the effect of the phase changes:
   ! - Di/Ra*rho*Rb*gamma*(T+T0)*DGamma/Dt = + Di/Ra* rho Rb gamma (T+T0) dGamma/dpi ( gamma u.gradT + w)
   ! first this term with u.gradT==0
   bigGamma=0.0_8
   temp=usol2(i) 
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
      if (subtract_bigGamma_norm) then
         temp_norm=pht0(iph)
         prespi_norm=get_prespi(z,temp_norm,iph)
         bigGamma_norm=0.5_8+0.5_8*tanh(prespi_norm/phdz(iph))
      else
         bigGamma_norm=0.0_8
      endif
 
      ! This is equivalent to adding the effective alpha term only
      ! work(i) = work(i) + rho*Rb1(iph)*gamma(iph)*fulltemp*dGdpi*w
      ! This is more complete
      if (rhobar_times_bigGamma) then
         work(i) = work(i) + rho*Rb1(iph)*gamma(iph)*fulltemp*dGdpi*(gamma(iph)*(u*dTdx+v*dTdy)+rho*w)
      else
         work(i) = work(i) + Rb1(iph)*gamma(iph)*fulltemp*dGdpi*(gamma(iph)*(u*dTdx+v*dTdy)+rho*w)
      endif
   enddo


   !if (y>0.45.and.y<0.55.and.x<1e-6) then
   !write(6,'(''rho etc.: '',11e13.6)') y,temp,pht0(1),gamma(1),w,phz0(1),phdz(1),prespi,prespi2,bigGamma,work(i)
   !endif

   work(i) = V1(i)+V2(i)+V3(i)

enddo


end subroutine compressw_cian

