real(kind=8) function func(ichoice,x,y,z)
use sepmodulecomio
use geometry
use control
use convparam
use coeff
use tracers
use eos
implicit none
integer,intent(in) :: ichoice
real(kind=8),intent(in) :: x,y,z
real(kind=8) :: r,drh,theta,argsin,funccf,th,vx,vy,xp,yp
real(kind=8) :: xr,wpi,zz,ah,rara,u0,v0,qq
real(kind=8) :: tu,tl,tr,ts,trarg,tsarg,dzero=0.0_8,gx,fxy
real(kind=8) :: theta0,cost0,sint0,sl,rerf,c1,c2,p,zdim,ztop_dim,zbot_dim
! for JGR97 Appendix C
real(kind=8) :: psimax=250.0_8/(pi*pi)
parameter(theta0=pi*0.25,cost0=sqrt(2.0_8)/2.0_8,sint0=cost0)
! specifically for axisymmetric spherical shell benchmark of van Keken PEPI 01
real(kind=8),parameter :: theta_pepi01=63.5d0*pi/180.0_8,wavel_pepi01=53.0_8*pi/180.0_8
real(kind=8) :: ybottom=0.0_8,ytop=1.0_8,adiabatic_difference=0.0_8,funky=0.0_8,funky2=0.0_8
real(kind=8) :: delta,pert

if (ichoice > 100) then
   func = funccf(ichoice-100,x,y,z)
   return
endif

if (.not.cyl) then
   ! Cartesian geometry
   func=0.0_8
   if (ichoice==12) then
      ! conductive temperature profile 
      func = t_bot*(1 - y)
      return
   else if (ichoice==1) then
      ! harmonic perturbation scaled with temperature difference across the domain
      func =  0.1*(T_bot-T_top)*cos(pi*x/wavel_perturb)*sin(pi*y)
      return
   else if (ichoice==2) then
      ! conductive temperature profile + harmonic perturbation
      if (use_varRa) then
        func = T_bot*(1 - y + 0.2*cos(pi*x/wavel_perturb)*sin(pi*y))
      else
        func = T_bot*(1 - y + 0.1*cos(pi*x/wavel_perturb)*sin(pi*y))
      endif
      return
   else if (ichoice==3) then
      ! density 
      func=funccf(3,x,y,z)
   else if (ichoice==4) then
     ! horizontal velocity from Appendix C, Van Keken et al., 1997
     func = psimax * sin (pi*x) * cos (pi*y)
   else if (ichoice==5) then
     ! vertical velocity from Appendix C, Van Keken et al., 1997
     func = - psimax * cos (pi*x) * sin (pi*y)
   else if (ichoice==6) then
     ! horizontal velocity from Nate Sime 2020 rotational test
   else if (ichoice==7) then
     ! adiabatic temperature difference
     !if (Kequivalent) then
     if (solve_for_Tperturb) then
        func=Tbars_nondimK*exp(Di*(1-y))
     else
        func=0.0
     endif
     return
     !endif
     !zdim=(1-y)*height_dim
     !func = get_Tbar_dim(zdim) -get_Tbar_dim(dzero)! dimensional
     !if (cartplume) func=get_Tbar_dim(zdim)-T0_dim
     !if (.not.delta1K) func=func/deltaT_dim
!    if (x<1e-6) write(6,'(''func7: '',5e15.7)') y,zdim,get_Tbar_dim(zdim),get_Tbar_dim(dzero),func
   else if (ichoice==10) then
     func = x*x + y*y
     return
   else if (ichoice==11) then
     ! boundary layer with perturbation for plume-like models
     if (cartplume) then
        ! use fixed non-dimensional parameters for the test case
        if (x>pi/16) then
           gx=0.0_8
        else 
           gx=cos(8*x)
        endif
        gx=gx*ampini_perturb ! fixed except for amplitude...
        fxy=1.0_8-rerf(y/0.045)+gx*max(0.0_8,1.0_8-y/0.045)
        func=750*min(fxy,1.0_8)
       !if (x<1e-2.and.print_node) then
       !   write(irefwr,'(''x,y, T: '',4e15.7)') x,y,func,pert
       !endif
     else
       pert=0.0_8
       drh = y
       func= (1.0_8-rerf(drh/deltah_nondim))
       if (x <= pi/kwave_perturb) then
          pert=ampini_perturb*cos(pi*x/rlamtemp)*max(0.0_8,1-0.5*drh/deltah_nondim)
       endif
       func=min(1.0_8,func+pert)
       if (delta1K) func=func*deltaT_dim
       !if (x<1e-2.and.print_node) then
       !   write(irefwr,'(''x,y, T: '',4e15.7)') x,y,func,pert
       !endif
     endif
   else if (ichoice==20) then
     ! initial condition for JGR97 benchmark 1
     if (y<=0.2+0.02*cos(pi*x/rlampix)) then
        func=1.0
     else
        func=0.0
     endif
     return

  else if (ichoice.eq.100) then
     ! temperature field for benchmark 2 from Van Keken et al., 1997
     xr = x
     if (xr.le.0.00001d0) xr=0.00001d0
     if (xr.ge.1.99999d0) xr=1.99999d0
     wpi = sqrt(pi)
     zz   = y
     c1=1d0
     c2=2d0
     p=1d0
     ah=2d0
     rara=3d5
     u0 = c1*(ah**(7d0/3d0))*((rara*0.5d0/wpi)**(2d0/3d0))
     u0 = u0 / (( 1d0+ah*ah*ah*ah)**(2d0/3d0))
     v0 = c2*u0/ah
     qq  = 2d0*(ah/(pi*u0))**0.5d0
     trarg = -    xr  *   xr  *v0/(4d0*( zz+p ))
     tsarg = - (ah-xr)*(ah-xr)*v0/(4d0*(1d0-zz+p))
     ts = 0.5d0 - 0.5d0*qq/wpi*sqrt(v0/(1d0-zz+p))*exp(tsarg)
     tr = 0.5d0 + 0.5d0*qq/wpi*sqrt(v0/( zz+p ))*exp(trarg)
     tu =         0.5d0*rerf(0.5d0*(1d0-zz)*(u0/   xr  )**0.5d0)
     tl =   1d0 - 0.5d0*rerf(0.5d0*  zz  *(u0/(ah-xr))**0.5d0)
     func = tu+tl+tr+ts - 1.5d0
     if (func.gt.1d0) func=1d0
     if (func.lt.0d0) func=0d0
     if (y.eq.0) func=1d0
     if (y.eq.1) func=0d0
   else if (ichoice==99) then
     ! temperature field for benchmark 2 from Van Keken et al., 1997
     ! modified for compressible convection
     ! Potential T follows analytical expression; add Tbar
     xr = x
     if (xr.le.0.00001d0) xr=0.00001d0
     if (xr.ge.1.99999d0) xr=1.99999d0
     wpi = sqrt(pi)
     zz   = y
     c1=1d0
     c2=2d0
     p=1d0
     ah=2d0
     rara=3d5
     u0 = c1*(ah**(7d0/3d0))*((rara*0.5d0/wpi)**(2d0/3d0))
     u0 = u0 / (( 1d0+ah*ah*ah*ah)**(2d0/3d0))
     v0 = c2*u0/ah
     qq  = 2d0*(ah/(pi*u0))**0.5d0
     trarg = -    xr  *   xr  *v0/(4d0*( zz+p ))
     tsarg = - (ah-xr)*(ah-xr)*v0/(4d0*(1d0-zz+p))
     ts = 0.5d0 - 0.5d0*qq/wpi*sqrt(v0/(1d0-zz+p))*exp(tsarg)
     tr = 0.5d0 + 0.5d0*qq/wpi*sqrt(v0/( zz+p ))*exp(trarg)
     tu =         0.5d0*rerf(0.5d0*(1d0-zz)*(u0/   xr  )**0.5d0)
     tl =   1d0 - 0.5d0*rerf(0.5d0*  zz  *(u0/(ah-xr))**0.5d0)
     func = tu+tl+tr+ts - 1.5d0

     ! truncate between 0 and 1
     if (func.gt.1.0_8) func=1.0_8
     if (func.lt.0d0) func=0.0_8
     if (y<1e-7_8) func=1.0_8
     if (y>0.999999999_8) func=0.0_8

     ! scale total amplitude by difference in Tbar(z=0) and Tbar(z=1)
     ybottom=0.0_8
     ytop=1.0_8
     zbot_dim=(1-ybottom)*height_dim
     ztop_dim=(1-ybottom)*height_dim
     adiabatic_difference=get_Tbar_dim(zbot_dim)-get_Tbar_dim(ztop_dim)
     funky2=func
     func=func*(1.0_8-adiabatic_difference)
     funky=func
     ! truncate between 0 and Tpot_max
     if (func.gt.1.0_8) func=1.0_8-adiabatic_difference
     if (func.lt.0d0) func=0.0_8
     if (y<1e-7_8) func=1.0_8-adiabatic_difference
     if (y>0.999999999_8) func=0.0_8
   else
      if (print_node) write(irefwr,*) 'PERROR(func): ichoice /= 1,2,3,4,5,7,10,100: ',ichoice
      call instop
   endif


else 

   ! cylinder geometry
   if (ichoice==1) then

!    thermal perturbation

     if (ibench_type == 2 .and. axi) then

!       initial condition specifically for axisymmetric
!       spherical shell model of van Keken, PEPI, 2001
!       rotate by 42 degrees, wavelength of 58 degrees
        r = sqrt(x*x+y*y)
        drh = (r-R1)/(R2-R1)
        theta = asin(y/r)
        th = (theta+theta_pepi01)/wavel_pepi01
        argsin = pi*(1-drh)
        func = 0.10*sin(argsin)*sin(th)
        if (delta1K) func=func*deltaT_dim

     else

        ntheta=4
        if (ibench_type==2) ntheta=2
        r = sqrt(x*x+y*y)
        drh = (r-R1)/(R2-R1)
        theta = asin(y/r)
        if (x<0) theta=2*pi-theta
        argsin = pi*(1-drh)
        !func = 0.1*sin(argsin)*cos((theta+0.5*pi)*ntheta)
        func=0.1*sin(argsin)*cos(theta*ntheta)
        if (r==r1 .or. r==r2) then
          if (print_node) write(irefwr,'(''func_cyl1: '',6f12.3)') x,y,drh,theta,th,func
        endif
        if (delta1K) func=func*deltaT_dim

     endif


   else if (ichoice==2) then
      ! linear profile with perturbation
        r = sqrt(x*x+y*y)
        drh = (r-R1)/(R2-R1)
        theta = asin(y/r)
        argsin = pi*(1-drh)
        func = 1-drh + 0.1*sin(argsin)*cos(theta*ntheta)
        if (r==r1 .or. r==r2) then
          if (print_node) write(irefwr,'(''func_cyl2: '',5f12.3)') x,y,drh,theta,func
        endif
        if (delta1K) func=func*deltaT_dim
        return
  else if (ichoice.eq.3) then

!    density 
     func = funccf(3,x,y,z)

  else if (ichoice.eq.4) then

!     horizontal velocity from Appendix C of Van Keken et al., JGR, 1997
      if (r1<0.5) then
         xp = x - 0.4192
      else
         xp = x - r1
      endif
      yp = y + 0.5
      func = psimax * sin (pi*xp) * cos (pi*yp)
      if (quart) then
!        rotate velocity field 45 degrees clockwise
         xp =  sint0*x + cost0*y
         yp = -cost0*x + sint0*y
         xp = xp - r1
         yp = yp + 0.5
         vx = psimax * sin (pi*xp) * cos (pi*yp)
         vy = - psimax * cos (pi*xp) * sin (pi*yp)
         func = vx*sint0-vy*cost0
      endif

   else if (ichoice.eq.5) then

!     vertical velocity from Appendix C of Van Keken et al., JGR, 1997
      if (r1<0.5) then
         xp = x - 0.4192
      else
         xp = x - r1
      endif
      yp = y + 0.5
      func = - psimax * cos (pi*xp) * sin (pi*yp)
      if (quart) then
!        rotate velocity field 45 degrees clockwise
         xp =  sint0*x + cost0*y
         yp = -cost0*x + sint0*y
         xp = xp - r1
         yp = yp  + 0.5
         vx = psimax * sin (pi*xp) * cos (pi*yp)
         vy = - psimax * cos (pi*xp) * sin (pi*yp)
         func = vx*cost0 + vy*sint0
      endif
   else if (ichoice == 7) then
!     adiabatic temperature difference
      r = sqrt(x*x+y*y)
      zdim=(R2-r)*height_dim
      func = get_Tbar_dim(zdim)-get_Tbar_dim(dzero) ! dimensional
     if (.not.delta1K) func=func/deltaT_dim
      return


   else if (ichoice.eq.10) then
     func = x*x + y*y
     return

   else
     if (print_node) write(irefwr,*) 'PERROR(func) :: ichoice/=1,2,3,4,5,10'
     call instop
   endif
endif


end function func

! Compute error function (from Numerical Recipes)
real(kind=8) function rerf(x)
implicit none
real(kind=8) :: x
integer :: iret
real(kind=8) :: delx,x0,x1,xx

if (x.gt.4d0) then
      rerf=1d0
      return
endif
   iret=0
   rerf=0.
   delx=0.0001d0
   x0=0.
10       x1=x0+delx
   if(x1.gt.x) then
      x1=x
      iret=1
   endif
   xx=(x0+x1)/2.
   rerf=rerf+dexp(-xx*xx)*(x1-x0)
   if(iret.eq.1) then
      rerf=rerf*2d0/dsqrt(3.14159265358979323846_8)
      return
   endif
   x0=x1
   goto 10
end function rerf
