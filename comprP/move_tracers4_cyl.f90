subroutine move_tracers4_cyl(ichoice,isolold,isolnew,tstep,tfach)
use sepmodulecomio
use sepmoduleoldrouts
use sepmodulekmesh
use mpetrac
use tracers
use mparallel
use geometry
use control
use sepran_interface
use sepran_arrays
implicit none
integer :: ichoice,isolold,isolnew
real(kind=8) :: tstep,tfach,rtest

real(kind=8) :: dt_little,p,q,frac1,frac2,factor,xm,ym,u,v,xi,eta
real(kind=8) :: rk1x,rk1y,rk2x,rk2y,rk3x,rk3y,rk4x,rk4y
real(kind=8) :: vth,vr,theta_init,r_init,thetal,costh,sinth
real(kind=8) :: rk1th,rk2th,rk3th,rk4th,rk1r,rk2r,rk3r,rk4r,r
real(kind=8) :: dt_tiny,xm0,ym0,oneminfrac1,oneminfrac2

integer :: isolav1,isolav2,ntimes,ipoint,ntot,nmark
integer :: itimes,i,itrac,k,ip,ip1,ip2,ip_newvel_start,CHUNK,ielh=0,called=0
real(kind=4) :: t01


logical :: fliprk1,fliprk2,fliprk3,flipF
save isolav1,isolav2,called


needed_neighbors=0
needed_neighbors_neighbors=0
missed=0
! Prepare for interpolation of velocity (see tdetvel.f)
if (ichoice.eq.1) then
!  Base interpolated on constant velocity stored in isolold
   if (10+2*npoint.gt.user_here(1)) then
      if (print_node) then
         write(irefwr,*) 'PERROR(move_tracers4_cyl): user_here is too small'
         write(irefwr,*) 'user_here(1): ',user_here(1)
         write(irefwr,*) 'needed: ',10+2*npoint
      endif
      call instop
   endif
else
   if (10+4*npoint.gt.user_here(1)) then
      if (print_node) then
         write(irefwr,*) 'PERROR(move_tracers4_cart): user_here is too small'
         write(irefwr,*) 'user_here(1): ',user_here(1)
         write(irefwr,*) 'needed: ',10+4*npoint
      endif
      call instop
   endif
endif

ntimes = nint(2*tfach)
dt_little = tstep/ntimes

if (itracoption.eq.2) then
!  ntot = total number of markers
   nmark=0
   do ichain=1,nochain
      nmark=nmark+imark(ichain)
   enddo
   ntot = nmark
else
!  ntot = total number of tracers
   ntot=0
   do idist=1,ndist
      ntot=ntot+ntrac(idist)
   enddo
endif
!called=called+1

factor=1d0
xi_eta_stored = .false.
! find velocity in tracers at initial time
call pecopy(2,user_here(10:10+npoint-1),isolold)
call pecopy(3,user_here(10+npoint:10+npoint*2-1),isolold)
CHUNK=ntot/4
do itrac=1,ntot
   tracer_number=itrac
   ip1 = 2*itrac-1
   ip2 = ip1+1
   xm = tracer(itrac)%x
   ym = tracer(itrac)%y
   tracer(itrac)%xhalf=xm
   tracer(itrac)%yhalf=ym

   call velintmark_cyl(0,xm,ym,u,v,ielh,xi,eta,'move_tracers4',itrac)
   tracer(itrac)%u = u
   tracer(itrac)%v = v
   ! write(irefwr,'(2i5,4f12.5)') called,0,xm,ym,u,v
enddo

!write(irefwr,'(''PINFO(move_tracers4): ntimes = '',i5,2e12.5,i5)') ntimes,tstep,dt_little,ichoice
do itimes=1,ntimes
   !write(irefwr,*) 'itimes = ',itimes

   if (ichoice==2) then
!     find velocity at time t=(i-0.5)*dt_little and t=i*dt_little
      frac1 = (itimes-0.5)*1d0/ntimes
      frac2 = itimes*1d0/ntimes
      oneminfrac1=1.0_8-frac1
      oneminfrac2=1.0_8-frac2
!     write(irefwr,*) 'frac1, frac2: ',frac1,frac2
      call algebr(3,0,isolold,isolnew,isolav1,kmesh1,kprob1,oneminfrac1,frac1,p,q,ipoint)
      call algebr(3,0,isolold,isolnew,isolav2,kmesh1,kprob1,oneminfrac2,frac2,p,q,ipoint)
!     copy new intermediate velocity
      call pecopy(2,user_here(10:10+npoint-1),isolav1)
      call pecopy(3,user_here(10+npoint:10+2*npoint-1),isolav1)
      call pecopy(2,user_here(10+2*npoint:10+3*npoint-1),isolav2)
      call pecopy(3,user_here(10+3*npoint:10+4*npoint-1),isolav2)
   endif

   ip_newvel_start = 1+(ichoice-1)*2*npoint

   CHUNK=ntot/4
   do itrac=1,ntot
!     find mid-point in current (sub)interval
      ip1=2*itrac-1
      ip2=ip1+1
      xm = tracer(itrac)%xhalf
      ym = tracer(itrac)%yhalf
      r = sqrt(xm*xm+ym*ym)
      if (r.lt.rtop_threshold.and.r.gt.rbot_threshold) then
          ! write(irefwr,*) 'cartesian tracing'
!         particle is sufficiently far away from top and
!         bottom boundary to use Cartesian coordinates
!         for particle tracing
          rk1x = dt_little*tracer(itrac)%u
          rk1y = dt_little*tracer(itrac)%v
          xm = tracer(itrac)%xhalf + 0.5*rk1x
          ym = tracer(itrac)%yhalf + 0.5*rk1y
          call checkbounds_cyl(0,xm,ym)

!         find average velocity at the midpoint
          call velintmark_cyl(0,xm,ym,u,v,ielh,xi,eta,'move_tracers4',itrac)
         ! write(irefwr,'(2i5,4f12.5)') called,1,xm,ym,u,v
          rk2x = dt_little*u
          rk2y = dt_little*v
          xm = tracer(itrac)%xhalf + 0.5*rk2x
          ym = tracer(itrac)%yhalf + 0.5*rk2y
          call checkbounds_cyl(0,xm,ym)

!         *** find average velocity at corrected mid-point
          call velintmark_cyl(0,xm,ym,u,v,ielh,xi,eta,'move_tracers4',itrac)
         ! write(irefwr,'(2i5,4f12.5)') called,2,xm,ym,u,v

!         *** find position of marker at t=t+dt
          rk3x = dt_little*u
          rk3y = dt_little*v
          xm = tracer(itrac)%xhalf + rk3x
          ym = tracer(itrac)%yhalf + rk3y
          call checkbounds_cyl(0,xm,ym)

          ! ioffset is 0 when ip_newvel_start is 1
          call velintmark_cyl(ip_newvel_start-1,xm,ym,u,v,ielh,xi,eta,'move_tracers4',itrac)
         ! write(irefwr,'(2i5,4f12.5)') called,3,xm,ym,u,v
          rk4x = dt_little*u
          rk4y = dt_little*v
          call checkbounds_cyl(0,xm,ym)
!         Combine all the estimates in the 4th order Runge Kutta step
          tracer(itrac)%xnew = tracer(itrac)%xhalf + rk1x/6d0 + rk2x/3d0 + rk3x/3d0 + rk4x/6d0
          tracer(itrac)%ynew = tracer(itrac)%yhalf + rk1y/6d0 + rk2y/3d0 + rk3y/3d0 + rk4y/6d0
          call checkbounds_cyl(0,tracer(itrac)%xnew,tracer(itrac)%ynew)
          !call checkbounds_cyl(0,xm,ym)

      else

!         Use cylindrical coordinates for tracing of particles
!         near the top surface
          ! write(irefwr,*) 'cylindrical tracing'
          xm0=xm
          ym0=ym
          thetal = acos(ym/r)
          if (xm.lt.0d0) thetal = 2.0d0*pi - acos(ym/r)

          theta_init = thetal
          r_init     = r

          u   = tracer(itrac)%u
          v   = tracer(itrac)%v
          costh = cos(thetal)
          sinth = sin(thetal)
          vth = costh*u - sinth*v
          vr  = sinth*u + costh*v

!         note: vth is not really the angular velocity yet
!         divide by r

          rk1th = dt_little*vth/r
          rk1r  = dt_little*vr
!         Calculate estimate for mid-point
          thetal = theta_init + 0.5*rk1th
          r     =  r_init     + 0.5*rk1r
          r = min(r,radius_max)
          r = max(r,radius_min)
!         convert to cartesian
          costh = cos(thetal)
          sinth = sin(thetal)
          xm = r*sinth
          ym = r*costh
          !write(irefwr,'(''rk1: '',4f15.9)') xm-xm0,ym-ym0,u,v

!         find new velocity in Cartesian coordinates
          call velintmark_cyl(0,xm,ym,u,v,ielh,xi,eta,'move_tracers4',itrac)
         ! write(irefwr,'(2i5,4f12.5)') called,4,xm,ym,u,v
!         convert to cylindrical
          costh = cos(thetal)
          sinth = sin(thetal)
          vth = costh*u - sinth*v
          vr  = sinth*u + costh*v
!         advance tracer to new position
          rk2th = dt_little*vth/r
          rk2r  = dt_little*vr
          thetal = theta_init + 0.5*rk2th
          r     =  r_init     + 0.5*rk2r
          r = min(r,radius_max)
          r = max(r,radius_min)
          costh = cos(thetal)
          sinth = sin(thetal)
          xm = r*sinth
          ym = r*costh
          !write(irefwr,'(''rk2: '',4f15.9)') xm-xm0,ym-ym0,u,v

          call velintmark_cyl(0,xm,ym,u,v,ielh,xi,eta,'move_tracers4',itrac)
         ! write(irefwr,'(2i5,4f12.5)') called,5,xm,ym,u,v
!         Find position of marker at t=t+dt
!         convert to cylindrical
          vth = costh*u - sinth*v
          vr  = sinth*u + costh*v
          rk3th = dt_little*vth/r
          rk3r  = dt_little*vr
          thetal = theta_init + rk3th
          r     =  r_init     + rk3r
          r = min(r,radius_max)
          r = max(r,radius_min)
!         convert new position to Cartesian
          costh = cos(thetal)
          sinth = sin(thetal)
          xm = r*sinth
          ym = r*costh
          !write(irefwr,'(''rk3: '',4f15.9)') xm-xm0,ym-ym0,u,v

          call velintmark_cyl(ip_newvel_start-1,xm,ym,u,v,ielh,xi,eta,'move_tracers4',itrac)
         ! write(irefwr,'(2i5,4f12.5)') called,6,xm,ym,u,v
!         convert to cylindrical
          vth = costh*u - sinth*v
          vr  = sinth*u + costh*v
          rk4th = dt_little*vth/r
          rk4r  = dt_little*vr

!         Combine all the estimates in the 4th order Runge Kutta step
          thetal  = theta_init + rk1th/6d0 + rk2th/3d0 + rk3th/3d0 + rk4th/6d0
          r   =  r_init + rk1r/6d0 + rk2r/3d0 + rk3r/3d0 + rk4r/6d0
          r = min(r,radius_max)
          r = max(r,radius_min)
!         final coordinates in Cartesian
          xm = r*sin(thetal)
          ym = r*cos(thetal)
          !write(irefwr,'(''all: '',4f15.9)') xm-xm0,ym-ym0,u,v
          !call checkbounds_cyl(1,xm,ym)
 
          tracer(itrac)%xnew = xm
          tracer(itrac)%ynew = ym
        endif ! end of (r,th) tracing
           
   enddo ! loop over tracers

!  markers have been updated, information in tracer%ielem is out of date
   xi_eta_stored = .false.


!  *** and, in case of the markerchain ...
   if (itracoption.eq.2) then
!      make sure markers at boundaries stay at boundary
       ip=0
       do ichain=1,nochain
          nmark=imark(ichain)
          ip1 = ip+1
          tracer(ip1)%xnew = xcmax
          ip2 = ip+nmark
          tracer(ip2)%xnew=xcmin
          ip=ip+nmark
       enddo
   endif

   if (itimes.lt.ntimes) then
      ! set up initial condition for new sub-iteration
      ! interpolate velocity in new points into tracer%(u,v)
      do itrac=1,ntot
         ip1 = 2*itrac-1
         ip2 = ip1+1
         xm = tracer(itrac)%xnew
         ym = tracer(itrac)%ynew

         call velintmark_cyl(ip_newvel_start-1,xm,ym,u,v,ielh,xi,eta,'move_tracers4',itrac)
         ! write(irefwr,'(2i5,4f12.5)') called,7,xm,ym,u,v
         tracer(itrac)%u = u
         tracer(itrac)%v = v
         tracer(itrac)%xhalf = xm
         tracer(itrac)%yhalf = ym
      enddo
   endif

enddo ! loop over itimes

end subroutine  move_tracers4_cyl
