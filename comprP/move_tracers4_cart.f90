subroutine move_tracers4_cart(ichoice,isolold,isolnew,tstep,tfach)
use sepmoduleoldrouts
use sepmodulekmesh
use mpetrac
use tracers
use mparallel
use geometry
use control
use sepran_arrays
use sepran_interface
implicit none
integer :: ichoice,isolold,isolnew
real(kind=8) :: tstep,tfach

real(kind=8) :: dt_little,p,q,frac1,frac2,factor,xm,ym,u,v,xi,eta
real(kind=8) :: rk1x,rk1y,rk2x,rk2y,rk3x,rk3y,rk4x,rk4y,oneminfrac1,oneminfrac2
integer :: isolav1=0,isolav2=0,ntimes,ipoint,ntot,nmark
integer :: itimes,i,j,k,ip,ip1,ip2,ip_newvel_start,ielh,itrac
save isolav1,isolav2

! Prepare for interpolation of velocity (see tdetvel.f)
call pefilxy(2)
!if (ichoice.eq.1) then
!!  Base interpolated on constant velocity stored in isolold
!   if (10+2*npoint.gt.user(1)) then
!      write(irefwr,*) 'PERROR(move_tracers4_cart): user is too small'
!      write(irefwr,*) 'user(1): ',user(1)
!      write(irefwr,*) 'needed: ',10+2*npoint
!      call instop
!   endif
!else 
!   if (10+4*npoint.gt.user(1)) then
!      write(irefwr,*) 'PERROR(move_tracers4_cart): user is too small'
!      write(irefwr,*) 'user(1): ',user(1)
!      write(irefwr,*) 'needed: ',10+2*npoint
!      call instop
!   endif
!endif

ntimes = int(2*tfach)
dt_little = tstep/ntimes
!write(irefwr,*) 'move_tracers4: tstep = ',tstep,dt_little,tfach,ntimes


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

factor=1d0
xi_eta_stored =.false.
! find velocity at initial time
call pecopy(2,user_here(10:10+npoint-1),isolold) ! isolold,user,kmesh,kprob,10,factor)
call pecopy(3,user_here(10+npoint:10+npoint*2-1),isolold) ! isolold,user,kmesh,kprob,10+npoint,factor)
do itrac=1,ntot
   ip1 = 2*itrac-1
   ip2 = ip1+1
   xm = tracer(itrac)%x
   ym = tracer(itrac)%y
   tracer(itrac)%xhalf=xm
   tracer(itrac)%yhalf=ym
   call velintmark_cart(xm,ym,u,v,ielh,xi,eta,user_here,'move_tracers4_cart')
   tracer(itrac)%u = u
   tracer(itrac)%v = v
enddo

do itimes=1,ntimes

   if (ichoice.eq.2) then
!     find velocity at time t=(i-0.5)*dt_little and t=i*dt_little
      frac1 = (itimes-0.5)*1d0/ntimes
      frac2 = itimes*1d0/ntimes
      oneminfrac1=1.0_8-frac1
      oneminfrac2=1.0_8-frac2
      call algebr(3,0,isolold,isolnew,isolav1,kmesh1,kprob1,oneminfrac1,frac1,p,q,ipoint)
      call algebr(3,0,isolold,isolnew,isolav2,kmesh1,kprob1,oneminfrac2,frac2,p,q,ipoint)
!     copy new intermediate velocity
      call pecopy(2,user_here(10:10+npoint-1),isolav1) ! ,user,kmesh,kprob,10,factor)
      call pecopy(3,user_here(10+npoint:10+2*npoint-1),isolav1) ! ,user,kmesh,kprob,10+npoint,factor)
      call pecopy(2,user_here(10+2*npoint:10+3*npoint-1),isolav2) ! ,user,kmesh,kprob,10+2*npoint,factor)
      call pecopy(3,user_here(10+3*npoint:10+4*npoint-1),isolav2) ! ,user,kmesh,kprob,10+3*npoint,factor)
   endif


   do itrac=1,ntot
!     find mid-point in current (sub)interval
      ip1=2*itrac-1
      ip2=ip1+1
      !if (itrac==1) write(6,*) 'movetracers 1'
      rk1x = dt_little*tracer(itrac)%u
      rk1y = dt_little*tracer(itrac)%v
      xm = tracer(itrac)%xhalf + 0.5*rk1x
      ym = tracer(itrac)%yhalf + 0.5*rk1y
      !write(6,*) 'checkbounds: ',xm,ym
      call checkbounds_cart(xm,ym)
      !write(6,*) 'velintmark_cart ',ielh
!     find average velocity at the midpoint 
      call velintmark_cart(xm,ym,u,v,ielh,xi,eta,user_here,'move_tracers4_cart')
!     if (j==1) write(irefwr,'(''move1: '',4f15.7)') xm,ym,u,v
 
      !if (itrac==1) write(6,*) 'movetracers 2'
      rk2x = dt_little*u
      rk2y = dt_little*v
      xm = tracer(itrac)%xhalf+0.5*rk2x
      ym = tracer(itrac)%yhalf+0.5*rk2y
      call checkbounds_cart(xm,ym)

!     find average velocity at corrected mid-point 
      call velintmark_cart(xm,ym,u,v,ielh,xi,eta,user_here,'move_tracers4_cart')
!     if (j==1) write(irefwr,'(''move1: '',4f15.7)') xm,ym,u,v

      !if (itrac==1) write(6,*) 'movetracers 3'
!     find position of marker at t=t+dt
      rk3x = dt_little*u
      rk3y = dt_little*v
      xm = tracer(itrac)%xhalf+rk3x
      ym = tracer(itrac)%yhalf+rk3y
      call checkbounds_cart(xm,ym)
      ip_newvel_start = 1+(ichoice-1)*2*npoint
      call velintmark_cart(xm,ym,u,v,ielh,xi,eta,user_here(ip_newvel_start:ip_newvel_start+npoint-1),'move_tracers4_cart')
!     if (j==1) write(irefwr,'(''move1: '',4f15.7)') xm,ym,u,v
      !if (itrac==1) write(6,*) 'movetracers 4'
      rk4x = dt_little*u
      rk4y = dt_little*v
      tracer(itrac)%xnew = tracer(itrac)%xhalf + rk1x/6d0 + rk2x/3d0 + rk3x/3d0 + rk4x/6d0
      tracer(itrac)%ynew = tracer(itrac)%yhalf+ rk1y/6d0 + rk2y/3d0 + rk3y/3d0 + rk4y/6d0
!     Make sure the markers stay in the computational domain
      call checkbounds_cart(tracer(itrac)%xnew,tracer(itrac)%ynew)
   enddo

   if (itracoption.eq.2) then
!      ... make sure markers at boundaries stay at boundary
       ip=0
       do ichain=1,nochain
          nmark=imark(ichain)
          ip1 = ip+1
          tracer(ip1)%xnew = xcmax
          ip2 = ip+nmark
          tracer(ip2)%xnew = xcmin
          ip=ip+nmark
       enddo
   endif

   if (itimes.lt.ntimes) then
!     set up initial condition for new sub-iteration
!     interpolate velocity in new points into tracer%(u,v)
      do itrac=1,ntot
         ip1 = 2*itrac-1
         ip2 = ip1+1
         xm = tracer(itrac)%xnew
         ym = tracer(itrac)%ynew
         call velintmark_cart(xm,ym,u,v,ielh,xi,eta,user_here(ip_newvel_start:ip_newvel_start+npoint-1),'move_tracers4_cart')
         tracer(itrac)%u = u
         tracer(itrac)%v = v
         tracer(itrac)%xhalf = xm
         tracer(itrac)%yhalf = ym
      enddo
   endif

enddo

end  subroutine move_tracers4_cart
