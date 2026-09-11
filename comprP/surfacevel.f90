subroutine surfacevel(kmesh,kprob,isol,veloc,veloc_d,noutje)
use sepmodulecomio
use convparam
use coeff
use mparallel
use geometry
use control
implicit none
integer :: kmesh,kprob,isol,noutje
real(kind=8) :: veloc(:),veloc_d(:)
integer,parameter :: NCMAX=10000
real(kind=8),dimension(NCMAX) :: funcx,funcy
integer :: ichoice,icurvs(10),ncoord,npoint,i
real(kind=8) :: velmax,avvel,x,y,vel,cost,r,theta,vx,vy,velint2
real(kind=8) :: velint,xold,yold,vold,dx,dy,dist,circ,rmidv,funcbc
character(len=80) :: fname

funcx(1)=1.0_8*NCMAX
funcy(1)=1.0_8*NCMAX
ichoice=0
icurvs(1)=0
icurvs(2)=itop

! note that we compute the tangential velocity component here
call compcr(ichoice,kmesh,kprob,isol,-2,icurvs,funcx,funcy)
ncoord=nint(funcx(5)/2)
npoint=nint(funcy(5))
if (ncoord/=npoint) then
   if (print_node) then
      write(irefwr,*) 'something weird in surfacevel'
      write(irefwr,*) 'ncoord, npoint: ',ncoord,npoint
   endif
   call instop
endif
velmax=0
avvel=0
do i=1,npoint
   vel = funcy(5+i)
   velmax = max(velmax,abs(vel))
   avvel = avvel+abs(vel)
enddo
veloc(1) = velmax

if (noutje>=0) then
   if (print_node) then
     write(fname,'(''GMT/surfacevel.'',i4.4)') noutje
     open(9,file=fname)
     if (cyl) then
        ! plot tangential velocity vs. azimuth (clockwise from
        ! north)
        do i=1,npoint
           x=funcx(4+2*i)
           y=funcx(5+2*i)
           r=sqrt(x*x+y*y)
           cost=y/r
           theta=acos(cost)
           if (x<0.0_8) theta=2*pi-theta
           write(9,*) theta*180.0/pi,funcy(5+i)
        enddo
     else if (.not.axi) then
!       if (its_CH94) then
!         do i=1,npoint
!            write(9,*) funcx(4+2*i),-funcy(5+i),funcbc(10,funcx(4+2*i),0.0_8,0.0_8)
!            ! negative because it's u_t and curve 3 goes to the left
!         enddo
!       else
          do i=1,npoint
             write(9,*) funcx(4+2*i),-funcy(5+i)  ! negative because it's u_t and curve 3 goes to the left
          enddo
!       endif
     endif
     close(9)
  endif ! print_node
endif ! noutje>=0


! compute average velocity
! cylindrical version
xold=funcx(6)
yold=funcx(7)
vold=abs(funcy(6))
velint=0
velint2=0
avvel=0
circ=0
do i=2,npoint
   x=funcx(4+2*i)
   y=funcx(5+2*i)
   vel = abs(funcy(5+i))
   dx = x-xold
   dy = y-yold
   rmidv = sqrt((x+0.5*dx)*(x+0.5*dx)+(y+0.5*dy)*(y+0.5*dy)) 
   if (icoor900.eq.0) then
!     If cylindrical use simple mean
      dist = sqrt(dx*dx+dy*dy)
      circ = circ+dist
      velint = velint+0.5*(vold+vel)*dist
      velint2 = velint2 + 0.5*(vold*vold+vel*vel)*dist
   else if (icoor900.eq.1) then
!     if axisymmetric, take the 2*pi*r factor into account
      dist = sqrt(dx*dx+dy*dy)
      circ = circ+dist*2*pi*rmidv
      velint = velint+0.5*(vold+vel)*dist*2*pi*rmidv
   endif
   yold=y
   xold=x
   vold=vel
enddo
avvel = velint/circ
veloc(2) = avvel
veloc(3) = sqrt(velint2/circ)

! dimensionalize velocities (cm/yr)
veloc_d(1) = veloc(1)*rkappa_dim/(R2_dim-R1_dim)*100*year_dim
veloc_d(2) = veloc(2)*rkappa_dim/(R2_dim-R1_dim)*100*year_dim

! note that we compute the tangential velocity component here
icurvs(1)=0
icurvs(2)=ibottom
call compcr(ichoice,kmesh,kprob,isol,-2,icurvs,funcx,funcy)
ncoord = nint(funcx(5)/2)
npoint = nint(funcy(5))
if (ncoord.ne.npoint) then
   if (print_node) then
      write(irefwr,*) 'something weird in surfacevel' 
      write(irefwr,*) 'ncoord, npoint: ',ncoord,npoint
   endif
   call instop
endif

velmax=0
avvel=0
do i=1,npoint
   vel = funcy(5+i)
   velmax = max(velmax,abs(vel))
   avvel = avvel+abs(vel)
enddo
veloc(4) = velmax
veloc(5) = avvel 

xold=funcx(6)
yold=funcx(7)
vold=abs(funcy(6))
velint=0
velint2=0
avvel=0
circ=0
velmax=0
avvel=0
do i=2,npoint
   x=funcx(4+2*i)
   y=funcx(5+2*i)
   vel = abs(funcy(5+i))
   dx = x-xold
   dy = y-yold
   rmidv = sqrt((x+0.5*dx)*(x+0.5*dx)+(y+0.5*dy)*(y+0.5*dy)) 
   if (icoor900.eq.0) then
!     If cylindrical use simple mean
      dist = sqrt(dx*dx+dy*dy)
      circ = circ+dist
      velint = velint+0.5*(vold+vel)*dist
      velint2 = velint2 + 0.5*(vold*vold+vel*vel)*dist
   else if (icoor900.eq.1) then
!     if axisymmetric, take the 2*pi*r factor into account
      dist = sqrt(dx*dx+dy*dy)
      circ = circ+dist*2*pi*rmidv
      velint = velint+0.5*(vold+vel)*dist*2*pi*rmidv
   endif
   yold=y
   xold=x
   vold=vel
enddo
avvel = velint/circ
veloc(5) = avvel
veloc(6) = sqrt(velint2/circ)

! dimensionalize velocities (cm/yr)
veloc_d(1:6) = veloc(1:6)*rkappa_dim/(R2_dim-R1_dim)*100*year_dim

end subroutine surfacevel
