!   PREDCOORT
!
!   Explicit Euler step for tracer/markerchain advance during
!   the predictor step.
!
!   tracer%(xnew,ynew):= tracer%(x,y) + tracer%(u,v)*tstep
!
!   Information about the number of tracers/markers is stored
!   in tracer_mod
!   
!   This version combines the cylindrical and cartesian box versions
subroutine predcoort()
use geometry
implicit none

if (cyl) then
   call predcoortcyl()
else
   call predcoortcart()
endif

return
end subroutine predcoort
  

!   PREDCOORTCYL
!
!   predcoor, adapted for use with tracers 
!
!   PvK 960214
!
!   PvK 970411: modified for cylindrical geometry.
!   PvK 073104: (re)adapted for markerchain
subroutine predcoortcyl()
use mtime
use tracers
use coeff
use mpetrac
use geometry
implicit none
!include 'SPcommon/ctimen'
integer :: itrac,j,ntot,nmark,ip
real(kind=8) :: xm,ym,r,th,cost
real(kind=8) :: thetal,costh,sinth,vr,vth,u,v
integer :: ip1,ip2

     
if (itracoption == 2) then
!  markerchain method
   nmark=0
   do ichain=1,nochain
      nmark = nmark+imark(ichain) 
   enddo
   ntot=nmark
else  ! includes itracoption=3 
   ntot=0
   do idist=1,ndist
      ntot=ntot+ntrac(idist)
   enddo
endif

do itrac = 1,ntot
   xm = tracer(itrac)%x
   ym = tracer(itrac)%y
   r = sqrt(xm*xm+ym*ym)

   if (r.lt.rtop_threshold.and.r.gt.rbot_threshold) then 
!     particle is sufficiently far away from the
!     top and bottom boundary to use Cartesian coordinates

      xm = xm+tstepp*tracer(itrac)%u
      ym = ym+tstepp*tracer(itrac)%v

   else

!     particle is deemed too close to top or bottom boundary;
!     advect particle in cylindrical coordinates such that
!     the rigid rotation of the surface is taken into
!     account in an accurate manner.
      thetal = acos(ym/r)
      u = tracer(itrac)%u
      v = tracer(itrac)%v
      costh = cos(thetal)
      sinth = sin(thetal)
      vth = costh*u - sinth*v
      vr  = sinth*u + costh*v
      thetal = thetal + tstepp*vth/r
      r      = r      + tstepp*vr
      r = min(r,radius_max)
      r = max(r,radius_min)
      xm = r*sin(thetal)
      ym = r*cos(thetal)

   endif

   call checkbounds(0,xm,ym)
   if (itracoption == 2) call stay_on_boundary()
   if (itracoption == 3) xm=0
   tracer(itrac)%xnew = xm
   tracer(itrac)%ynew = ym
enddo

if (itracoption == 2) then
!  make sure markers at boundaries stay on boundary
   if (quart) then
      ip=0
      do ichain=1,nochain
         nmark = imark(ichain)
!        x coordinate of first tracer should be zero
         ip1 = ip+1
         tracer(ip1)%xnew = 0
!        y coordinate of last tracer should be zero
         ip2 = ip+nmark
         tracer(ip2)%ynew = 0
         ip = ip+nmark
      enddo
   else if (half) then
      ip=0
      do ichain=1,nochain
         nmark = imark(ichain)
!        x coordinate of first tracer should be zero
         ip1 = ip+1
         tracer(ip1)%xnew = 0
!        x coordinate of last tracer should be zero
         ip2 = ip+nmark
         tracer(ip2)%xnew = 0
         ip = ip+nmark
      enddo
   else if (eighth) then
      ip=0
      do ichain=1,nochain
         nmark = imark(ichain)
!        x coordinate of first tracer should be zero
         ip1 = ip+1
         tracer(ip1)%xnew = 0
!        x and y coordinates of last tracer should be equal
         ip2 = ip+nmark
         tracer(ip2)%xnew = tracer(ip2)%ynew
         ip = ip+nmark
      enddo
   endif
endif

end subroutine predcoortcyl

!   PREDCOORTCART
!
!   PvK 10-8-89
subroutine predcoortcart()
use mtime
use tracers
use mpetrac
use coeff
use geometry
implicit none

integer :: ntot,nmark,itrac,j,ip1,ip2,ip

    
if (itracoption == 2) then
   nmark=0
   do ichain=1,nochain
      nmark = nmark + imark(ichain)
   enddo
   ntot = nmark
else ! includes itracoption=3
   ntot = 0
   do idist=1,ndist
      ntot = ntot+ntrac(idist)
   enddo
endif
do itrac = 1,ntot
   tracer(itrac)%xnew = tracer(itrac)%x + tstepp*tracer(itrac)%u
   tracer(itrac)%ynew = tracer(itrac)%y + tstepp*tracer(itrac)%v
   call checkbounds_cart(tracer(itrac)%xnew,tracer(itrac)%ynew)
enddo

if (itracoption == 2) then
!  make sure markers at boundaries stay at boundary
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

end subroutine predcoortcart
