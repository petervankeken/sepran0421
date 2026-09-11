! CHECKBOUNDS
! 
! Make sure the tracers don't get advected past the mesh boundaries.
! PvK 2000
subroutine checkbounds(ichoice,xm,ym)
use control
use coeff
use tracers
use geometry
implicit none
integer :: ichoice
real(kind=8) :: xm,ym,r,cost,sint,th

if (ichoice.eq.0) then 
  r = sqrt(xm*xm+ym*ym)
  if (r.gt.radius_max+2e-7) then
     cost = ym/r
     if (xm.gt.0) then
        th = acos(cost)
     else
        th = 2*pi - acos(cost)
     endif
     ! if (print_node) write(irefwr,'(''modify 1: '',4f10.7,$)') xm,ym,r,th
     xm = xm - (r-radius_max+eps_tracer_bound(1))*sin(th)
     ym = ym - (r-radius_max+eps_tracer_bound(1))*cos(th)
     r = sqrt(xm*xm+ym*ym)
     ! if (print_node) write(irefwr,'(3x,3f10.7)') xm,ym,r
  endif

  if (r.lt.(radius_min-2e-7)) then
     cost = ym/r
     if (xm.gt.0) then
        th = acos(cost)
     else
        th = 2*pi - acos(cost)
     endif
     ! if (print_node) write(irefwr,'(''modify 2: '',4f12.9,$)') xm,ym,r,radius_min
     xm = xm + (radius_min-r+eps_tracer_bound(1))*sin(th)
     ym = ym + (radius_min-r+eps_tracer_bound(1))*cos(th)
     r = sqrt(xm*xm+ym*ym)
     ! if (print_node) write(irefwr,'(3f12.9)') r,xm,ym
  
  endif
endif

! for the half geometry
if (half) then    
   if (xm.lt.0) then
      xm=eps_tracer_bound(1)
   endif
else if (quart) then
!  and for the quart ...
   if (xm.lt.0) xm=eps_tracer_bound(1)
   if (ym.lt.0) ym=eps_tracer_bound(1)
else if (eighth) then
!  for an eight section
   if (xm.lt.0) xm=eps_tracer_bound(1)
   if (xm.ge.ym) xm = ym-eps_tracer_bound(1)
endif

end subroutine checkbounds

subroutine checkbounds_cyl(ichoice,xm,ym)
use tracers
use coeff
use geometry
implicit none
integer :: ichoice
real(kind=8) :: xm,ym,r,cost,sint,th

r = sqrt(xm*xm+ym*ym)
if (r.gt.radius_max+2e-7) then
     cost = ym/r
     if (xm.gt.0) then
        th = acos(cost)
     else
        th = 2*pi - acos(cost)
     endif
     xm = xm - (r-radius_max+eps_tracer_bound(1))*sin(th)
     ym = ym - (r-radius_max+eps_tracer_bound(1))*cos(th)
     r = sqrt(xm*xm+ym*ym)
endif

if (r.lt.(radius_min-2e-7)) then
     cost = ym/r
     if (xm.gt.0) then
        th = acos(cost)
     else
        th = 2*pi - acos(cost)
     endif
     xm = xm + (radius_min-r+eps_tracer_bound(1))*sin(th)
     ym = ym + (radius_min-r+eps_tracer_bound(1))*cos(th)
     r = sqrt(xm*xm+ym*ym)
endif

! for the half geometry
if (half) then    
   if (xm.lt.(2.0e-7)) then
         xm = eps_tracer_bound(1)
   endif
else if (quart) then
!  and for the quart ...
   if (xm.lt.0) xm=eps_tracer_bound(1)
   if (ym.lt.0) ym=eps_tracer_bound(1)
else if (eighth) then
!  for an eight section
   if (xm.lt.0) xm=eps_tracer_bound(1)
   if (xm.ge.ym) xm = ym-eps_tracer_bound(1)
endif

return
end subroutine checkbounds_cyl

subroutine checkbounds_cylp(ichoice,xm,ym,flip)
use tracers
use coeff
use geometry
implicit none
integer ichoice
real*8 xm,ym,r,cost,sint,th

logical flip

  r = sqrt(xm*xm+ym*ym)
  if (r.gt.radius_max+2e-7) then
     cost = ym/r
     if (xm.gt.0) then
        th = acos(cost)
     else
        th = 2*pi - acos(cost)
     endif
     xm = xm - (r-radius_max+eps_tracer_bound(1))*sin(th)
     ym = ym - (r-radius_max+eps_tracer_bound(1))*cos(th)
     r = sqrt(xm*xm+ym*ym)
  endif

  if (r.lt.(radius_min-2e-7)) then
     cost = ym/r
     if (xm.gt.0) then
        th = acos(cost)
     else
        th = 2*pi - acos(cost)
     endif
     xm = xm + (radius_min-r+eps_tracer_bound(1))*sin(th)
     ym = ym + (radius_min-r+eps_tracer_bound(1))*cos(th)
     r = sqrt(xm*xm+ym*ym)
  
  endif

! for the half geometry
if (half) then    
   if (periodic) then
      if (xm.lt.(0.0)) then
!         mirror the particle across the periodic b.c.
!         This is a really crude way to handle it, and leads to 
!               very inaccurate tracing near the boundary
!               xm = -1.0d0*xm
!               ym = -1.0d0*ym
          flip = .true.
      else
          flip = .false.
      endif

   else if (.not.periodic) then
       if (xm.lt.(2.0e-7)) xm = 2.0e-6
   endif
else if (quart) then
!  and for the quart ...
   if (xm.lt.0) xm=eps_tracer_bound(1)
   if (ym.lt.0) ym=eps_tracer_bound(1)
else if (eighth) then
!  for an eight section
   if (xm.lt.0) xm=eps_tracer_bound(1)
   if (xm.ge.ym) xm = ym-eps_tracer_bound(1)
endif

end subroutine checkbounds_cylp

subroutine checkbounds_cart(xm,ym)
use tracers
use geometry
use control
use sepmodulecomio
implicit none
real(kind=8) :: xm,ym,epsh=0.0_8
!epsh=eps_tracer_bound(1)
!if (epsh>1e-7) then
!   if (print_node) then 
!      write(irefwr,*) 'PERROR(checkbounds_cart): set eps_tracer_bound: ',eps_tracer_bound(1)
!   endif
!   call instop
!endif

xm = max(xm,xcmin+epsh)
xm = min(xm,xcmax-epsh)
ym = max(ym,ycmin+epsh)
ym = min(ym,ycmax-epsh)
 
end subroutine checkbounds_cart


subroutine stay_on_boundary()
use mpetrac
use tracers
use coeff
use geometry
implicit none
integer :: ip,ip1,ip2,nmark

if (itracoption == 2) then
!  make sure markers at boundaries stay on boundary
   if (quart) then
      ip=0
      do ichain=1,nochain
         nmark = imark(ichain)
!        x coordinate of first tracer should be zero
         ip1 = 2*ip+1
         tracer(ip)%xnew = 0
!        y coordinate of last tracer should be zero
         ip2 = 2*(ip+nmark)
         tracer(ip+nmark)%ynew = 0
         ip = ip+nmark
      enddo
   else if (half) then
      ip=0
      do ichain=1,nochain
         nmark = imark(ichain)
!        x coordinate of first tracer should be zero
         ip1 = 2*ip+1
         tracer(ip)%xnew = 0
!        x coordinate of last tracer should be zero
         ip2 = 2*(ip+nmark)-1
         tracer(ip+nmark)%xnew = 0
         ip = ip+nmark
      enddo
   else if (eighth) then
      ip=0
      do ichain=1,nochain
         nmark = imark(ichain)
!        x coordinate of first tracer should be zero
         ip1 = 2*ip+1
         tracer(ip)%xnew = 0
!        x and y coordinates of last tracer should be equal
         ip2 = 2*(ip+nmark)-1
         tracer(ip+nmark)%xnew = tracer(ip+nmark)%ynew
         ip = ip+nmark
      enddo
   endif
endif
    
end subroutine stay_on_boundary
