!   INITTRAC
!   Create initial tracer distribution
!   Combines cylindrical and cartesian versions
!
!   itracoption=1
!     Adapted to give uniform distribution in cylindrical (r,theta) space
!     Tracer coordinates are in (r,theta)
!   itracoption=2
!     Provide markerchain from theta0-theta1 radius r1
!   itracoption=9
!     Test with JGR97 benchmark
!   itracoption=10
!     Stokes spheres benchmark 2019
!   For itracoption=3 the single tracer coordinate is set in compr_start
!
!   PvK 071505
subroutine inittrac()
use sepmodulecomio
use mpetrac
use coeff
use geometry
use tracers
use control
implicit none
integer :: alloc_status

if ((compress.and.itracoption==1).and..not.stretch_tracers) then
   if (.not.ratio_method) then
      ! ratio method cares far less than Stokeslet of field method
      if (print_node) write(irefwr,*) 'PERROR(inittrac): really should stretch them tracers for compressible convection!'
      call instop
   endif
endif
if (print_node) write(irefwr,*) 'pre_allocate_tracers: ',pre_allocate_tracers
if (pre_allocate_tracers) then
   call allocate_tracers(NTRACMAX,1,'inittrac') ! default number of tracers (needs fixing)
endif
if (print_node) write(irefwr,*) 'tracers_allocated: ',tracers_allocated

if (cyl) then
   call inittrac_cyl()
else
   call inittrac_cart()
endif
return
end subroutine inittrac

