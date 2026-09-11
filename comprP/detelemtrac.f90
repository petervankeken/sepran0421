! DETELEMTRAC
! Determine the element number for each tracer
! 
! PvK 021304
subroutine detelemtrac(ic,calling_routine)
use msper01
use coeff
use geometry
implicit none
integer :: ic
character(len=*) :: calling_routine

if (cyl) then
   call detelemtrac_cyl(ic,calling_routine)
else
   call detelemtrac_cart(ic,calling_routine)
endif
end

