subroutine makecrust(ichoice)
use coeff
use tracers
use mparallel
use geometry
use brandenburg
implicit none
integer :: ichoice

if (cyl) then
   call makecrust_cyl(ichoice)
else
   call makecrust_cart(ichoice)
endif

end subroutine makecrust
