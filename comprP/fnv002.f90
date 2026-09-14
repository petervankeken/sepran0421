! special kind called for weakzones
real(kind=8) function fnv002(x,y,z,u,v,w,secinv)
use control
use tracers
implicit none
real(kind=8) :: x,y,z,u,v,w,secinv
real(kind=8) :: temp=0.0,adia=0.0,pefvis,A0
integer :: ival=1

fnv002=1.0_8
if (weak_zones) then
   ! hack
   A0=1.0_8
   if (y>=0.95) then
      A0=eta_plate! lithosphere
      if (x<=0.05.or.x>=rlampix-0.05) A0=eta_wz
   endif
   fnv002=A0
endif

end function fnv002
