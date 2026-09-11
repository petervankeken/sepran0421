! Ramp function for smoothing velocities across boundaries
! The actual ramp function
real(kind=8) function ramp(x,iplate)
use geometry
use brandenburg
use coeff
implicit none

real(kind=8) :: lb,ub,fac,xr
real(kind=8) :: x,tx1,tx2,rm
real(kind=8) :: xuse
integer :: i,j,iplate,ipuse

! slope of ramp function:
rm = 1.0_8/whalf
ramp = 1.0_8

ipuse = iplate
xuse = x


lb = plate_boundaries(ipuse) - (whalf/2.0_8)
ub = plate_boundaries(ipuse+1) + (whalf/2.0_8)

tx1 = DABS(lb - xuse)
tx2 = DABS(ub - xuse)

if (tx1.lt.whalf) then
   ramp = rm*tx1
endif

if (tx2.lt.whalf) then
   ramp = rm*tx2
endif

if (xuse.lt.lb) then
   ramp = 0.0_8
endif

if (xuse.gt.ub) then
   ramp = 0.0_8
endif

if (cyl) then
   if (full) then
      ! catch the extra parts for plate 1 and plate n
      ! that are less than 0 or more than 2pi in the
      ! full annulus
      if (ipuse.eq.1) then
         if (xuse.ge.(2.0_8*pi - (whalf/2.0_8)) ) then
            ramp = rm*(xuse - (2.0_8*pi - (whalf/2.0_8)) )
         endif
      endif

      if (ipuse.eq.nplates) then
         if (xuse.le.(whalf/2.0_8)) then
            ramp = -rm*xuse + 0.5_8
         endif
      endif
   endif ! full
endif ! cyl

end function ramp
