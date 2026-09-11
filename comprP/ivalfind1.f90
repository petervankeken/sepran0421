!   Determine layer in which (x,y) sits
integer function ivalfind1(x,y)
use geometry
implicit none
real(kind=8) :: r,sint,cost,x,y
integer :: ilay

r = sqrt(x**2 + y**2)

ilay=1
100   continue
  if (r.gt.r_zint(ilay)) then
     ivalfind1=ilay
     return
  else
     ilay=ilay+1
     if (ilay.lt.nlay) then
        goto 100
     endif
     ivalfind1=ilay
endif

end function ivalfind1
