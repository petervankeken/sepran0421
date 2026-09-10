integer function ivalfind(x,y)
use geometry
implicit none
real(kind=8) :: x,y,r
integer :: ilay

r=y
if (cyl) r=sqrt(x*x+y*y)

ilay=1
100   continue
  if (r.gt.r_zint(ilay)) then
     ivalfind = ilay
     return
  else
     ilay=ilay+1
     if (ilay.lt.nlay) then
       goto 100
     endif
     ivalfind=ilay
endif

return
end
