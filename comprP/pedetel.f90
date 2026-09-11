! PEDETEL
! 
! Determine in which element a point with given coordinates lies.
! ICHOICE
! 1     Mesh is rectangular, coordinates of the primary
!       nodal points are given by xc(nx) and yc(ny) 
!       The element are triangles where the diagonal of
!       the rectangle has an angle with the positive 
!       x-axis larger then 90 degrees.
!  
! 2     Mesh is rectangular, coordinates of the primary
!       nodal points are given by xc(nx) and yc(ny)
!       The elements are rectangles.
! PvK 10-8-89
subroutine pedetel(ichoice,xm,ym,ielh)
use geometry
implicit none
integer ichoice,ielh
real*8 xm,ym,ydia
integer ix,iy,iel2
real*8 x1,x2,y1,y2
  
! Find (ix,iy) such that  xc(ix-1) <= xm < xc(ix)
! and  yc(iy-1) <= ym < yc(iy)
ix = 1
100   continue
   ix = ix + 1
   if (xc(ix-1).le.xm.and.xm.le.xc(ix)) then
      continue
   else
      if (ix.ge.nx) then
         continue
      else
         goto 100
      endif
   endif
 
iy = 1
200   continue
   iy = iy + 1
   if (yc(iy-1).le.ym.and.ym.le.yc(iy)) then
      continue
   else
      if (iy.ge.ny) then
         continue
      else
         goto 200
      endif
   endif

if (ichoice.eq.2) then
   ielh = (nx-1)*(iy-2) + (ix-1)
   return
endif

! The element number is given by iel2 or iel2+1
! ix and iy are >= 2. The elements in the first row are
! numbered from 1 .. 2*(nx-1), in the second row from
! 2*(nx-1)+1 .. 4*(nx-1) etc.
iel2 = 2*(nx-1)*(iy-2) + 2*(ix-1) - 1

! Determine position of marker with respect to diagonal.
! Here it is assumed that the diagonal has a negative
! angle with the positive x-axis.
x1 = xc(ix-1)
x2 = xc(ix)
y1 = yc(iy-1)
y2 = yc(iy)
ydia = y2 + (xm-x1)/(x2-x1)*(y1-y2)
if (ym.le.ydia) then
   ielh = iel2 + 1
else
   ielh = iel2 
endif

end subroutine pedetel
