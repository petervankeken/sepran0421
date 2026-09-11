!   PEDTCF
!
!   Determines CFL timestep; for use with elements 100 or 800
!   for the heat equation.
!   This version determines size of element using all three sides 
!   of the element, which is necessary for irregular grids (like
!   that used in axisymmetric and cylindrical modeling.
!
!   PvK 040390/981106
subroutine pedtcf(kmesh,iuser,user,dtcfl)
use sepmodulecomio
use sepmodulekmesh
use coeff
use control
implicit none
integer :: kmesh,iuser(*)
real(kind=8) :: dtcfl,user(*)
integer :: ipvelx,ipvely

dtcfl = 1e30
ipvelx = iuser(8)
ipvely = iuser(9)
if (ipvelx == 0 .or. ipvely == 0) then
   if (print_node) then
      write(irefwr,*) 'PERROR(pedtcf): ipvelx or ipvely is 0: ',ipvelx,ipvely
   endif
   call instop
endif
if (npelm == 3) then
   ! linear elements for heat equation
   call pedtc1(nelem,kmeshc,coor,user(ipvelx),user(ipvely),dtcfl,npelm)
else if (npelm == 6) then
   ! quadratic elements
   call pedtc2(nelem,kmeshc,coor,user(ipvelx),user(ipvely),dtcfl,npelm)
else
   if (print_node) write(irefwr,*) 'PERROR(pedtcf): npelm <> 3,6: ',npelm
   call instop
endif
 
end subroutine pedtcf

subroutine pedtc1(nelem,kmeshc,coor,xvel,yvel,dtcfl,npelm)
implicit none
integer :: kmeshc(*),nelem,npelm
real(kind=8) :: coor(2,*),xvel(*),yvel(*),dtcfl
integer :: index(4),i,ipc,j
real(kind=8) :: u(4),v(4),dx1,dy1,dx2,dy2,dx,dy,x,y,udx,vdy
real(kind=8) :: ubiggest,vbiggest,xsmallest,ysmallest,dx3,dy3
integer :: ismallest

do i=1,nelem
   ipc = 3*(i-1)+1
   do j=1,3
      index(j) = kmeshc(ipc+j-1)
   enddo
   dx1 = abs(coor(1,index(1)) - coor(1,index(2)))
   dy1 = abs(coor(2,index(2)) - coor(2,index(3)))
   dx2 = abs(coor(1,index(2)) - coor(1,index(3)))
   dy2 = abs(coor(2,index(1)) - coor(2,index(2)))
   dx3 = abs(coor(1,index(3)) - coor(1,index(1)))
   dy3 = abs(coor(2,index(3)) - coor(2,index(1)))
   dx = max(dx1,dx2,dx3)
   dy = max(dy1,dy2,dy3)
   do j=1,3
      u(j) = xvel(index(j))
      v(j) = yvel(index(j))
      x = coor(1,index(j))
      y = coor(2,index(j))
      if (abs(u(j)).gt.1e-7.and.abs(v(j)).gt.1e-7) then
          udx = dx/abs(u(j))
          vdy = dy/abs(v(j))
          dtcfl = min(dtcfl,udx,vdy)
          ismallest = i
          xsmallest = dx
          ysmallest = dy
          ubiggest = u(j)
          vbiggest = v(j)
      endif
   enddo
enddo

end subroutine pedtc1

subroutine pedtc2(nelem,kmeshc,coor,xvel,yvel,dtcfl,npelm)
implicit none
integer :: kmeshc(*),nelem,npelm
real(kind=8) :: coor(2,*),xvel(*),yvel(*),dtcfl
integer :: index(7),i,ipc,j
real(kind=8) :: u(7),v(7),dx1,dy1,dx2,dy2,dx,dy,x,y,udx,vdy
real(kind=8) :: ubiggest,vbiggest,xsmallest,ysmallest,dx3,dy3
integer :: ismallest

do i=1,nelem
   ! find the nodal point numbers in this element from kmesh part c
   ipc = npelm*(i-1)+1
   do j=1,npelm
      index(j) = kmeshc(ipc+j-1)
   enddo
   ! find the length of the sides of the quadratic element
   dx1 = abs(coor(1,index(1)) - coor(1,index(3)))
   dy1 = abs(coor(2,index(3)) - coor(2,index(5)))
   dx2 = abs(coor(1,index(3)) - coor(1,index(5)))
   dy2 = abs(coor(2,index(1)) - coor(2,index(3)))
   dx3 = abs(coor(1,index(5)) - coor(1,index(1)))
   dy3 = abs(coor(2,index(5)) - coor(2,index(1)))
   ! find maximum size of element in x and y
   ! then divide by two for quadratic elements
   dx = max(dx1,dx2,dx3)/2
   dy = max(dy1,dy2,dy3)/2
   do j=1,npelm
      ! get the velocity in the nodal points
      u(j) = xvel(index(j))
      v(j) = yvel(index(j))
      x = coor(1,index(j))
      y = coor(2,index(j))
      ! exclude regions with very velocity
      if (abs(u(j)).gt.1e-7.and.abs(v(j)).gt.1e-7) then
         udx = dx/abs(u(j))
         vdy = dy/abs(v(j))
         dtcfl = min(dtcfl,udx,vdy)
         ismallest = i
         xsmallest = dx
         ysmallest = dy
         ubiggest = u(j)
         vbiggest = v(j)
      endif
   enddo
enddo

end subroutine pedtc2
