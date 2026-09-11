module msper01
contains
   subroutine sper01(kmeshc,coor,nodno,xn,yn,iel)
   implicit none
   integer :: kmeshc(*),nodno(*),iel
   real(kind=8) :: coor(2,*),xn(*),yn(*)
   integer :: inpelm,ip,i

   inpelm = 6
   ip = (iel-1)*inpelm
   do i=1,inpelm
      nodno(i) = kmeshc(ip+i)
   enddo
   do i=1,inpelm
      xn(i) = coor(1,nodno(i))
      yn(i) = coor(2,nodno(i))
   enddo

   end subroutine sper01

   subroutine sper06(user,un,nodno)
   implicit none
   real(kind=8) :: user(:),un(*)
   integer :: nodno(*),i

   do i=1,6
      un(i) = user(nodno(i))
   enddo

   end subroutine sper06

! DETSHAPE
! 
! Determine value of quadratic shapefunction of quadratic
! triangular element (xn,yn) at point (xm,ym)
! 
! PvK 950508
subroutine detshape(xn,yn,xm,ym,shapef)
implicit none
real(kind=8) :: xn(6),yn(6),xm,ym,shapef(6)
real(kind=8) :: a1,a2,a3,b1,b2,b3,c1,c2,c3,delta,rl1,rl2,rl3

a1 = xn(3)*yn(5) - yn(3)*xn(5)
a2 = -xn(1)*yn(5) + yn(1)*xn(5)
a3 = xn(1)*yn(3) - yn(1)*xn(3)
b1 = yn(3) - yn(5)
b2 = yn(5) - yn(1)
b3 = yn(1) - yn(3)
c1 = xn(5) - xn(3)
c2 = xn(1) - xn(5)
c3 = xn(3) - xn(1)
delta = -c3*b1 + b3*c1
a1 = a1/delta
a2 = a2/delta
a3 = a3/delta
b1 = b1/delta
b2 = b2/delta
b3 = b3/delta
c1 = c1/delta
c2 = c2/delta
c3 = c3/delta

rl1 = a1 + b1*xm + c1*ym     
rl2 = a2 + b2*xm + c2*ym     
rl3 = a3 + b3*xm + c3*ym     
   
shapef(1) = rl1*(2*rl1-1)
shapef(3) = rl2*(2*rl2-1)
shapef(5) = rl3*(2*rl3-1)
shapef(2) = 4*rl1*rl2
shapef(4) = 4*rl2*rl3
shapef(6) = 4*rl3*rl1
 
end subroutine detshape

! DETSHAPE7
! 
! Determine value of quadratic shapefunction of 
! extended quadratictriangular element (xn,yn) at point (xm,ym)
! This assumes that the sides of the triangle (1,3,5) are straight!
! 
! PvK 990408
subroutine detshape7(xn,yn,xm,ym,shapef)
implicit none
real(kind=8) :: xn(7),yn(7),xm,ym,shapef(7)
real(kind=8) :: a1,a2,a3,b1,b2,b3,c1,c2,c3,delta,rl1,rl2,rl3
real(kind=8) :: rl123

a1 = xn(3)*yn(5) - yn(3)*xn(5)
a2 = -xn(1)*yn(5) + yn(1)*xn(5)
a3 = xn(1)*yn(3) - yn(1)*xn(3)
b1 = yn(3) - yn(5)
b2 = yn(5) - yn(1)
b3 = yn(1) - yn(3)
c1 = xn(5) - xn(3)
c2 = xn(1) - xn(5)
c3 = xn(3) - xn(1)
delta = -c3*b1 + b3*c1
a1 = a1/delta
a2 = a2/delta
a3 = a3/delta
b1 = b1/delta
b2 = b2/delta
b3 = b3/delta
c1 = c1/delta
c2 = c2/delta
c3 = c3/delta

rl1 = a1 + b1*xm + c1*ym     
rl2 = a2 + b2*xm + c2*ym     
rl3 = a3 + b3*xm + c3*ym     
rl123 = rl1 * rl2 * rl3
   
shapef(1) = rl1*(2*rl1-1) + 3 * rl123
shapef(3) = rl2*(2*rl2-1) + 3 * rl123
shapef(5) = rl3*(2*rl3-1) + 3 * rl123
shapef(2) = 4*rl1*rl2 - 12 * rl123
shapef(4) = 4*rl2*rl3 - 12 * rl123
shapef(6) = 4*rl3*rl1 - 12 * rl123
shapef(7) = 27*rl123
 
end subroutine detshape7

end module msper01
