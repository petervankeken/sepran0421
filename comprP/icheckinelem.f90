!   ICHECKINELEM
!
!   Check of point (xm,ym) is in the quadratic triangular
!   element (xn,yn).  For this, the barycentric coordinates
!   should all be positive
!
!   PvK 950508
integer function icheckinelem(xn,yn,xm,ym)
implicit none
real(kind=8) :: xn(6),yn(6),xm,ym
real(kind=8) :: a1,a2,a3,b1,b2,b3,c1,c2,c3,delta,rl1,rl2,rl3
logical :: out

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

out = .true.
rl1 = a1 + b1*xm + c1*ym
rl2 = a2 + b2*xm + c2*ym
rl3 = a3 + b3*xm + c3*ym
out = (rl1.ge.-1d-5).and.(rl2.ge.-1d-5).and.(rl3.ge.-1d-5)
icheckinelem=-999
if (out) icheckinelem=1

end function icheckinelem
