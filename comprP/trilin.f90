! *************************************************************
! *    TRILIN
! *
! *    Determine coefficients of shapefunctions of linear triangle
! *    
! *    PvK 170590
! *************************************************************
subroutine trilin(x,y,a,b,c,area_element)
use sepmodulecomio
use control
implicit none
real(kind=8) ::  x(3),y(3),a(3),b(3),c(3)
real(kind=8) :: delta,area_element,inv_delta
integer :: i
 
a(1) = x(2)*y(3) - y(2)*x(3)
a(2) = -x(1)*y(3) + y(1)*x(3)
a(3) = x(1)*y(2) - y(1)*x(2)
b(1) = y(2)-y(3)
b(2) = y(3)-y(1)
b(3) = y(1)-y(2)
c(1) = x(3)-x(2)
c(2) = x(1)-x(3)
c(3) = x(2)-x(1)
delta = -c(3)*b(1) + b(3)*c(1)
area_element = abs(delta)*0.5
if (delta.eq.0.and.print_node) then
   write(irefwr,*) 'PERROR(trilin): surface of element is zero'
   call instop
endif
inv_delta=1d0/delta
do i=1,3
   a(i) = a(i)*inv_delta
   b(i) = b(i)*inv_delta
   c(i) = c(i)*inv_delta
enddo

end subroutine trilin
