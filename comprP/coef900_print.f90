subroutine coef900_print(npoint,coor,rho,f1,f2,eta)
use sepmodulecomio
use geometry
use coeff
implicit none
integer :: npoint
real(kind=8) :: coor(2,npoint),rho(npoint),f1(npoint),f2(npoint)
real(kind=8) :: eta(npoint),x,y,r
integer :: i

if (cyl) then
 write(irefwr,*) 'cylindrical coordinates'
!i=npoint
!x = coor(1,i)
!y = coor(2,i)
!r = sqrt(x*x+y*y)
!write(irefwr,'(''         '',7f12.3)') x,y,r,rho(i),f1(i),f2(i),eta(i)

! do i=1,npoint
!   x = coor(1,i)
!   y = coor(2,i)
!   r = sqrt(x*x+y*y)
!   if (r>2.205) then
!      write(irefwr,'(''         '',7f12.3)') x,y,r,rho(i),f1(i),f2(i),eta(i),i*1.0/npoint
!   endif
! enddo
! call instop

else
 write(irefwr,*) 'Cartesian coordinates'
endif
write(irefwr,'(''coef900: '',7a12)') 'x','y','r','rho','f1','f2','eta'
do i=1,npoint,npoint/10
   x = coor(1,i)
   y = coor(2,i)
   r = sqrt(x*x+y*y)
   write(irefwr,'(''         '',7e12.4)') x,y,r,rho(i),f1(i),f2(i),eta(i)
enddo
   
end subroutine coef900_print
