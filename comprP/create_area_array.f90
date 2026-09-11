! CREATE_AREA_ARRAY
! 
! Create an array that contains the combined area of
! the elements that a particular nodal point is in.
!
! PvK 2000s
subroutine create_area_array()
use sepmodulecomio
use sepmodulekmesh
use sepmodulevecs
use sepran_arrays
use control
implicit none
 
! kmesh(18) contains memory management sequence number for kmesh
! part D. Part D contains the elements corresponding to the nodal
! points
! check to see that the mesh is composed of linear tetrahedra
if (kelgrp.ne.3 .and. kelgrp.ne.6) then
   if (print_node) then
      write(irefwr,*) 'PERROR(create_area_array): the mesh seems to be'
      write(irefwr,*) 'composed of odd elements as kmesh(14) /= 3,6: '
      write(irefwr,*) 'kelgrp = ',kelgrp
   endif
   call instop
endif
      
if (kelgrp==3) then
   call create_area_lin(ks(iarea)%sol,kmeshc,kmshd1,kmshd2,coor,npoint)
else if (kelgrp==6) then
!  call create_area_quad(ks(iarea)%sol,kmeshc,kmshd1,kmshd2,coor,npoint)
   !write(irefwr,*) 'iarea: ',iarea,kmeshc,kmshd1,kmshd2,npoint,nelem
   call create_area_for_lin_on_quad(ks(iarea)%sol,kmeshc,kmshd1,kmshd2,coor,npoint,nelem)
endif


end subroutine create_area_array

subroutine create_area_lin(area,kmeshc,kmeshd1,kmeshd2,coor,npoint)
use sepmodulecomio
implicit none
real(kind=8) :: area(*),coor(2,*)
integer :: kmeshc(*),kmeshd1(*),kmeshd2(*),npoint
integer :: i,j,k,nelem_here,nodno,ielem
real(kind=8) :: xn(3),yn(3),a(3),b(3),c(3),area_element

!if (print_node) write(irefwr,*) 'npoint'
do i=1,npoint
   area(i) = 0d0
   nelem_here = kmeshd1(i+1)-kmeshd1(i)
   do j=1,nelem_here
      ielem = kmeshd2(kmeshd1(i)+j)
!     write(irefwr,'(i10)') ielem
      do k=1,3
         nodno = kmeshc(3*ielem-3+k)
         xn(k) = coor(1,nodno)
         yn(k) = coor(2,nodno)
!        write(irefwr,'(''    '',2f12.3)') xn(k),yn(k)
      enddo
      call trilin(xn,yn,a,b,c,area_element)
      area(i)=area(i)+area_element
   enddo
!  write(irefwr,'(''area: '',i10,2e15.7,i5)') i,area(i),area(i)/nelem_here,nelem_here
enddo
!call instop

end subroutine create_area_lin

subroutine create_area_quad(area,kmeshc,kmeshd1,kmeshd2,coor,npoint)
use sepmodulecomio
use control
use msper01
implicit none
real(kind=8) :: area(*),coor(2,*)
integer :: kmeshc(*),kmeshd1(*),kmeshd2(*),npoint,nodlin(3)
integer :: i,j,k,nelem_here,nodno(6),ielem,ip,findsubelem,isub
real(kind=8) :: xn(6),yn(6),a(3),b(3),c(3),area_element
real(kind=8) :: xl(3),yl(3),xm,ym,rl(3)
integer,dimension(12) :: lindex

lindex=(/1,2,6, 2,3,4, 4,5,6, 2,6,4/)

write(irefwr,*) 'npoint'
do i=1,npoint
   area(i) = 0d0
   nelem_here = kmeshd1(i+1)-kmeshd1(i)
   xm=coor(1,i)
   ym=coor(2,i)
   do j=1,nelem_here
      ielem = kmeshd2(kmeshd1(i)+j)
      ! first grab coordinates of this quadratic element
      call sper01(kmeshc,coor,nodno,xn,yn,ielem)

      ! if we really want to use a P2 triangle we should use something like
      ! this:
      ! now subdivide into 4 linear triangles
      ! this should be replaced by numerical integration for
      ! general isoparametric P2 triangles as the below is
      ! only accurate for straight-sided P2 triangles
      ip=0
      do k=1,4
         xl(1)=xn(lindex(ip+1))
         yl(1)=yn(lindex(ip+1))
         xl(2)=xn(lindex(ip+2))
         yl(2)=yn(lindex(ip+2))
         xl(3)=xn(lindex(ip+3))
         yl(3)=yn(lindex(ip+3))
         call trilin(xl,yl,a,b,c,area_element)
         area(i)=area(i)+area_element
         ip=ip+3
       enddo ! k

      ! but until that time find neighboring subtriangle
      ! isub = findsubelem(xn,yn,xm,ym,xl,yl,rl,nodno,nodlin)
      ! and use that area
      ! call trilin(xl,yl,a,b,c,area_element)
      ! area(i)=area(i)+area_element
   enddo ! nelem_here
   if (i<10.and.print_node) write(irefwr,*) 'area: ',area(i)
enddo
call instop

end subroutine create_area_quad

! special case when we want to compute the area of the surrounding linear (subdivided) triangle
! when we have a quadratic mesh. 
! Loop over elements: 
!    find area of subtriangle
!    add this area to the points that belong to this subtriangle
!call create_area_for_lin_on_quad(ks(iarea)%sol,kmeshc,kmshd1,kmshd2,coor,npoint,nelem)
subroutine create_area_for_lin_on_quad(area,kmeshc,kmeshd1,kmeshd2,coor,npoint,nelem)
use sepmodulecomio
use msper01
implicit none
real(kind=8) :: area(*),coor(2,*)
integer :: kmeshc(*),kmeshd1(*),kmeshd2(*),npoint,nodlin(3),nelem
integer :: i,j,k,nelem_here,nodno(6),ielem,ip,findsubelem,isub
real(kind=8) :: xn(6),yn(6),a(3),b(3),c(3),area_element
real(kind=8) :: xl(3),yl(3),xm,ym,rl(3)
integer,dimension(12) :: lindex

lindex=(/1,2,6, 2,3,4, 4,5,6, 2,6,4/)

do ielem=1,nelem
   !write(irefwr,*) 'ielem: ',ielem
   area(ielem) = 0d0
   ! first grab coordinates of this quadratic element
   call sper01(kmeshc,coor,nodno,xn,yn,ielem)

   ! if we really want to use a P2 triangle we should use something like
   ! this:
   ! now subdivide into 4 linear triangles
   ! this should be replaced by numerical integration for
   ! general isoparametric P2 triangles as the below is
   ! only accurate for straight-sided P2 triangles
   ip=0
   do k=1,4
      xl(1)=xn(lindex(ip+1))
      yl(1)=yn(lindex(ip+1))
      xl(2)=xn(lindex(ip+2))
      yl(2)=yn(lindex(ip+2))
      xl(3)=xn(lindex(ip+3))
      yl(3)=yn(lindex(ip+3))
      call trilin(xl,yl,a,b,c,area_element)
      area(nodno(lindex(ip+1)))=area(nodno(lindex(ip+1)))+area_element
      area(nodno(lindex(ip+2)))=area(nodno(lindex(ip+2)))+area_element
      area(nodno(lindex(ip+3)))=area(nodno(lindex(ip+3)))+area_element
       ip=ip+3
   enddo ! k
!  do i=1,6
!     if (nodno(i)<10) then
!         write(irefwr,'(''point subelement: '',3i5,2e15.3)') i,k,nodno(i),area(nodno(i))
!     endif
!  enddo

enddo

end subroutine create_area_for_lin_on_quad
