!     ICHECKINELEM
!  
!     Check of point (xm,ym) is in the quadratic triangular
!     element (xn,yn).  
!     ichoice   i  1  the barycentric coordinates of the linear triangle 
!                      made by the vertices should be positive
!                  2  test the barycentric coordinates of the 4 
!                     linear subtriangles separately
!                  3  assume the element is curved; use Newton iteration
!                     to obtain barycentric coordinates 
!     xn,yn     i  Coordinates of the nodal points
!     xm,ym     i  Coordinates of the tracer
!     nodno     i  Nodal point numbers
!     nodlin    i  nodal point numbers of linear sub element (if ichoice=2)
!     rl        o  barycentric coordinates 
!     xi,eta    o  barycentric coordinates
!     phiq      o  quadratic shapefunctions 
!     isub      o  number of linear sub element (if ichoice=2)
!  
!     PvK 950508
!  
!     PvK 970413: modified for cylindrical coordinates. Make sure
!                 to test each of the 4 subtriangles individually
!     PvK 021004: modified to allow for curved elements
logical function checkinelem(ichoice,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,ielh)
use sepmodulecomio
use control
implicit none
real(kind=8) :: xn(6),yn(6),xm,ym,rl(3),xi,eta,phiq(6)
integer :: ichoice,isub,nodno(6),nodlin(3),i,ielh,ludcmp_error
real(kind=8),parameter :: epsh=1.0e-4_8
logical :: out,fail,okay


checkinelem=.false.


if (ichoice.eq.1) then

   ! check just the barycentric coordinates
   isub=0
   call findrl(xn,yn,xm,ym,rl,isub,nodno,nodlin)
   out = (rl(1).ge.-epsh).and.(rl(2).ge.-epsh).and.(rl(3).ge.-epsh)  
   if (out) then
      checkinelem=.true.
   else
      checkinelem=.false.
   endif
   if (verbose) write(irefwr,*) 'rl: ',rl(1),rl(2),rl(3)
   xi = rl(1)
   eta = rl(2)
   return

else if (ichoice.eq.2) then          

   ! check each four of the sub triangles separately
   isub=1
100  continue
     call findrl(xn,yn,xm,ym,rl,isub,nodno,nodlin)
     out = (rl(1).ge.-epsh).and.(rl(2).ge.-epsh).and.(rl(3).ge.-epsh)
     if (out) then
         checkinelem=.true.
         xi = rl(1)
         eta = rl(2)
         return
      else if (isub.lt.4) then
         isub=isub+1
         goto 100
      else
         checkinelem=.false.
         return
      endif

else if (ichoice.eq.3) then

! find local coordinates in reference triangle
! with corners (0,0) (1,0) (0,1)
      ! First check to make sure the point is not very far from the element
      ! This would cause find_xi_eta to fail
      isub=0
      call findrl(xn,yn,xm,ym,rl,isub,nodno,nodlin)
      okay = (rl(1).ge.-0.5).and.(rl(2).ge.-0.5).and.(rl(3).ge.-0.5)
      okay = okay.and.(rl(1).le.1.5).and.(rl(2).le.1.5).and.(rl(3).le.1.5)
      fail = .not.okay
      if (fail) then 
         checkinelem = .false.
      else 
         if (verbose) then
            write(irefwr,*) 'before find_xi_eta:  x= ',xn(1),xn(3),xn(5),xm
            write(irefwr,*) 'before find_xi_eta:  y= ',yn(1),yn(3),yn(5),ym
            write(irefwr,*) '                    rl= ',rl(1),rl(2),rl(3)
         endif
         call find_xi_eta(xn,yn,xm,ym,xi,eta,phiq,fail,ludcmp_error)
         if (ludcmp_error/=0) then
            if (print_node) write(irefwr,*) 'find_xi_eta failed on ludcmp_error=',ludcmp_error
            call instop
         endif
         if (xi.ge.-epsh.and.eta.ge.-epsh.and.(xi+eta).le.1+epsh) checkinelem = .true.
         if (verbose) then
             write(irefwr,'(''find_xi_eta: '',4e15.7,2L5)') xm,ym,xi,eta,okay,checkinelem
             write(irefwr,'(''             '',6e15.7)') xn(1),yn(1),xn(3),yn(3),xn(5),yn(5)
             write(irefwr,'(''         rl: '',3e15.7)') rl(1:3)
             !open(9,file='elem.dat') 
             !do i=1,6
             !   write(9,*) xn(i),yn(i)
             !enddo
             !open(9,file='coor.dat')
             !write(9,*) xm,ym
             !close(9)
             write(irefwr,'('' xi, eta, xi+eta: '',3e15.7)') xi,eta,xi+eta
         endif
      endif
      return
else 
 
    if (print_node) then
       write(irefwr,*) 'PERROR(checkinelem): unkwown option for ichoice ' 
       write(irefwr,*) 'ichoice = ',ichoice
    endif
    call instop
endif

end function checkinelem

subroutine findrl(xn,yn,xm,ym,rl,isub,nodno,nodlin)
implicit none
real(kind=8) :: xn(6),yn(6),xm,ym,rl(3)
real(kind=8) :: a(3),b(3),c(3),xl(3),yl(3),area
integer :: isub,i,nodno(*),nodlin(*)
integer :: inod(4,3)
data (inod(1,i),i=1,3)/1,2,6/
data (inod(2,i),i=1,3)/2,3,4/
data (inod(3,i),i=1,3)/2,4,6/
data (inod(4,i),i=1,3)/4,5,6/

if (isub.eq.0) then
!        *** Use vertices of quadratic element
   xl(1) = xn(1)
   yl(1) = yn(1)
   xl(2) = xn(3)
   yl(2) = yn(3)
   xl(3) = xn(5)
   yl(3) = yn(5)
else
!        *** use the ISUBth linear subtriangle
   xl(1) = xn(inod(isub,1))
   yl(1) = yn(inod(isub,1))
   xl(2) = xn(inod(isub,2))
   yl(2) = yn(inod(isub,2))
   xl(3) = xn(inod(isub,3))
   yl(3) = yn(inod(isub,3))
endif

call trilin(xl,yl,a,b,c,area)

rl(1) = a(1) + b(1)*xm + c(1)*ym
rl(2) = a(2) + b(2)*xm + c(2)*ym
rl(3) = a(3) + b(3)*xm + c(3)*ym

if (isub.eq.0) then
   nodlin(1) = nodno(1)
   nodlin(2) = nodno(3)
   nodlin(3) = nodno(5)
else
   nodlin(1) = nodno(inod(isub,1))
   nodlin(2) = nodno(inod(isub,2))
   nodlin(3) = nodno(inod(isub,3))
endif

end subroutine findrl
