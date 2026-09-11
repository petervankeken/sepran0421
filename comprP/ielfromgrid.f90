!   IELFROMGRID
!
!   Determine the element that contains (x,y). IELTOGRID
!   should have been used to store element numbers in the
!   equidistant array ielgrid(nxel,nyel).
!
!   If the correct element is found the element & interpolation info
!   is returned in ielh, nodno, xn, yn, phiq, xi, eta.
!
!   /pieltogrid/
!      dxel  - grid spacing in x
!      dyel  - grid spacing in y
!      iel   - array containing element numbers 
!      nxel  - x-dimension of iel
!      nyel  - y-dimension of iel
!
!   PvK 970411
subroutine ielfromgrid(ichoice,x,y,kmeshc,coor,nodno,xn,yn,iel_now,ielh,phiq,xi,eta)
use sepmodulecomio
use geometry
use control
implicit none
integer :: ichoice,ielh,kmeshc(*),nodno(*),iel_now(*)
real(kind=8) :: x,y,coor(2,*),xn(*),yn(*),phiq(*)
logical :: out

integer :: i,j,inpelm,ip,ix1,ix2,iy1,iy2
integer :: ith0,ith3,ir0,ir3
integer :: isub,nodlin(3)
logical :: checkinelem
real(kind=8) :: rl(3),xi,eta

!if (print_node) verbose=.true.
 
! estimate location in regular grid based on (x,y)
ix1 = int((x+xoff)/dxel) + 1
if (ix1.eq.nxel) then
   ix1=ix1-1
endif
ix2 = ix1+1
iy1 =int((y+yoff)/dyel) + 1
if (iy1.eq.nyel) then
   iy1=iy1-1
endif
iy2 = iy1 + 1
if (verbose.and.print_node) write(irefwr,'(''ix1,etc: '',6i10)') ix1,ix2,iy1,iy2,nxel,nyel

! find the elements indicated by iel
iel_now(1) = ielgrid(ix1,iy1)
iel_now(2) = ielgrid(ix2,iy1)
iel_now(3) = ielgrid(ix2,iy2)
iel_now(4) = ielgrid(ix1,iy2)
if (verbose.and.print_node) then
   write(irefwr,'(''iel_now before duplicates: '',4i10)') (iel_now(i),i=1,4)
endif
! Check for duplicates
if (iel_now(2).eq.iel_now(1)) then
   iel_now(2)=0
endif
if (iel_now(3).eq.iel_now(1).or.iel_now(3).eq.iel_now(2)) then 
   iel_now(3)=0
endif
if (iel_now(4).eq.iel_now(1).or.iel_now(4).eq.iel_now(2).or. iel_now(4).eq.iel_now(3)) then 
   iel_now(4)=0
endif


! test coordinates to see which element 
inpelm=6
if (verbose.and.print_node) then
   write(irefwr,'(''iel_now: '',4i10)') (iel_now(i),i=1,4)
endif
i=1
100   continue
    if (iel_now(i).gt.0) then
       ip = (iel_now(i)-1)*inpelm
       do j=1,inpelm
          nodno(j) = kmeshc(ip+j)
          xn(j) = coor(1,nodno(j))
          yn(j) = coor(2,nodno(j))
       enddo
       if (verbose.and.print_node) write(irefwr,*) 'curved_elem: ',iel_now(i),curved_elem(iel_now(i))
       if (curved_elem(iel_now(i))) then
          !if (iel_now(i)==8284) write(irefwr,*) 'checkinelem(3)'
          out = .not.checkinelem(3,xn,yn,x,y,nodno,nodlin,rl,xi,eta,phiq,isub,iel_now(i))
       else
          !if (iel_now(i)==8284) write(irefwr,*) 'checkinelem(1)'
          out = .not.checkinelem(1,xn,yn,x,y,nodno,nodlin,rl,xi,eta,phiq,isub,iel_now(i))
       endif
       if (verbose.and.print_node) then
          write(irefwr,'(''iel_now: '',2i10)') i,iel_now(i)
          write(irefwr,'(''xm,ym  : '',10f12.3,L5)') x,y,xn(1),yn(1),xn(3),yn(3),xn(5),yn(5),xi,eta,out
       endif
    else
       out = .true.
    endif
    if (out) then
       if (i.lt.4) then
          i=i+1
          goto 100
       else
          ielh = -1
       endif
    else
       ielh=iel_now(i)
       return
    endif

end subroutine ielfromgrid
