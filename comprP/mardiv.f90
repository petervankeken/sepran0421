! MARDIV
! Create a pixel representation of the layer below a markerchain
! This subroutine combines the cylindrical and Cartesian box versions
! PvK 071504
subroutine mardiv(rname,iout)
use sepmodulecomio
use control
use coeff
use geometry
use tracers
implicit none
character(len=*) :: rname
integer :: iout
     
if (cyl) then
   !call mardiv_cyl(coormark,rname,iout)
   if (print_node) write(irefwr,*) 'PERROR(mardiv): cyl version not yet adjusted for coormark->tracer'
   call instop
else
   call mardiv_cart(rname,iout)
endif
     
end subroutine mardiv

! MARDIV_CART (formerly known as MARRAS)
! 
! Create a rasterfile in which the value of the pixels
! indicate the layer in which it is positioned.
! The layerboundaries are represented by the boundaries
! of the geometry ([0,rlampix]x[0,1]) and the markerchains
! The algorithm is described in PVK 040690
! 
! PvK 050690/300790
! PvK 921013 : correction for y-direction: midpoint
! PvK 072704: change pixel grid coordinates to 'normal' x,y orientation
!                 with iy=1 corresponding to y=0.
! PvK 071705: reverse pixel values: 1=bottom layer; 0=top layer.
!                 This is in general more efficient since the bottom layer
!                 is typically thinner. 
subroutine mardiv_cart(rname,iout)
use sepmodulecomio
use tracers
use coeff
use mpetrac
use control
implicit none
character(len=*) :: rname
integer, parameter :: NYMAX=10000
real(kind=8) :: cy(NYMAX)
integer :: icy(NYMAX),iout
real(kind=8) :: dyl,xm,ym,x1,x2,y1,y2,temp
integer :: i,ip,nmark,ix,iy,j,it,iregio,icut,im,ic,ntot,iadd
integer :: icol(256)
data icol(1),icol(2),icol(3),icol(4),icol(5)/128,1,180,200,255/

if (rlampix.le.1e-7) then
   if (print_node) write(irefwr,*) 'PERROR(mardiv_cart): rlam=0: ',rlampix
   call instop
endif
if (nxpix.eq.0) then
   if (print_node) write(irefwr,*) 'PERROR(mardiv_cart): nxpix=0'
   call instop
endif
if (nypix>NYMAX) then
    if (print_node) write(irefwr,*) 'PERROR(mardiv_cart: nypix>NYMAX: ',nypix,NYMAX
    call instop
endif
dyl = 1d0/nypix
! Initialize
do i=1,nxpix*nypix
   ipix(i) = 0
enddo
ip = 0
do ichain=1,nochain
   nmark = imark(ichain)
!  Loop over columns
   do ix=1,nxpix
!     Determine in which pixels the markerchain
!     intersects the column. Start with the lower
!     boundary (y=0,iy=1) and finish with the
!     upper one (y=1,iy=nypix) 
      xm = (0.5+(ix-1))*rlampix/nxpix
      icut=1
      cy(1)=0d0
      do im=2,nmark
         if (ibuoy_trac==2) then
           x1 = tracer(ip+im-1)%xnew
           x2 = tracer(ip+im)%xnew
         else
           x1 = tracer(ip+im-1)%x
           x2 = tracer(ip+im)%x
         endif
         if (x1.le.xm.and.x2.gt.xm.or.x1.ge.xm.and.x2.lt.xm) then
!           the markerchain is intersecting the current column
            if (ibuoy_trac==2) then
               y1 = tracer(ip+im-1)%ynew
               y2 = tracer(ip+im)%ynew
            else
               y1 = tracer(ip+im-1)%y
               y2 = tracer(ip+im)%y
            endif
            ym = y1+(xm-x1)/(x2-x1)*(y2-y1)
            icut = icut+1
            cy(icut) = ym
         endif
      enddo
      icut=icut+1
      cy(icut)=1d0
!     sort this array
      do i=2,icut-1
         do j=3,icut-1
            if (cy(j-1).gt.cy(j)) then
               temp = cy(j-1)
               cy(j-1) = cy(j)
               cy(j) = temp
            endif
         enddo
      enddo
!     Figure out which pixels are below the markerchain
      iadd=0
!     go from the top down; top layer will be indicated by '0', bottom layer '1'
!     Change the above line to 'iadd=1' for the reverse
      do iy=nypix,1,-1
         ym = (iy-0.5)*dyl
111      continue
            if (ym.le.cy(icut).and.ym.ge.cy(icut-1)) then
!               pixel lies in current interval
                icy(iy) = iadd
            else
!               Try next interval
                icut = icut-1
                iadd = mod(iadd+1,2)
                if (icut.ge.2) then
                   goto 111
                else 
                   if (print_node) write(irefwr,*) 'problem in mardiv'
                   call instop
                endif
            endif
 
       enddo
!      Now add this to IPIX
       do iy=1,nypix
          ic = (iy-1)*nxpix + ix
          ipix(ic)=ipix(ic)+icy(iy)
       enddo
    enddo
    ip = ip+2*nmark
enddo

if (iout.eq.-1) then
   ntot = nxpix*nypix
!  Represents the layers 1,2,3 ... with colors stored in icol
   do ip=1,ntot
      ipix(ip) = icol(ipix(ip)+1)
   enddo
   call rasout(ntot,ipix,rname)
else if (iout==-2.and.print_node) then
   open(99,file=rname)
   do ix=1,nxpix
      do iy=1,nypix
         ic=(iy-1)*nxpix+ix
         if (ipix(ic)==1) then
            write(99,*) (ix-0.5)*dxpix,(iy-0.5)*dypix
         endif
      enddo
   enddo
endif

end subroutine mardiv_cart


 
