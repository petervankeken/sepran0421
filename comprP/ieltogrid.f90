! IELTOGRID
! 
! Find the mapping of element numbers to a regular Cartesian grid.
! This greatly facilitates tracer interpolation  (see ielfromgrid)
! PvK 970406
! 
! Rewritten for curved elements. This should at some point be made
! more efficient by an iterative scheme: 1) fill regular grid using
! assumption of straight elements; 2) update regular grid using
! curved elements.
! 
! PvK 020304
! 
! Modified for cylindrical lookup table (for itracoption=2).
! Still keep Cartesian one around for tracer integration.
! Put (r,th) pixel grid together for buoyancy computation.
! PvK 092204
! *************************************************************
subroutine ieltogrid()
use sepmodulecomio
use sepmodulecpack
use sepmodulemain
use sepmodulekmesh
use sepmodulemeshinf
use sepmodulemesh
use convparam
use coeff
use tracers
use geometry
use control
use msper01
implicit none
integer :: ikelmc,ikelmi,ielem

real(kind=8) :: xn(6),yn(6)
real(kind=8) :: a(3),b(3),c(3),x,y,rl(3),r
real(kind=8) :: xcminh,xcmaxh,ycminh,ycmaxh
real(kind=8) :: dx1,dy1,dx2,dy2,dx3,dy3,dr1,dr2,dr3,rmin,rmax
integer :: ipx,ipy,i,ip
integer :: nodno(6),isub
integer :: iniget,inidgt,nodlin(3)
logical :: checkinelem,correct
integer :: imissed,ifound,ix,iy
integer :: ipxmin1,ipxmax1,ipymin1,ipymax1
real(kind=8) :: xi,eta,phiq(6),dr_elem,un(6),vn(6),user,shapef(7)
logical :: out,check,flip,guess_first
logical :: on_mesh,off_mesh,holes=.false.
integer :: ielx,iely,iela(4,2),ielnow,ir,ith,ncircle
real(kind=8) :: xm,ym,rm,theta,prefac,rpix,thpix
real(kind=8) :: th1,th2,th3,th4,rm1,rm2,rm3,rm4,thmin,thmax
real(kind=8) :: xmin,xmax,ymin,ymax,dthcircle


! *** prescribe the dimensions of the regular grid
if (axi) then
   if (quart.or.eighth) then
      NXEL = NXELM
      NYEL = NYELM
      dxel = r2/(NXEL-1)
      dyel = r2/(NYEL-1)
   else
      NXEL = (NXELM-1)/2+1
      NYEL = NYELM
      dxel = r2/(NXEL-1)
      dyel = 2*r2/(NYEL-1)
   endif
   if (print_node) write(irefwr,'(''  PINFO(ieltogrid): axi: '',2i5,2f8.3)') NXEL,NYEL,dxel,dyel
else ! cyl
   if (quart.or.eighth) then
      NXEL = NXELM
      NYEL = NYELM
      dxel = r2/(NXEL-1)
      dyel = r2/(NYEL-1)
   else if (half) then
      NXEL = (NXELM-1)/2+1
      NYEL = NYELM
      dxel = r2/(NXEL-1)
      dyel = 2*r2/(NYEL-1)
   else
      NXEL = NXELM
      NYEL = NYELM
      dxel = 2*r2/(NXELM-1)
      dyel = 2*r2/(NYELM-1)
   endif
   if (print_node) write(irefwr,'(''  PINFO(ieltogrid): cyl: '',2i5,2e12.5)') NXEL,NYEL,dxel,dyel
endif

rmin = r2
rmax = 0
if (axi) then
   if (quart.or.eighth) then
      xoff = 0
      yoff = 0
      grid_height=r2
      grid_width =r2
   else
      xoff=0
      yoff=r2
      grid_height=2*r2
      grid_width =r2
   endif
else ! cyl
   if (quart.or.eighth) then
      xoff = 0
      yoff = 0
      grid_height=r2
      grid_width =r2
   else if (half) then
      xoff=0
      yoff=r2
      grid_height=2*r2
      grid_width =r2
   else
      xoff=r2
      yoff=r2
      grid_height=2*r2
      grid_width =2*r2
   endif
endif
if (print_node) then
   write(irefwr,'(''  PINFO(ieltogrid): grid width/height   : '',2f8.3)') grid_width,grid_height
   write(irefwr,'(''                    xoff yoff           : '',2f8.3)') xoff,yoff
   write(irefwr,'(''                   axi eighth quart half: '',4L8)') axi,eighth,quart,half
endif

if (periodic) then
   if (print_node) then
      write(irefwr,*) 'PERROR(ieltogrid): needs update for periodic'
   endif
   call instop
endif

! initialize element numbers in regular grid to aid error check
do ipy = 1,NYEL
   do ipx = 1,NXEL
      ielgrid(ipx,ipy) = -1
   enddo
enddo

if (nrpix.le.60) open(99,file='elements.dat')


! loop over elements
do ielem=1,nelem
!  determine coordinates for nodal points of this element
   call sper01(kmeshc,coor,nodno,xn,yn,ielem)
   if (nrpix<=60) then
      do i=1,6
         write(99,*) xn(i),yn(i)
      enddo
      write(99,*) xn(1),yn(1)
      write(99,'(''>'')') 
   endif
!  Cartesian lookup grid
   xcminh = minval(xn(1:6))
   xcmaxh = maxval(xn(1:6))
   ycminh = minval(yn(1:6))
   ycmaxh = maxval(yn(1:6))
   dx1 = xn(1)-xn(3)
   dy1 = yn(1)-yn(3)
   dx2 = xn(3)-xn(5)
   dy2 = yn(3)-yn(5)
   dx3 = xn(5)-xn(1)
   dy3 = yn(5)-yn(1)
   dr1 = sqrt(dx1*dx1+dy1*dy1)
   dr2 = sqrt(dx2*dx2+dy2*dy2)
   dr3 = sqrt(dx3*dx3+dy3*dy3)
   rmin = min(rmin,dr1,dr2,dr3)
   rmax = max(rmax,dr1,dr2,dr3)
!  translate this to grid element values
   ipxmin(ielem) = int(max(1.0,(xcminh+xoff)/dxel+1.0))
   ipxmax(ielem) = int(min(NXEL*1.0,(xcmaxh+xoff)/dxel+1.0))
   ipymin(ielem) = int(max(1d0,(ycminh+yoff)/dyel+1.5))
   ipymax(ielem) = int(min(NYEL*1.0,(ycmaxh+yoff)/dyel+1.5))

   if (curved_elem(ielem)) then
     dr_elem =  min(dr1,dr2,dr3)
     ipxmin(ielem) = ipxmin(ielem) - 0.1*dr_elem/dxel
     ipxmax(ielem) = ipxmax(ielem) + 0.1*dr_elem/dxel
     ipymin(ielem) = ipymin(ielem) - 0.1*dr_elem/dyel
     ipymax(ielem) = ipymax(ielem) + 0.1*dr_elem/dyel
     ipxmin(ielem) = max(1,ipxmin(ielem))
     ipxmax(ielem) = min(NXEL,ipxmax(ielem))
     ipymin(ielem) = max(1,ipymin(ielem))
     ipymax(ielem) = min(NYEL,ipymax(ielem))
!    The ipxmin/ipxmax etc. values are overestimates.
!    use ipxmin1 etc. to reduce the values.
     ipxmin1 = ipxmin(ielem)
     ipxmax1 = ipxmax(ielem)
     ipymin1 = ipymin(ielem)
     ipymax1 = ipymax(ielem)
     do ipy = ipymin(ielem),ipymax(ielem)
        do ipx = ipxmin(ielem),ipxmax(ielem)
           x = (ipx-1)*dxel - xoff
           y = (ipy-1)*dyel - yoff
           r = sqrt(x*x+y*y)
           if (r<r1.or.r>r2) then
              correct=.false.
           else if (curved_elem(ielem)) then
              correct = checkinelem(3,xn,yn,x,y,nodno,nodlin,rl,xi,eta,phiq,isub,ielem)
           else
              correct = checkinelem(1,xn,yn,x,y,nodno,nodlin,rl,xi,eta,phiq,isub,ielem)
           endif
           if (correct) then
!             *** this grid point is in the element
              ielgrid(ipx,ipy) = ielem
!             *** improve estimate of min/max pixel extent
              ipxmin1=min(ipx,ipxmin1)
              ipxmax1=max(ipx,ipxmax1)
              ipymin1=min(ipy,ipymin1)
              ipymax1=max(ipy,ipymax1)
           endif
         enddo
     enddo
     ipxmin(ielem) = ipxmin1
     ipxmax(ielem) = ipxmax1
     ipymin(ielem) = ipymin1
     ipymax(ielem) = ipymax1
   else

!    Element is straight
     do ipy = ipymin(ielem),ipymax(ielem)
        do ipx = ipxmin(ielem),ipxmax(ielem)
           x = (ipx-1)*dxel - xoff
           y = (ipy-1)*dyel - yoff
           r = sqrt(x*x+y*y)
           if (r<r1.or.r>r2) then
              correct=.false.
           else
              correct = checkinelem(1,xn,yn,x,y,nodno,nodlin,rl,xi,eta,phiq,isub,ielem)
           endif
           if (correct) then
!             *** this grid point is in the element
              ielgrid(ipx,ipy) = ielem
           endif
         enddo
      enddo
   endif
enddo

if (NYELM<60) then
  do ipy=NYELM,5,-5
     if (print_node) write(irefwr,'(i5,'': '',40i4,:)') ipy,(ielgrid(ipx,ipy),ipx=5,NXELM,5)
  enddo
endif

if (print_node) write(irefwr,*) 'check on holes in grid'
do ipy=1,NYELM
   do ipx=1,NXELM
      if (ielgrid(ipx,ipy)==0) then
         holes=.true.
         if (print_node) write(irefwr,*) ipx,ipy
      endif
   enddo
enddo
if (holes) then
   if (print_node) write(irefwr,*) 'look up grid has holes'
   call instop
endif

if (print_node) then
  write(irefwr,'(''PINFO(ieltogrid): '')')
  write(irefwr,'(''  Minimum size of quadratic elements  : '',f8.5,'' ('',f8.1,'' km)'')')rmin,rmin*height_dim*1e-3
  write(irefwr,'(''  Maximum size of quadratic elements  : '', f8.5,'' ('',f8.1,'' km)'')') rmax,rmax*height_dim*1e-3
  write(irefwr,'(''  Resolution of regular grid          : '', f8.5,'' ('',f8.1,'' km)'')') dxel,dxel*height_dim*1e-3
endif

if (print_node) then
   open(9,file='top_boundary.dat')
   ncircle=2*pi*r2/rmin*5
   dthcircle=2*pi/ncircle
  do i=1,ncircle+1
     write(9,*) r2*cos((i-1)*dthcircle),r2*sin((i-1)*dthcircle)
  enddo
endif
close(9)
if (print_node) then
   open(9,file='bottom_boundary.dat')
   ncircle=2*pi*r1/rmin*5
   dthcircle=2*pi/ncircle
  do i=1,ncircle+1
     write(9,*) r1*cos((i-1)*dthcircle),r2*sin((i-1)*dthcircle)
  enddo
  close(9)
endif

if (itracoption==2) then
   if (print_node) write(irefwr,*) 'PERROR(ieltogrid): needs updating for itractoption==',itracoption
   call instop
endif



end subroutine ieltogrid
