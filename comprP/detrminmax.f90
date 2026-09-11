! For a cylindrical mesh with inner boundary IBOT and
! outer boundary ITOP, determine the minimum/maximum radius
! of any point on the curve. 
! 
! PvK 970414
! PvK 020304
subroutine detrminmax(kmesh,kprob,isol)
use sepmodulekmesh
use tracers
use convparam
use coeff
use geometry
use control
implicit none
integer :: kmesh,kprob,isol
integer :: NUM,NC
parameter(NUM=10000,NC=10)
real(kind=8) :: funcx(NUM),funcy(NUM),x1,y1,x2,y2,r,th,x,y
integer :: i,npcurv
integer :: len_nodals
parameter(len_nodals=10000)
integer :: nodals(len_nodals),ielem,ipn
real(kind=8) :: xn(6),yn(6),rmax,rmin
real(kind=8) :: ds,dsbot_min,dstop_min,dstop_av,dsbot_av
integer :: nodno,icurvs(NC),numcurvs,num_nodals(NC),ip_nodals(NC) 
integer :: ielh

funcx(1)=NUM
funcy(1)=NUM

icurvs(1)=0
icurvs(2)=ibottom
call compcr(-1,kmesh,kprob,isol,0,icurvs,funcx,funcy)

npcurv = nint(funcx(5)/2)
radius_min=r1
do i=1,npcurv
   x = funcx(4+2*i)
   y = funcx(5+2*i)
   r = sqrt(x*x+y*y)
   radius_min= max(r,radius_min)
enddo
radius_min = radius_min+1d-6

dsbot_min=1e3
dsbot_av=0
do i=1,npcurv-1
   x1 = funcx(4+2*i)
   y1 = funcx(5+2*i)
   x2 = funcx(6+2*i)
   y2 = funcx(7+2*i)
   ds = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
   dsbot_min = min(dsbot_min,ds)
   dsbot_av  = dsbot_av + ds
enddo
dsbot_av = dsbot_av/(npcurv-1)

icurvs(1)=0
icurvs(2)=itop
call compcr(-1,kmesh,kprob,isol,0,icurvs,funcx,funcy)

npcurv = nint(funcx(5)/2)
radius_max= 0
do i=1,npcurv
   x = funcx(4+2*i)
   y = funcx(5+2*i)
   r = sqrt(x*x+y*y)
   radius_max = max(r,radius_max)
enddo
radius_max = radius_max-1d-6

dstop_min=1e3
dstop_av=0
do i=1,npcurv-1
   x1 = funcx(4+2*i)
   y1 = funcx(5+2*i)
   x2 = funcx(6+2*i)
   y2 = funcx(7+2*i)
   ds = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
   dstop_min = min(dstop_min,ds)
   dstop_av  = dstop_av + ds
enddo
dstop_av = dstop_av/(npcurv-1)

! Check which elements are curved
if (print_node) write(irefwr,*) 'nelem, NELEM_TOPO: ',nelem,NELEM_TOPO
if (nelem.gt.NELEM_TOPO) then
  if (print_node) then
     write(irefwr,*) 'PERROR(detrminmax): nelem>NELEM_TOPO'
     write(irefwr,*) 'increase NELEM_TOPO in elem_topo.inc'
     write(irefwr,*) 'to at least ',nelem
     write(irefwr,*) 'and recompile'
  endif
  call instop
endif

curved_elem(1:nelem)=.false.
if (check_on_curved_elements) call check_curved_elements(nelem,npoint,kmeshc,coor)



! We will keep track of top and bottom row of elements by raising a 
! flag in the corresponding element numbers in boundary_elem(NELEM)
! First, set all to false
do i=1,nelem
   boundary_elem(i)=.false.
enddo

! Get topology of top and bottom row of elements. Save for more
! accurate particle tracking.
numcurvs=2
icurvs(1)=itop
icurvs(2)=ibottom
call det_nods_boundary(kmeshm,ncurvs,numcurvs,icurvs,nodals,len_nodals,num_nodals,ip_nodals) 
! KMESH part D contains the elements associated with nodal points
call det_elems_boundary(nodals,num_nodals,kmeshd,npoint,nelem,numcurvs,ip_nodals)

if (print_node) then 
   write(irefwr,*) 'top boundary has # points    : ',num_nodals(1)
   ! write(irefwr,'(100i6,$)') (nodals(i),i=1,num_nodals(1))
   ! write(irefwr,*) 
   write(irefwr,*) 'bottom boundary has # points : ',num_nodals(2)
   ! ipn = ip_nodals(2)
   ! write(irefwr,'(100i6,$)') (nodals(i+ipn-1),i=1,num_nodals(2))
endif

! Check how large the elements at top and bottom are and
! figure out the threshold at which the tracing switches
! to cylindrical coordinates (see predcoort and mark4c).

! PvK May 2010
! Something goes wrong here with large number of elements
! the do loop is not terminated correctly leading to a seg fault
! once boundary_elem(NELEM_TOPO+1) is accessed. 
! the algorithm to find the boundary elements also looks wrong (finds too few 
! Let's just skip this for now and just assume that we switch to cylindrical
! tracing within 5 km from the boundary.

rtop_max = r2
rtop_min = r2 - 16.2d0/height_dim
rbot_min = r1
rbot_max = r1 + 16.2d0/height_dim

! Need to check on JGR97 tracing benchmark
! something wrong in cylindrical tracing in movetracers4
!rtop_threshold=1000
!rbot_threshold=0
! rtop_threshold is the radius above which we will use
! cylindrical coordinates for particle tracing
! the line below makes sense only if 
!rtop_threshold = rtop_min+(rtop_max-rtop_min)*0.5d0
! rbot_threshold is the radius above which we will use
! cylindrical coordinates for particle tracing
!rbot_threshold = rbot_min+(rbot_max-rbot_min)*0.5d0
rbot_threshold=rbot_max
rtop_threshold=rtop_min
if (print_node) then
   write(irefwr,'(''Min/max radius of elements of top element row estimated   :'', 2f12.7)') rtop_max,rtop_min
   write(irefwr,'(''     Threshold used for switch to cylindrical tracing  : '',f12.7)') rtop_threshold
   write(irefwr,'(''Min/max radius of elements of bottom element row :'', 2f12.7)') rbot_max,rbot_min
   write(irefwr,'(''     Threshold used for switch to cylindrical tracing  : '',f12.7)') rbot_threshold

   write(irefwr,'(''     Minimum size of elements along top boundary '',f12.7,'' ('',f12.7,'' km)'')') &
      & dstop_min,dstop_min*height_dim*1e-3
   write(irefwr,'(''     Average size of elements along top boundary '',f12.7,'' ('',f12.7,'' km)'')') & 
      & dstop_av,dstop_av*height_dim*1e-3
   write(irefwr,'(''     Minimum size of elements along bottom boundary '',f12.7,'' ('',f12.7,'' km)'')') & 
      & dsbot_min,dsbot_min*height_dim*1e-3
   write(irefwr,'(''     Average size of elements along bottom boundary '',f12.7,'' ('',f12.7,'' km)'')') & 
      & dsbot_av,dsbot_av*height_dim*1e-3
endif

end subroutine detrminmax

! GET_ELEM_COORDS
! Get nodal point numbers and coordinates for element ielem
! PvK 020304
subroutine get_elem_coords(ielem,kmeshc,coor,nodno,xn,yn)
implicit none
integer :: ielem 
integer :: kmeshc(*),nodno(*)
real(kind=8) :: coor(2,*),xn(*),yn(*)
integer :: i,ip

ip = (ielem-1)*6
do i=1,6
   nodno(i) = kmeshc(ip+i)
   xn(i)    = coor(1,nodno(i))
   yn(i)    = coor(2,nodno(i))
enddo

end subroutine get_elem_coords
 
! DET_ELEMS_BOUNDARY
! Determine the elements that are associated with the
! nodal points in nodals. Store in /elem_topo/ boundary_elems
! 
! nodals       i   nodal point numbers associated with curves
! num_nodals   i   array containing number of nodal points per curve
! kmeshd       i   KMESH part D (see Sepran PG)
! npoint       i   number of points in mesh
! nelem        i   number of elements in mesh
! numcurvs     i   number of curves associated with nodals
! ip_nodals    i   starting infomation for each curve
! 
! PvK 020304
subroutine det_elems_boundary(nodals,num_nodals,kmeshd,npoint,nelem,numcurvs,ip_nodals)
use sepmodulecomio
use control
use geometry
implicit none
integer :: num_nodals(*),npoint,nelem,nodals(*),kmeshd(*)
integer :: ip_nodals(*),numcurvs
integer :: i,nn,ne,ip,id2,ielem,j,ic

! point to kmeshd part 2
id2 = npoint+1
! ip keeps track of the number of elements found
! loop over nodal points at this boundary
do ic=1,numcurvs
   ip = 0
   do i=1,num_nodals(ic)
      ip=ip+1
!     nodal point number         
      nn = nodals(ip)
!     associated number of elements
      ne = kmeshd(nn+1)-kmeshd(nn)
      do j=1,ne
         ielem = kmeshd(id2+kmeshd(nn)+j)
         boundary_elem(ielem)=.true.
      enddo
   enddo
   if (print_node) write(irefwr,*) 'found ',ip,' elements along boundary curve ',ic
enddo

end subroutine det_elems_boundary

! DET_NODS_BOUNDARY
! Determine the nodal points associated with the the curves 
! icurvs.  Store in nodals(num_nodals).
! 
! kmeshm     i    kmesh part m (see Sepran PG)
! ncurvs     i    Number of curves stored in icurvs
! icurvs     i    array containing curve numbers
! nodals     o    nodal point numbers for all curves in icurvs
! num_nodals o    array containing number of nodal points in each curve
! 
! PvK 020304
subroutine det_nods_boundary(kmeshm,ncurvs,numcurvs,icurvs,nodals,len,num_nodals,ip_nodals)
use sepmodulecomio
use control
implicit none
integer :: kmeshm(*),ncurvs,icurvs(*),numcurvs,len,nodals(len)
integer :: ninner,ip3,ip4,len4,num_nodals(*),ip
integer :: i,ic,ip_nodals(*),ipn,num_nodals_tot,icurv

ninner = kmeshm(1)
ip3 = kmeshm(2)
ip4 = kmeshm(3)
len4 = kmeshm(4)
ip_nodals(1)=0
num_nodals_tot=0
ipn = 1

do ic=1,numcurvs
   ip_nodals(ic) = ipn
   icurv = icurvs(ic)
!  total number of points along ICURV
   num_nodals(ic) = kmeshm(4+icurv+1)-kmeshm(4+icurv)
!  position of first point of ICURV in kmeshm
   ip = ip3+kmeshm(4+icurv)-1
!  check to make sure the total number of nodal points 
!  doesn't exceed declared length of nodals(len)
   num_nodals_tot = num_nodals_tot + num_nodals(ic)
   if (num_nodals_tot.gt.len) then
      if (print_node) then
        write(irefwr,*) 'PERROR(det_topo_boundary): num_nodals>len'
        write(irefwr,*) 'num_nodals: ',num_nodals_tot
        write(irefwr,*) 'len       : ',len
      endif
      call instop
   endif
   do i=1,num_nodals(ic)
      nodals(i+ipn-1) = kmeshm(ip+i-1)
   enddo
   ipn = ipn + num_nodals(ic)
enddo

end subroutine det_nods_boundary

! CHECK_CURVED_ELEMENTS
! Check how many elements are curved. Set flags in 
! /elem_topo/ curved_elem(NELEM)
! 
! PvK 020304
subroutine check_curved_elements(nelem,npoint,nodno,coor)
use sepmodulecomio
use control
use geometry
implicit none
integer :: nelem,npoint,nodno(*)
real(kind=8) :: coor(2,*),x(6),y(6)
integer :: ielem,inpelm,i,icurved,ip
real(kind=8) :: dx,dy,fx,fy,epsh
logical :: curved

epsh = 1d-8
inpelm = 6
icurved = 0
ip=0
do ielem=1,nelem
   curved=.false.
   do i=1,inpelm
      x(i) = coor(1,nodno(ip+i))
      y(i) = coor(2,nodno(ip+i))
   enddo
   ip = ip+inpelm

   fx = abs( x(1)-2*x(2)+x(3) )
   fy = abs( y(1)-2*y(2)+y(3) )
   dx = abs(x(1)-x(3))
   dy = abs(y(1)-y(3))
   if (max(fx,fy).gt.max(dx,dy)*epsh) curved=.true.
   fx = abs( x(3)-2*x(4)+x(5) )
   fy = abs( y(3)-2*y(4)+y(5) )
   dx = abs(x(3)-x(5))
   dy = abs(y(3)-y(5))
   if (max(fx,fy).gt.max(dx,dy)*epsh) curved=.true.
   fx = abs( x(5)-2*x(6)+x(1) )
   fy = abs( y(5)-2*y(6)+y(1) )
   dx = abs(x(5)-x(1))
   dy = abs(y(5)-y(1))
   if (max(fx,fy).gt.max(dx,dy)*epsh) curved=.true.

   if (curved) then
      curved_elem(ielem) = .true.
      icurved = icurved + 1
   else
      curved_elem(ielem) = .false.
   endif
enddo

if (print_node) write(irefwr,*) 'Percentage of curved elements: ',icurved/nelem*100,' out of  ',nelem

  
end subroutine check_curved_elements
      
