!   TRACDENS
!
!   Calculate density of tracers in nodal points
!   Parameters:
!        kmesh1   i  mesh definition for velocity
!        kprob1   i  problem definition for velocity
!        coormark i  array containing tracers
!        idens   o  vector of special structure containing density
!                   Can contain multiple density structures
!                   (e.g., for tracer ratio method)
!
!   PvK 080300/151105
!   Version updated in as in comprJP 
subroutine tracdens()
use sepmodulecomio
use sepmodulevecs
use sepmodulekmesh
use sepran_arrays
use geometry
use control
implicit none

if (cyl) then
   !call tracdens_cyl()
   if (print_node) write(irefwr,*) 'PERROR(tracdens): not yet suited for cyl'
   call instop
else 
   !write(irefwr,*) 'before tracdens_cart'
   !call instop
   call tracdens_cart()
   !if (print_node) call printdens(coor,ks(idens)%sol,ks(iarea)%sol,npoint)
   !write(irefwr,*) 'stop after printdens'
   !call instop
endif
end  subroutine tracdens

subroutine printdens(coor,dens,area,npoint)
use sepmodulecomio
use control
implicit none
integer :: npoint
real(kind=8) :: dens(*),area(*),coor(2,*)
integer :: i

if (.not.print_node) return
write(irefwr,'(4a15)') 'i','y','C','area'
do i=1,npoint,npoint/100
   write(irefwr,'(i15,3e15.7)') i,coor(2,i),dens(i),area(i)
enddo
!write(irefwr,*) 'stop in printdens'
!call instop

end subroutine printdens

!   TRACDENS_CART
!   Compute density of the tracers averaged in the nodal points
!   Cartesian version
subroutine tracdens_cart()
use sepmodulecomio
use sepmodulekmesh
use sepmodulevecs
use tracers
use coeff
use geometry
use mpetrac
use msper01
use sepran_arrays
use control
implicit none

integer :: ntot,i,j,k,ielh
real(kind=8) :: xm,ym,xn(6),yn(6),un(6),vn(6),u,v,done
integer :: nodno(6),imissed,isub,ikelmo,ifound
integer :: findsubelem,icorrect,index,nodlin(3)
real(kind=8) :: r,th,rl(3),xl(3),yl(3),sumsol,user(1),shapef(6)
real(kind=8) :: xi,eta,inv_Rb,phil(3),phiq(6),dphidksi(6),dphideta(6)
logical :: guess_first
integer :: itrac_start,itrac_end,itrac,ip,nl,nq
real(kind=8) :: area_per_trac
real(kind=8),parameter :: dp_three=3.0_8

inv_Rb=-1d0/Rb_local  ! negative sign to make composition C positive when it is heavy

if (itracoption.ne.1) then
   if (print_node) write(irefwr,*) 'PERROR(tracdens): not suited for markerchain'
   call instop
endif

!   Find extent of tracers for which to compute the density
!   We've implemented two cases so far:
!      1) JGR97 benchmark, absolute ratio
!      2) CH94 with duplicate tracers
ntot = 0
do idist=1,ndist
   ntot=ntot+ntrac(idist)
enddo
itrac_start = 1
itrac_end = ntrac(1)

if (periodic) then
   if (print_node) write(irefwr,*) 'PERROR(tracdens): not yet suited for periodic b.c.'
   call instop
endif

if (kelgrp/=6) then
   if (print_node) write(irefwr,*) 'PERROR(tracdens_cart):: current routine only suited for P2 triangles'
   if (print_node) write(irefwr,*) 'but kelgrp = ',kelgrp
   call instop
endif
 
! Make sure to have the up to date nodal point info in /pexcyc/
call pefilxy(2,kmesh1,kprob1,idens)



! check whether idens exists
if (idens <= 0) then
   if (print_node) write(irefwr,*) 'PERROR(tracdens): vector idens does not exist'
   call instop
endif

! We need to have a second array if we use the ratio method
if (ratio_method) then
   if (idens2<=0) then
      if (print_node) write(irefwr,*) 'PERROR(tracdens): vector idens 2 does not exist'
      if (print_node) write(irefwr,*) '  but ratio_method = ',ratio_method
      call instop
   endif
endif

!  Check whether array with area of surrounding elements has
!  been created
if (fieldC.and.iarea<=0) then
   if (print_node) write(irefwr,*) 'PERROR(tracdens): vector iarea does not exist'
   call instop
endif

!  Set density vector to zero
call setsoltozero(ks(idens)%sol,npoint)
if (ratio_method) call setsoltozero(ks(idens2)%sol,npoint)

!     write(irefwr,*) 'volume, ntot: ',volume,ntot

guess_first=.false.
imissed = 0
ifound=0
nl=3
nq=6
do i=itrac_start,itrac_end
   if (ibuoy_trac==2) then
      xm = tracer(i)%xnew
      ym = tracer(i)%ynew
   else
      xm = tracer(i)%x
      ym = tracer(i)%y
   endif
!  Find the element in which (xm,ym) lies
!  Look it up in the table 
   call pedetel(1,xm,ym,ielh)
   call sper01(kmeshc,coor,nodno,xn,yn,ielh)
!  NB Subdivision is only valid when area array is created using P1 triangles
!  Now that we have the right element, subdivide it into
!  four linear ones and find correct subtriangle
   ! find subelement
   isub = findsubelem(xn,yn,xm,ym,xl,yl,rl,nodno,nodlin)
   call tracdens01(ks(idens)%sol,tracer(i)%density,nodlin,xl,yl,xm,ym,rl,3)


   !call getshape_quad_xy(xm,ym,xn,yn,phil,phiq,dphidksi,dphideta,nl,nq)
   !call addtodensquad(ks(idens)%sol,tracer(i)%density,nodno,phiq)
!  if (ym<0.1834e-2) then
!     write(irefwr,'(''ielh: '',i5,2f12.3,6i5:)') ielh,xm,ym,nodno(1:6)
!  endif


enddo


if (ratio_method.and.fieldC) then
   !if (compress) then 
   !   write(irefwr,*) 'PERROR(tracdens): need to figure out how to do ratio method with compress'
   !   call instop 
   !endif
!  Set up density of the neutral tracers
   ip=ntrac(1)
   do i=1,ntrac(2)
      if (ibuoy_trac==2) then
         xm = tracer(ip+i)%xnew
         ym = tracer(ip+i)%ynew
      else
         xm = tracer(ip+i)%x
         ym = tracer(ip+i)%y
      endif
!     Find the element in which (xm,ym) lies
!     Look it up in the table
      call pedetel(1,xm,ym,ielh)
      call sper01(kmeshc,coor,nodno,xn,yn,ielh)

!     Now that we have the right element, subdivide it into
!     four linear ones and find correct subtriangle
      isub = findsubelem(xn,yn,xm,ym,xl,yl,rl,nodno,nodlin)
   
      call tracdens01(ks(idens2)%sol,tracer(i)%density,nodlin,xl,yl,xm,ym,rl,3)
   enddo
endif
      
if (bilinearC.and.ratio_method) then
   !do i=1,npoint,100
   !    write(irefwr,'(''dens   : '',i5,2f12.7)') i,buffr(ipdens+i-1),buffr(ipdens2+i-1)
   !enddo
   call tracdens06(ks(idens)%sol,ks(idens2)%sol,npoint)
   return

   !do i=1,npoint,100
   !    write(irefwr,'(''densr  : '',i5,f12.7)') i,buffr(ipdens+i-1)
   !enddo
   !call instop
endif

! normalize by integrated surface area (Tackley and King, eq. 9)
if (volume_per_tracer(1) <= 1e-12) then
   if (print_node) write(irefwr,*) 'PERROR(tracdens): volume_per_tracer(1) < 1e-12'
   if (print_node) write(irefwr,*) 'probably means it is zero: ',volume_per_tracer(1)
   call instop
endif
! shouldn't have to multiply by area of tracer as this info is already in tracer%density
!call tracdens05(ks(idens)%sol,ks(iarea)%sol,3d0*inv_Rb,npoint,truncateC)

! sepran 1020: we'll multiply with Rb_local in pefilbuoy()
!write(irefwr,*) 'idens, iarea: ',idens,iarea,dp_three,npoint,truncateC
call tracdens05(ks(idens)%sol,ks(iarea)%sol,dp_three,npoint,truncateC)
if (ratio_method) then
   if (volume_per_tracer(2) <= 1e-12) then
     if (print_node) write(irefwr,*) 'PERROR(tracdens): volume_per_tracer(2) < 1e-12'
     if (print_node) write(irefwr,*) 'probably means it is zero: ',volume_per_tracer(2)
     call instop
   endif
   call tracdens05(ks(idens2)%sol,ks(iarea)%sol,dp_three,npoint,truncateC)
!  find the ratio (Tackley and King, eq. 10)
   call tracdens06(ks(idens)%sol,ks(idens2)%sol,npoint)
!  write(irefwr,*) 'sumdensratio: ',sum(buffr(ipdens:ipdens+npoint-1))
endif

!if (compress) call multiplydensbyrhobar(ks(idens)%sol,npoint)


end subroutine tracdens_cart

subroutine multiplydensbyrhobar(dens,coor,npoint)
implicit none
integer :: npoint
real(kind=8) :: dens(*),coor(2,*)
real(kind=8) :: funccf,x,y
integer :: i

do i=1,npoint
   dens(i)=dens(i)*funccf(3,coor(1,i),coor(2,i),coor(2,i))
enddo

end subroutine multiplydensbyrhobar
