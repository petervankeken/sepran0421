subroutine tracdens05(dens,area,area_per_trac,npoint,truncateC)
implicit none
real(kind=8) :: dens(*),area(*),area_per_trac
integer :: npoint,i
logical :: truncateC

do i=1,npoint
   dens(i) = dens(i)*area_per_trac/area(i)
   if (truncateC) then
     dens(i) = min(dens(i),1d0)
     dens(i) = max(dens(i),0d0)
   endif
enddo
!write(irefwr,*) 'tracdens05: ',minval(dens(1:npoint)),maxval(dens(1:npoint))

end subroutine tracdens05

subroutine tracdens06(dens,dens2,npoint)
use coeff
implicit none
real(kind=8) :: dens(*),dens2(*),ratio
integer :: npoint,i

if (its_CH94) then
  ! dens  = dense eclogite tracers
  ! dens2 = harzburgite tracers with 7x the volume
! write(irefwr,*) 'dens : ',minval(dens(1:npoint)),maxval(dens(1:npoint)),dens(npoint-1)
! write(irefwr,*) 'dens2: ',minval(dens2(1:npoint)),maxval(dens2(1:npoint)),dens2(npoint-1)
  do i=1,npoint
     if (dens(i)+dens2(i).gt.0) dens(i) = dens(i)/(dens(i)+7*dens2(i))
  enddo
!  write(irefwr,*) 'dens: ',minval(dens(1:npoint)),maxval(dens(1:npoint))
!  !call instop
  return
endif

if (ratio_method_overlay) then
  ! neutral tracers are box filling with same original density
  ! so maximum of N/(N+M)=0.5
  do i=1,npoint
     if (dens(i)+dens2(i).gt.0) dens(i) = 2.0_8*dens(i)/(dens(i)+dens2(i))
  enddo
else
! write(irefwr,*) 'dens : ',minval(dens(1:npoint)),maxval(dens(1:npoint))
! write(irefwr,*) 'dens2: ',minval(dens2(1:npoint)),maxval(dens2(1:npoint))
  do i=1,npoint
     if (dens(i)+dens2(i).gt.0) dens(i) = dens(i)/(dens(i)+dens2(i))
  enddo
! write(irefwr,*) 'dens: ',minval(dens(1:npoint)),maxval(dens(1:npoint))
! call instop
endif


end subroutine tracdens06
 
!   TRACERPIX_TOVEC
subroutine tracerpix_tovec(dens,coor,npoint)
use tracers
use coeff
implicit none
real(kind=8) :: dens(*),coor(2,*)
integer :: npoint,i,ip,ix,iy,ndens
real(kind=8) :: xm,ym,denstot

denstot = 0
ndens=0
do i=1,npoint
   xm = coor(1,i)
   ym = coor(2,i)
   ix = min(nint(xm/rlampix*nxpix+0.5),nxpix)
   iy = min(nint(ym*nypix+0.5),nypix)
   ip = (iy-1)*nxpix+ix
   dens(i) = pix(ip)
   denstot = denstot+dens(i)
   if (dens(i).gt.0.1) ndens=ndens+1
    
enddo

end subroutine tracerpix_tovec

!   FINDSUBELEM
!   Each quadratic element contains 4 linear subelements
!   Find the subelement number, based on the standard subdivision
!   that is used in nstmsh().
!   5
!   | \
!   |  \
!   |IV \
!   6----4
!   |\III/\
!   |I\ /II\
!   1---2---3
integer function findsubelem(xn,yn,xm,ym,xl,yl,rl,nodno,nodlin) 
use sepmodulecomio
use control
implicit none
real(kind=8) :: xn(6),yn(6),xm,ym
real(kind=8) :: xl(3),yl(3),a(3),b(3),c(3)
real(kind=8) :: epsh,rl(3),r
integer :: isub,nodlin(3),nodno(6),iprint
parameter(epsh=1d-4)
logical :: in_this_one
real(kind=8) :: area

iprint=0
isub=1

xl(1)=xn(1)
yl(1)=yn(1)
xl(2)=xn(2)
yl(2)=yn(2)
xl(3)=xn(6)
yl(3)=yn(6)
nodlin(1)=nodno(1)
nodlin(2)=nodno(2)
nodlin(3)=nodno(6)

call trilin(xl,yl,a,b,c,area)
rl(1) = a(1)+b(1)*xm+c(1)*ym
rl(2) = a(2)+b(2)*xm+c(2)*ym
rl(3) = a(3)+b(3)*xm+c(3)*ym
in_this_one = (rl(1).ge.-epsh).and.(rl(2).ge.-epsh).and.(rl(3).ge.-epsh)
if (iprint.eq.1.and.print_node) write(irefwr,'(i5,3f8.3,L10)') isub,rl(1),rl(2),rl(3),in_this_one
if (in_this_one) then
   findsubelem = isub
   return
endif

isub=2

xl(1)=xn(2)
yl(1)=yn(2)
xl(2)=xn(3)
yl(2)=yn(3)
xl(3)=xn(4)
yl(3)=yn(4)
nodlin(1)=nodno(2)
nodlin(2)=nodno(3)
nodlin(3)=nodno(4)

call trilin(xl,yl,a,b,c,area)
rl(1) = a(1)+b(1)*xm+c(1)*ym
rl(2) = a(2)+b(2)*xm+c(2)*ym
rl(3) = a(3)+b(3)*xm+c(3)*ym
in_this_one = (rl(1).ge.-epsh).and.(rl(2).ge.-epsh).and.(rl(3).ge.-epsh)
if (iprint.eq.1.and.print_node) write(irefwr,'(i5,3f8.3,L10)') isub,rl(1),rl(2),rl(3),in_this_one
if (in_this_one) then
   findsubelem = isub
   return
endif

isub=3

xl(1)=xn(4)
yl(1)=yn(4)
xl(2)=xn(5)
yl(2)=yn(5)
xl(3)=xn(6)
yl(3)=yn(6)
nodlin(1)=nodno(4)
nodlin(2)=nodno(5)
nodlin(3)=nodno(6)
call trilin(xl,yl,a,b,c,area)
rl(1) = a(1)+b(1)*xm+c(1)*ym
rl(2) = a(2)+b(2)*xm+c(2)*ym
rl(3) = a(3)+b(3)*xm+c(3)*ym
in_this_one = (rl(1).ge.-epsh).and.(rl(2).ge.-epsh).and.(rl(3).ge.-epsh)
if (iprint.eq.1.and.print_node) write(irefwr,'(i5,3f8.3,L10)') isub,rl(1),rl(2),rl(3),in_this_one
if (in_this_one) then
   findsubelem = isub
   return
endif

isub=4

xl(1)=xn(2)
yl(1)=yn(2)
xl(2)=xn(4)
yl(2)=yn(4)
xl(3)=xn(6)
yl(3)=yn(6)
nodlin(1)=nodno(2)
nodlin(2)=nodno(4)
nodlin(3)=nodno(6)

call trilin(xl,yl,a,b,c,area)
rl(1) = a(1)+b(1)*xm+c(1)*ym
rl(2) = a(2)+b(2)*xm+c(2)*ym
rl(3) = a(3)+b(3)*xm+c(3)*ym
in_this_one =  (rl(1).ge.-epsh).and.(rl(2).ge.-epsh).and.(rl(3).ge.-epsh)
if (iprint.eq.1.and.print_node) write(irefwr,'(i5,3f8.3,L10)') isub,rl(1),rl(2),rl(3),in_this_one
if (in_this_one) then
   findsubelem = isub
   return
endif
    
if (print_node) write(irefwr,*) 'PERROR(findsubelem): tracer is not in any'
if (print_node) write(irefwr,*) ' of the sub elements - assuming linear shapes.'
if (print_node) write(irefwr,*) ' This can happen when elements are curved.'
if (print_node) write(irefwr,*) ' Rewrite code for this situation.'
if (print_node) write(irefwr,*) xm,ym
if (print_node) write(irefwr,*) xn(1),yn(1),xn(3),yn(3),xn(5),yn(5)
  
end function findsubelem

!   TRACDENS01
!   Add subelement contribution to the correct nodal points
subroutine tracdens01(dens,densmark,nodno,xl,yl,xm,ym,phi,npelm)
use sepmodulecomio
use coeff
use control
implicit none
real(kind=8) :: dens(*),phi(*),densmark,xl(*),yl(*),xm,ym
integer :: nodno(*),npelm
integer :: k
real(kind=8) :: dx,dy,Sx,Sy,funccf

if (bilinearC) then
  if (compress) then
     if (print_node) write(irefwr,*) 'PERROR(tracdens01): not yet suited for bilinearC and compress'
     call instop
  endif
  
  ! map tracer onto nodal points using bilinear shape functions
  ! works only if element sides are parallel to x and y, of course
  !dx=max(abs(xl(1)-xl(2)),abs(xl(2)-xl(3)),abs(xl(3)-xl(1)))
  !dy=max(abs(yl(1)-yl(2)),abs(yl(2)-yl(3)),abs(yl(3)-yl(1)))
  dx=max(abs(xl(1)-xl(2)),abs(xl(2)-xl(3)))
  dy=max(abs(yl(1)-yl(2)),abs(yl(2)-yl(3)))
  !write(irefwr,'(''dx,dy: '',2f12.3,5x,3f12.3,5x,3f12.3)') dx,dy,xl(1:3),yl(1:3)
  do k=1,npelm
     Sx=xm-xl(k)
     Sy=ym-yl(k)
     Sx=abs(Sx/dx)
     Sy=abs(Sy/dy)
     Sx=max(0.0_8,1.0_8-Sx)
     Sy=max(0.0_8,1.0_8-Sy)
     dens(nodno(k))=dens(nodno(k)) + Sx*Sy
  enddo
else
  ! Use linear shapefunctions on triangle
! if (compress) then
!    do k=1,npelm
!         dens(nodno(k)) = dens(nodno(k)) + phi(k)*densmark*funccf(3,xm,ym,ym)
!    enddo
! else
  ! when forming C vector for field method do not multiply with rhobar
  ! in the case of compressibility - sepran multiplies buoyancy for with rhobar
     do k=1,npelm
          dens(nodno(k)) = dens(nodno(k)) + phi(k)*densmark
     enddo
! endif
endif

end subroutine tracdens01

!call addtodensquad(ks(idens)%sol,tracer(i)%density,nodno,phiq)
subroutine addtodensquad(dens,densmark,nodno,phiq)
implicit none
integer :: k
real(kind=8) :: dens(*),densmark,phiq(*)
integer :: nodno(*)

! quadratic shape function on 6 point triangle
do k=1,6
   dens(nodno(k))=dens(nodno(k))+phiq(k)*densmark
enddo

end subroutine addtodensquad

!   TRACDENS02
!   Initialize density vector
subroutine setsoltozero(dens,npoint)
implicit none
real(kind=8) :: dens(*)
integer :: npoint,i

do i=1,npoint
   dens(i)=0
enddo
  
return
end subroutine setsoltozero

!   TRACDENS03
!   Print density vector
subroutine printsol(dens,npoint)
use sepmodulecomio
use control
implicit none
real(kind=8) :: dens(*)
integer :: npoint,i,ip,irow

if (.not.print_node) return

do ip=1,npoint,54
   irow=ip/54+1
   write(irefwr,'(i5,21f4.1,:)') irow, (dens(i),i=ip,min(npoint,ip+54),4)
enddo
end subroutine printsol

!   TRACDENS04
!   Compute sum of density vector
real(kind=8) function sumsol(dens,npoint)
use sepmodulecomio
use control
implicit none
real(kind=8) :: dens(*)
integer :: npoint,i

sumsol=0
do i=1,npoint
   sumsol= sumsol+ dens(i)
enddo

end function sumsol

!   TRACERPIX_DENSITY
subroutine tracerpix_density(coormark,densmark)
use sepmodulecomio
use coeff
use tracers
use control
implicit none
real(kind=8) :: coormark(*),densmark(*)
integer :: ntot,itrac_start,itrac_end
integer :: i,ielh,ix,iy,ip,npixtot,ip1
real(kind=8) :: xm,ym,pixvolume,pixtot,dxpix_inv,dypix_inv
real(kind=8) :: average_tracer_per_pix,volume_layer
integer :: ipixtot

ntot = 0
do idist=1,ndist
   ntot = ntot+ntrac(idist)
enddo

! Find extent of tracers for which to compute the density
! We've implemented two cases so far:
!   1) JGR97 benchmark, absolute ratio
!   2) CH94 with duplicate tracers
if (ibench_type == 5) then
!  JGR 97 2 thermochemical benchmark
   itrac_start = 1
   itrac_end = ntot
endif
if (itrac_start.eq.0.or.itrac_end.eq.0) return

! Initialize pixel array
pix=0
ipix=0

! For each tracer, find corresponding pixel. Add value
! of density in tracer
dxpix_inv = 1d0/dxpix
dypix_inv = 1d0/dypix
do i=itrac_start,itrac_end
   xm = coormark(2*i-1)
   ym = coormark(2*i)
   ix = min(nint(xm*dxpix_inv + 0.5),nxpix)
   iy = min(nint(ym*dypix_inv + 0.5),nypix)
!  write(irefwr,*) 'tracer: ',i,xm,ym,ix,iy,densmark(i)
   ip = (iy-1)*nxpix+ix
   pix(ip)=pix(ip)+densmark(i)
!  keep track of number of tracers in ipix(1:nxpix*nypix)
   ipix(ip)=ipix(ip)+1
enddo

if (ibench_type == 5) then
   ! JGR97  2 thermochemical benchmark      
!  How many tracers do you need for C=1?
   average_tracer_per_pix=(itrac_end-itrac_start+1)/volume_dist(1)
   average_tracer_per_pix=average_tracer_per_pix/(nxpix*nypix)
   pixtot = 0
   ipixtot = 0
   ip=0
   do iy=1,nypix
      do ix=1,nxpix
         ip=ip+1
         pix(ip)=pix(ip)/average_tracer_per_pix
         pixtot = pixtot+pix(ip)
         ipixtot = ipixtot + ipix(ip)
      enddo
   enddo
   if (print_node) then
      write(irefwr,*) 'nxpix, nypix: ',nxpix,nypix,ntot,itrac_end-itrac_start+1,volume_dist(1)
      write(irefwr,*) 'Sum of pixel values: ', pixtot,ipixtot,average_tracer_per_pix
   endif
   return
endif

end subroutine tracerpix_density
