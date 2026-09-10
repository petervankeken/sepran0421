! Find surface mobility and plateness following Tackley (2000a)
! Mobility is rms velocity of surface over total rms velocity
! Plateness is a bit arbitrary. First, find the area where > 80% of average surface strainrate occurs (f_80)
! Then compute plateness P =  1 - f_80/0.6
! The 0.6 is the arbitrary bit.
subroutine find_mobility_plateness_cart(vrms,vsurf_rms,plate_mobility,plateness)
use sepmoduleoldrouts
use sepmodulecomio  ! provides lu's for reading & writing
use sepmodulecpack  ! parallel computing interface
use sepmodulekmesh
use sepmodulevecs
use sepran_arrays ! access to kmesh, kprob, isol, etc. and interfaces to subroutines
!use sepran_interface
use geometry
use control
use coeff
use convparam
!use dtm_elem  ! controls building of pressure mass matrix
implicit none
real(kind=8),intent(in) :: vrms,vsurf_rms
real(kind=8),intent(out) :: plate_mobility,plateness
integer :: ichoice,icurvs(10),ncoord,i,ndof,allocate_status
character(len=80) :: fname
real(kind=8) :: u,v,u2v2,dx,surface_length,IIsurf_int,secinv,f80,secinv1,secinv2
real(kind=8),dimension(:),allocatable :: funcx,funcy


if (.not.allocated(funcx)) then
   allocate(funcx(5+2*npoint),stat=allocate_status)
   if (allocate_status /= 0) then
      if (print_node) write(irefwr,*) 'PERROR(find_mobility_plateness_cart) trouble allocation funcx: ',allocate_status
      call instop
   endif
endif
if (.not.allocated(funcy)) then
   allocate(funcy(5+2*npoint),stat=allocate_status)
   if (allocate_status /= 0) then
      if (print_node) write(irefwr,*) 'PERROR(find_mobility_plateness_cart) trouble allocation funcy: ',allocate_status
      call instop
   endif
endif
funcx(1)=1.0_8*(5+2*npoint)
funcy(1)=1.0_8*(5+2*npoint)
write(6,*) 'find_mobility_plateness: start'

plate_mobility=vsurf_rms/vrms
! Get the second invariant of the strainrate tensor at the top
ichoice=0
icurvs(1)=0
icurvs(2)=itop
call compcr(ichoice,kmesh,kprob,isecinv,0,icurvs,funcx,funcy)
ncoord=nint(funcx(5)/2)
ndof=nint(funcy(5))
if (ndof /= ncoord) then
   if (print_node) then
      write(irefwr,*) 'something weird in fluid_mobility_plateness_cart'
      write(irefwr,*) 'ndof, ncoord= ',ndof,ncoord
   endif
   call instop
endif
! Find integrated surface second invariant
IIsurf_int=0
! Trapezoid rule
! Deal with end points first
dx = abs(funcx(8)-funcx(6))
secinv=funcy(5+i)
IIsurf_int=IIsurf_int+secinv*dx*0.5
dx = abs(funcx(5+2*ncoord-1)-funcx(5+2*ncoord-3))
secinv=funcy(5+ncoord)
IIsurf_int=IIsurf_int+0.5*secinv*dx
surface_length=dx
! Interior points
do i=2,ncoord-1
   secinv=funcy(5+i)
   dx = abs(funcx(5+2*i+1)-funcx(5+2*i-1))
   IIsurf_int= IIsurf_int+ secinv*dx
   surface_length=surface_length+dx
enddo
plateness=0.0_8

write(fname,'(''GMT/IIsurf.dat'')') 
open(9,file=fname)
do i=1,ncoord
   write(9,'(2f15.7)') funcx(5+2*i-1),funcy(5+i)
enddo
close(9)
! for plotting line at 80% of the integrated strainrate
write(fname,'(''GMT/IIsurf80.dat'')')
open(9,file=fname)
write(9,'(2f15.7)') funcx(6),IIsurf_int*0.8
write(9,'(2f15.7)') funcx(5+2*ncoord-1),IIsurf_int*0.8
close(9)

! find f80 by checking each element
f80=0
secinv1=funcy(6)
do i=1,ncoord-1
   secinv2=funcy(5+i+1)
   if ((secinv2+secinv1)*0.5 > IIsurf_int) f80=f80+dx
   secinv1=secinv2
enddo
! rather arbitrary definition from Tackley 2000a
plateness = 1 - f80/0.6

write(irefwr,*) 'vrms,surf, vrms = ',vsurf_rms,vrms,plate_mobility,IIsurf_int,f80,plateness

if (allocated(funcx)) deallocate(funcx)
if (allocated(funcy)) deallocate(funcy)

end subroutine find_mobility_plateness_cart





