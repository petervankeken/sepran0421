!************************************************************
!   TOGMT
!   Interpolate solution to regular grid for output with GMT
!   PvK 012809
!************************************************************
subroutine toGMT(jchoice,rlam,nres,fname)
use coeff
use sepran_arrays
use geometry
use control
implicit none
character(len=*) :: fname
integer :: jchoice,nres
real(kind=8) :: rlam

if (.not.print_node) return

if (cyl) then
   !write(irefwr,*) 'PERROR(toGMT): needs work for cyl'
   call toGMT_cyl(jchoice,rlam,nres,fname)
else if (.not.axi) then
   call toGMT_cart(jchoice,rlam,nres,fname)
endif

end

   
subroutine toGMT_cart(jchoice,rlam,nres,fname)
use sepmodulecomio
use eos
use convparam
use control
use sepran_arrays
implicit none
character(len=*) fname
integer :: jchoice,nres,nresx,nresy
real(kind=8) :: rlam
integer :: NXH,NYH
character(len=80) fname2
real(kind=8) :: x,y,xmax,ymax,average,dtheta,dr,dx,dy,yd,zd
integer :: iinmap(10),ncoor,i,j,ip,ichoice
integer :: nunks,ndim,allocate_status,map(5)=0
logical :: do_averages,first=.true.
real(kind=8),dimension(:),allocatable :: coorint2,solint2
save map,iinmap,first

do_averages = jchoice.ge.10
ichoice = jchoice - (jchoice/10)*10
!     *** figure out x/ymax and size of interpolation grid
xmax = rlam  ! rlam is aspect ratio determined from mesh
ymax = 1d0
!     *** set up coordinates of interpolation grid
dx = 1d0/nres
dy = dx
nresx=nint(rlam*nres+1)
nresy=nres+1
ncoor = nresx*nresy
if (.not.allocated(solint2)) then
   allocate(solint2(ncoor),stat=allocate_status)
   if (allocate_status /=0) then
      write(irefwr,*) 'PERROR(toGMT_cyl): could not allocate solint2: ',allocate_status
      call instop
   endif
endif
if (.not.allocated(coorint2)) then
   allocate(coorint2(2*ncoor),stat=allocate_status)
   if (allocate_status /=0) then
      write(irefwr,*) 'PERROR(toGMT_cyl): could not allocate coor: ',allocate_status
      call instop
   endif
endif
write(irefwr,*) 'ncoor: ',ncoor,nresx,nresy

ip=0
do j=1,nresy
   y=(j-1)*dy
   do i=1,nresx
      x = (i-1)*dx
      ip=ip+1
      coorint2(2*ip-1)=x
      coorint2(2*ip)=y
   enddo
enddo
open(11,file='GMT/coor.dat')
do ip=1,ncoor
      write(11,*) coorint2(2*ip-1),coorint2(2*ip)
enddo
close(11)


iinmap(1)=2
if (ichoice == 0) then
   iinmap(2)=1
else
   iinmap(2)=2
endif
!iinmap(2)=0
!map=0
nunks=1
ndim=2
!write(irefwr,*) 'map=',map(1:5)
!write(irefwr,*) 'isol2: ',isol2,kmesh1,kprob1
!write(irefwr,*) coorint2(1),coorint2(2)
!write(irefwr,*) nunks,ncoor,ndim
call intcoor(kmesh1,kprob1,isol2,solint2,coorint2,nunks,ncoor,ndim,iinmap,map)
!write(irefwr,*) 'map=',map(1:5)
!write(irefwr,*) 'isol2: ',isol2,kmesh1,kprob1
!write(irefwr,*) coorint2(1),coorint2(2)
!write(irefwr,*) nunks,ncoor,ndim

fname2='GMT/'//fname
open(11,file=fname2)
do i=1,ncoor
   write(11,*) solint2(i)
enddo
close(11)

if (cartplume.and.first) then
   first=.false.
   open(11,file='GMT/cartplume_profile.dat')
   ip=0
   do j=1,nresy
      ip=ip+nresx
      yd=coorint2(2*ip)*height_dim
      zd=(1-coorint2(2*ip))*height_dim
      write(11,'(5e15.7)') coorint2(2*ip),solint2(ip),get_Tbar_dim(zd)-Ts_dimK,get_alpha_dim(zd,rho_dim)/alpha_dim,&
        & get_rhobar_dim(zd)/rho_dim
   enddo
   close(11)
endif
      

if (do_averages) then
   fname2 = 'GMT/AVG_'//fname
   open(11,file=fname2)
   ip=0
   do j=1,nresy
      average=0d0
      do i=1,nresx
         ip=ip+1
         average = average + solint2(ip)
      enddo
      write(11,*) average/nresx,(j-1)*1d0/nresy
   enddo
 
   close(11)
endif
   
 
!if (allocated(coorint2)) deallocate(coorint2)
!if (allocated(solint2)) deallocate(solint2)

return
end subroutine toGMT_cart

subroutine toGMT_cyl(jchoice,rlam,nres,fname)
use sepmodulecomio
use coeff
use sepran_arrays
use geometry
implicit none
character(len=*) :: fname
integer :: jchoice,nres,nresx,nresy
real(kind=8) :: rlam
character(len=80) :: fname2
real(kind=8),dimension(:),allocatable :: coorint2,solint2
real(kind=8) :: x,y,xmax,ymax,average,dtheta,dr,r,theta,dx,dy,ratio_of_full,volume_full
integer :: iinmap(10),ncoor,i,j,ip,ichoice,map(5)=0
integer :: nunks,ndim,allocate_status
logical :: do_averages
save iinmap,map

dr=1d0/nres
NY=nres+1
NX=pi*(r1+r2)*nres*volume/(4.0_8/3*pi*(r2*r2-r1*r1)) ! make sure to account for geometries other than 'full'
volume_full=4.0_8/3*pi*r2*r2-4.0_8/3*pi*r1*r1
volume_full=volume_full/4*3
ratio_of_full=volume/volume_full
dtheta=2*pi/(NX-1)*ratio_of_full
write(irefwr,*) 'ratio: ',ratio_of_full,dtheta,(NX-1)*dtheta*360/(2*pi)
write(irefwr,*) 'volume: ',volume,volume_full,volume/volume_full

do_averages = jchoice.ge.10
ichoice = jchoice - (jchoice/10)*10

nresx=NX  
nresy=NY
ncoor = (nresx-1)*nresy  ! skip 1 in theta to avoid overlap at north pole
write(irefwr,*) 'nres: ',nresx,nresy,ncoor
if (.not.allocated(solint2)) then
   allocate(solint2(ncoor),stat=allocate_status)
   if (allocate_status /=0) then
      write(irefwr,*) 'PERROR(toGMT_cyl): could not allocate solint: ',allocate_status
      call instop
   endif
   allocate(coorint2(2*ncoor),stat=allocate_status)
   if (allocate_status /=0) then
      write(irefwr,*) 'PERROR(toGMT_cyl): could not allocate coor: ',allocate_status
      call instop
   endif
endif

ip=0
do j=1,nresy
   r=r1+(j-1)*dr
   do i=1,nresx-1 ! skip 1
      theta = (i-1)*dtheta
      ip=ip+1
      coorint2(2*ip-1)=r*sin(theta)
      coorint2(2*ip)=r*cos(theta)
   enddo
enddo
open(11,file='GMT/coor.dat')
do ip=1,ncoor
      write(11,*) coorint2(2*ip-1),coorint2(2*ip)
enddo
close(11)
 

iinmap(1)=0
if (ichoice == 0) then
   iinmap(2)=0
else
   iinmap(2)=0
endif
nunks=1
ndim=2
call intcoor(kmesh1,kprob1,isol2,solint2,coorint2,nunks,ncoor, ndim,iinmap,map)

fname2='GMT/'//fname
open(11,file=fname2)
do i=1,ncoor
   write(11,*) solint2(i)
enddo
close(11)

if (do_averages) then
   fname2 = 'GMT/AVG_'//fname
   open(11,file=fname2)
   ip=0
   do j=1,nresy
      average=0d0
      do i=1,nresx-1
         ip=ip+1
         average = average + solint2(ip)
      enddo
      write(11,*) average/nresx,r1+(j-1)*dr
   enddo
 
   close(11)
endif

!if (allocated(coorint2)) deallocate(coorint2)
!if (allocated(solint2)) deallocate(solint2)

end subroutine toGMT_cyl
