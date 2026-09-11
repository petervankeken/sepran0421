subroutine inittrac_cart_1(dxtrac,dytrac)
use sepmodulecomio
use convparam
use blobs
use mpetrac
use tracers
use coeff
use mparallel
use geometry
use control
#ifdef MPI
use mpi
#endif
implicit none
integer :: ip,i,j,ntot,iy,ix,itrac,nxtrac,nytrac,ip1
real(kind=8) ::  dxtrac,dytrac,x,y,local_volume
real(kind=8) ::  piw,ym,xm,yoffset,ymin,ymax,row_intvl,col_intvl
real(kind=8) ::  dx,xstart,xend,xhere,xb,yb,theta,xbl,dz
character(len=80) :: tname,fdummy
integer :: nxtrac_here,nyhere,ierr=0
real(kind=8) :: density_ratio,factor_a,sumhi,scale_sumhi,htop,mass_now,massgoal,htot,funccf,stretchcorrect
integer :: nytractot,iy1,iy2
real(kind=8),dimension(:),allocatable :: htrac,ystretch,zstretch
logical :: my_turn

! tracer method for normal cases (not benchmark type 4 or 101)

if (print_node) write(irefwr,*) 'PINFO(inittrac_cart_1)::  mpi_partrac,numprocs= ',mpi_partrac,numprocs
!  one dense layer without amplitude (for use in 2nd JGR97
!  benchmark or Tackley&King 2003 new benchmark
ip1=0
ndist = 1
if (ratio_method.and..not.ratio_method_overlay) then
   ! create second set of tracers that is in the complement of the domain
   ndist = 2
   x0tr(2)=x0tr(1)
   x1tr(2)=x1tr(1)
   y0tr(2)=y1tr(1) ! note second set of tracers sit on top of first set
   y1tr(2)=1.0_8
else if (ratio_method.and.ratio_method_overlay) then
   ! fill entire box with duplicate set
   ndist=2
   x0tr(2)=xcmin
   x1tr(2)=xcmax
   y0tr(2)=ycmin
   y1tr(2)=ycmax
endif
 
if (its_CH94.and.ratio_method) then
   ! note second set of layer overlaps with first
   ndist = 2
   x0tr(2)=x0tr(1)
   x1tr(2)=x1tr(1)
   y0tr(2)=y0tr(1)
   y1tr(2)=y1tr(1)
endif
        

! set up 1 or 2 distributions
!  nxtrac = rlampix/dr
nxtrac = int((x1tr(1)-x0tr(1))/dr)  ! rlampix/dr
dr = (x1tr(1)-x0tr(1))/nxtrac ! rlampix/nxtrac
dxtrac = dr
dytrac = dxtrac
eps_tracer_bound(1:ndist) = 0.5*dxtrac
ip = 0
ntot = 0

! Always use array ystretch (even for incompressible)
! Makes for slightly more work in initialization but leads to more compact and coherent code
nytractot=nint(4.0_8/dytrac)
write(irefwr,*) 'nytractot: ',nytractot,dytrac,dxtrac
allocate(htrac(nytractot),ystretch(nytractot),zstretch(nytractot))
call prepare_ystretch(htrac,zstretch,ystretch,nytractot,dytrac)

do idist=1,ndist
   if (compress.and.stretch_tracers) then
      iy1=-1
      iy2=-1
      do i=1,nytractot
         if (ystretch(i) >= y0tr(idist)+0.5*dytrac-1e-5.and.iy1<0) iy1=i
         if (ystretch(i) >= y1tr(idist)-0.5*dytrac+1e-5.and.iy2<0) iy2=i-1
      enddo
   else
      ymin = y0tr(idist)
      ymax = y1tr(idist)
      nytrac = int((ymax-ymin)/dytrac)
      iy1=1
      iy2=nytrac
   endif
   !if (mpi_partrac.and.numprocs>1) then
      ! replace this by striping over ix
      !
      !nxtrac_here = nxtrac/numprocs
      !dx = rlampix/numprocs
      !xstart = x0tr(idist)+myid*dx
      !xend   = x0tr(idist)+(myid+1)*dx
      !if (nxtrac_here * numprocs /= nxtrac) then
      !   PvK 020121: again, not sure why this needs to be enforced
      !   if (print_node) then
      !    write(irefwr,*) 'PERROR(inittrac): nxtrac_here*numprocs'
      !    write(irefwr,*) '    should be equal to nxtrac, but: '
      !    write(irefwr,*) '    nxtrac_here * numprocs = ', nxtrac_here*numprocs
      !    write(irefwr,*) '    nxtrac : ',nxtrac
      !   endif ! myid==0
      !   call instop
      !endif
   !else
   !   ! just one process to manage
   !   nxtrac_here = nxtrac
   !   dx = x1tr(idist)-x0tr(idist) ! rlampix
   !   xstart = x0tr(idist) ! 0
   !   xend = x1tr(idist) ! rlampix
   !endif
   !if (print_node) write(irefwr,'(''idist: '',3i10,4f12.3)') &
   !      & idist,nxtrac_here,nytrac,x0tr(idist),x1tr(idist),y0tr(idist),y1tr(idist) ! ,0d0,rlampix,ymin,ymax
   itrac=0
   if (print_node) open(9,file='firstrow.dat')
   do ix=1,nxtrac
#ifdef MPI
     my_turn=mod(ix+numprocs-1,numprocs)==myid    
     if (my_turn) then
#endif
        xhere = xstart+(ix-0.5)*dxtrac
        do iy=iy1,iy2
           itrac=itrac+1
           if (itrac>NTRACMAX) then
              write(irefwr,*) 'PERROR(inittrac_cart): ntot>NTRACMAX; myid=',myid
              write(irefwr,*) 'ntot     = ',ntot
              write(irefwr,*) 'NTRACMAX = ',NTRACMAX
              call instop
           endif
           tracer(itrac+ip)%x = xhere
           tracer(itrac+ip)%y = ystretch(iy) !ymin+(iy-0.5)*dytrac
           if (ix==1) write(9,*) iy,tracer(itrac+ip)%y
        enddo !iy
#ifdef MPI
     endif
#endif
   enddo !ix 
   close(9)
   ntrac(idist) = itrac
   ip = ip+itrac
   itrac = 0
enddo ! idist=1,ndist
ntot=0
do idist=1,ndist
   ntot=ntot+ntrac(idist)
enddo
end subroutine inittrac_cart_1
     
