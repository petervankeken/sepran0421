! INITTRAC_CART_1_BENCH
! Initialize tracers for Cartesian geometry specifically for the RT instability benchmark
subroutine inittrac_cart_1_bench4(dxtrac,dytrac)
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
integer :: ip,i,j,ntot,iy,ix,itrac,nxtrac,nytrac,ip1,ierr
real(kind=8) ::  dxtrac,dytrac,x,y,local_volume
real(kind=8) ::  piw,ym,xm,yoffset,ymin,ymax,row_intvl,col_intvl
real(kind=8) ::  dx,xstart,xend,xhere,xb,yb,theta,xbl,dz
character(len=80) :: tname,fdummy
integer :: nxtrac_here,nyhere
real(kind=8) :: density_ratio,factor_a,sumhi,scale_sumhi,htop,mass_now,massgoal,htot,funccf,stretchcorrect
integer :: nytractot
real(kind=8),dimension(:),allocatable :: htrac,ystretch,zstretch
logical :: my_turn

if (ainittr(1) <= 1e-7_8 .and. ibench_type==4) then
   if (print_node) then
      write(irefwr,*) 'PERROR(inittrac_cart): ainittr(1) < 0: ', ainittr(1)
      write(irefwr,*) 'which is inconsistent with ibench_type=4'
   endif
   call instop
endif
ip1=0
idist = 1
nxtrac = int(rlampix/dr)
if (weak_scaling_test) then
   nxtrac=nxtrac*sqrt(numprocs*1.0_8)
endif
dr = rlampix/nxtrac
dxtrac = dr
dytrac = dxtrac
eps_tracer_bound(1:ndist) = 0 ! 0.5*dxtrac
itrac=0
!if (mpi_parallel.and. (nxtrac/numprocs)*numprocs /= nxtrac) then
!   if (print_node) then
!      write(irefwr,*) 'PERROR(inittrac): tracers for RT benchmark '
!      write(irefwr,*) 'need to be set for parallel computing such '
!      write(irefwr,*) 'that nxtrac divides without remainder by '
!      write(irefwr,*) 'the number of processors, but here '
!      write(irefwr,*) '  nxtrac   = ',nxtrac
!      write(irefwr,*) '  numprocs = ',numprocs
!   endif
!   call instop
!endif
if (print_node) write(irefwr,'(''inittrac_cart: nxtrac= '',i10,f12.3,e15.7)') nxtrac,rlampix,dr

! For RT benchmark the compressible and incompressible formulations are slightly
! different due to the stack stretching in compressible specifically for this benchmark

if (compress) then
   nytractot=nint(2.0_8/dytrac) ! generous amount for general application
   allocate(htrac(nytractot),zstretch(nytractot),ystretch(nytractot))
   if (print_node) write(irefwr,*) 'allocate for nytractot'
   ! should we stretch at every x position for RB instability?????
   call prepare_ystretch(htrac,zstretch,ystretch,nytractot,dytrac)
#ifdef MPI
   write(fdummy,'(''xtracer.'',i3.3)') myid
   open(907,file=fdummy)
#endif

   ! now fill the tracer array by picking all stretched particles in a column below the interfac
   !do ix=myid+1,nxtrac,numprocs  ! stripe tracers in x direction
   do ix=1,nxtrac
#ifdef MPI
      my_turn=mod(ix+numprocs-1,numprocs)==myid     
      if (my_turn) then
#endif
        xm = (ix-0.5)*dxtrac
        ymax = y1tr(idist)+ainittr(idist)*cos(pi*xm/rlampix)
        iy=0
        do  ! iy
           iy=iy+1
           if (ystretch(iy)<=ymax.and.iy<=nyhere) then
              ! particle is still in dense layer
              ym=ystretch(iy)
              itrac = itrac+1
              if (itrac.gt.NTRACMAX) then   
                 write(irefwr,*) 'PERROR(inittrac_cart): ',myid
                 write(irefwr,*) '    itrac>NTRACMAX '
                 write(irefwr,*) '    itrac,NTRACMAX = ',itrac,NTRACMAX
                 call instop
              endif
              tracer(ip1+itrac)%x = xm
              tracer(ip1+itrac)%y = ym
           else
              exit
           endif
        enddo ! loop over y positions
#ifdef MPI
     endif
#endif
     if (ix==1.and.print_node.and..not.mpi_partrac) then 
        open(9,file='firstcolumn.dat')
        do i=1,itrac
           write(9,*) tracer(ip1+i)%x,tracer(ip1+i)%y
        enddo
        close(9)
        !write(irefwr,*) 'stop after firstcolumn' 
        !call instop
     endif
   enddo ! loop over x positions
#ifdef MPI
   close(907)
#endif
   ntrac(1)=itrac
   if (allocated(zstretch)) deallocate(zstretch)
   if (allocated(ystretch)) deallocate(ystretch)
   if (allocated(htrac)) deallocate(htrac)


else ! not compress

    ! first the set of dense tracers
    do ix=myid+1,nxtrac,numprocs  ! stripe tracers in x direction
       xm = (ix-0.5)*dxtrac
       ymax = y1tr(idist)+ainittr(idist)*cos(pi*xm/rlampix)
       ! make sure tracers fill the stack below the interface
       nytrac = int((ymax-0.5*dytrac)/dytrac + 1)
       do iy=nytrac,1,-1
          ym = (iy-0.5)*dytrac
          ! if (ix==1) write(irefwr,*) ym
          itrac = itrac+1
          if (itrac.gt.NTRACMAX) then   
             write(irefwr,*) 'PERROR(inittrac_cart): ',myid
             write(irefwr,*) '    itrac>NTRACMAX '
             write(irefwr,*) '    itrac,NTRACMAX = ',itrac,NTRACMAX
             call instop
          endif
          tracer(ip1+itrac)%x = xm
          tracer(ip1+itrac)%y = ym
       enddo ! iy
    enddo !ix
    ntrac(1)=itrac
    if (ratio_method) then
    ! set volume of distribution
    x0tr(2)=x0tr(1)
    x1tr(2)=x1tr(1)
    y0tr(2)=y1tr(1) ! note second set of tracers sit on top of first set
    y1tr(2)=1.0_8
    write(irefwr,*) 'neutral tracers'
    ip1=ip1+ntrac(1)
    ! now the neutral tracers
    idist=2 
    itrac=0
    do ix=myid+1,nxtrac,numprocs  ! stripe tracers in x direction
       xm = (ix-0.5)*dxtrac
       ymax = y1tr(1)+ainittr(1)*cos(pi*xm/rlampix)
       nytrac = int((1d0-ymax+0.5*dytrac)/dytrac )
       do iy=nytrac,1,-1
          ym = 1d0-(iy-0.5)*dytrac
          ! if (ix==1) write(irefwr,*) ym
          itrac = itrac+1
          if (ip1+itrac.gt.NTRACMAX) then   
             if (myid==0) then
                write(irefwr,*) 'PERROR(inittrac_cart): '
                write(irefwr,*) '    itrac>NTRACMAX '
                write(irefwr,*) '    itrac,NTRACMAX = ',itrac,NTRACMAX
             endif
             call instop
          endif
          tracer(ip1+itrac)%x = xm
          tracer(ip1+itrac)%y = ym
       enddo !iy
     enddo !ix
     ntrac(2)=itrac
     ndist=2
  endif ! tracer_ratio
endif ! compress

if (print_node) then
   write(irefwr,*) 'made RT initial condition'
endif

end subroutine inittrac_cart_1_bench4
