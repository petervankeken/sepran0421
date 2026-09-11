! INITTRAC_CART
subroutine use_tracers_from_nate_or_cian
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
integer :: nytractot
real(kind=8),dimension(:),allocatable :: htrac,ystretch,zstretch
logical my_turn

if (use_tracers_from_nate) then
      !tracers are preallocated
      if (ntrac_allocated<ntracers_from_nate) then
         if (print_node) write(irefwr,*) 'PERROR(inittrac_cart): ntrac_allocated is too small: ',ntrac_allocated,ntracers_from_nate
         call instop
      endif
      !call allocate_tracers(ntrac_allocated,1,'inittrac')
      open(9,file='black.dat')
      do i=1,ntracers_from_nate
         read(9,*) tracer(i)%x,tracer(i)%y
      enddo
      ip1=0
      idist=1
      ndist=1
      ntrac(idist)=ntracers_from_nate
      volume_per_tracer(1)=0.2*rlampix/ntracers_from_nate
      tracer(ip1+1:ip1+ntrac(idist))%density = -Rb_local*volume_per_tracer(idist)
      if (print_node) write(irefwr,*) 'ntracers_from_nate & density = ',ntracers_from_nate,-Rb_local*volume_per_tracer(idist)
      !write(irefwr,*) 'stopping in inittrac_cart'
      !call instop
      return

else if (use_tracers_from_cian) then
      !tracers are preallocated
      if (ntrac_allocated<ntracers_from_cian) then
         if (print_node) write(irefwr,*) 'PERROR(inittrac_cart): ntrac_allocated is too small: ',ntrac_allocated,ntracers_from_cian
         call instop
      endif
      open(9,file='cian.dat')
      do i=1,ntracers_from_cian
         read(9,*) tracer(i)%x,tracer(i)%y
      enddo
      ip1=0
      idist=1
      ndist=1
      ntrac(idist)=ntracers_from_cian
      volume_per_tracer(1)=0.2*rlampix/ntracers_from_cian
      if (ibench_type==4.and.nint(100*Di)==50) then
         volume_per_tracer(1)=0.31371856/ntracers_from_cian
         if (print_node) write(irefwr,*) 'cian volume per tracer: ',volume_per_tracer(1)
      else
         volume_per_tracer(1)=0.001_8*0.001_8
      endif
      tracer(ip1+1:ip1+ntrac(idist))%density = -Rb_local*volume_per_tracer(idist)
      if (print_node) write(irefwr,*) 'ntracers_from_cian & density = ',ntracers_from_cian,-Rb_local*volume_per_tracer(idist)
      return
   
endif

end subroutine use_tracers_from_nate_or_cian
