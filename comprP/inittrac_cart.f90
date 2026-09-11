! INITTRAC_CART
subroutine inittrac_cart()
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
logical :: my_turn,old_sub=.false.

! PvK Feb 2 2021: initrac_cart_old is from before Feb 1, 2021. It does not stretch the tracers
! in the general case of (T)ALA even if stretch_tracers=T.
!old_sub=.true.
!if (old_sub) then
!   call inittrac_cart_old()
!   return
!endif

if (use_tracers_from_nate.or.use_tracers_from_cian) then
   if (mpi_partrac) then
      if (print_node) write(6,*) 'PERROR(inittrac_cart): not suited for use_tracers_from_XXX and mpi_partrac'
      call instop
   else
      call use_tracers_from_nate_or_cian()
   endif
   return  

else if (itracoption == 9) then
   if (mpi_partrac) then
      ! bit of a dumb way to do a parallel tracing test but it works
      call inittrac_cart_1(dxtrac,dytrac)
      write(lu_out,*) 'tracer1: ',tracer(1)%x,tracer(1)%y
      write(lu_out,*) 'tracerN: ',tracer(ntrac(1))%x,tracer(ntrac(1))%y,ntrac(1)
   else
      ! initial condition JGR97 App C
      ntot = 1
      ndist = 1
      ntrac(1) = 1
      tracer(1)%x = 0.5 
      tracer(1)%y = 0.0159221 
   endif
   return
else if (itracoption == 2) then
   ! markerchain method
   if (nochain.gt.1) then
      if (print_node) write(irefwr,*) 'PERROR(inittrac_cart): nochain>1'
      call instop
   endif

   ntot = 0
   ichain=1
   ip=0
   piw = pi/rlampix
   if (print_node) write(irefwr,*) 'rlampix, dm: ',rlampix,dm(1)
   imark(ichain) = int(rlampix/dm(ichain)+1)
   ntrac(ichain) = imark(ichain)
   dm(ichain) = rlampix/(ntrac(ichain)-1)
   ntot = ntot+ntrac(ichain)
   if (ntot.gt.NTRACMAX) then
      if (print_node) write(irefwr,*) 'PERROR(inittrac_cart): ntot>NTRACMAX'
      if (print_node) write(irefwr,*) 'ntrac, NTRACMAX: ',ntot,NTRACMAX
      call instop
   endif
   do i=1,ntrac(ichain)
      xm=dm(ichain)*(ntrac(ichain)-i)
      ym=y0m(ichain)+ainit(ichain)*cos(piw*xm)
      tracer(ip+i)%x = xm
      tracer(ip+i)%y = ym
   enddo
   if (print_node) write(irefwr,'(''PINFO(inittrac_cart): nmark = '',i10)') ntot
   if (print_node) write(irefwr,'(''PINFO(inittrac_cart): first marker: '',2f15.4)') tracer(1)%x,tracer(1)%y
   if (print_node) write(irefwr,'(''PINFO(inittrac_cart): last marker : '',2f15.4)') tracer(ntot)%x,tracer(ntot)%y

else if (itracoption == 1) then

   ! tracer method

   if (print_node) write(irefwr,*) 'mpi_partrac: ',mpi_partrac,numprocs
   if (ibench_type == 4) then
      call inittrac_cart_1_bench4(dxtrac,dytrac)

   else if (ibench_type==101) then
 
      if (print_node) write(irefwr,*) 'update for ibench_type==101'
      call instop
      !! "Stokes sphere" benchmark
      call inittrac_cart_blob()

   else ! not ibench_type=4,101

     call inittrac_cart_1(dxtrac,dytrac)

   endif ! ibench_type /= 4, 101

! make sure everyone knows what the global number of tracers (ntracg) is...
ntracg=ntrac
#ifdef MPI
   if (mpi_partrac.and.numprocs>1) then
      call MPI_REDUCE(ntrac,ntracg,NDISTMAX,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if (ierr/=0) then 
         write(irefwr,*) 'PERROR(inittrac_cart) myid : ',myid,' reports error ',ierr,' after reduce'
         call instop
      endif
      call MPI_BCAST(ntracg,NDISTMAX,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (ierr/=0) then 
         write(irefwr,*) 'PERROR(inittrac_cart) myid : ',myid,' reports error ',ierr,' after broadcast'
         call instop
      endif
      write(irefwr,'(''Processor '',i4,'' has '',i10,'' tracers or '',f8.2,''%'')') myid,ntrac(1),ntrac(1)*100.0_8/ntracg(1)
      if (print_node) write(irefwr,'(''    Total number of tracers is : '',2i10)') ntracg(1),ntracg(2)
      ! call instop
   endif
#endif

   ! Tracer distributions have been set up for all cases
   if (ibench_type<100) then         
      ! set densmark
      ip1=0
      if (its_CH94 .and. CH94_blob) then
           idist=1
           ! test with a thin layer. Compress particles in y direction to
           ! simulate effect of OC extraction
!          write(irefwr,*) 'ainittrac_CH94 = ',ainittrac_CH94,pi
           if (print_node) then
              open(99,file='tracer_init.dat')
              do i=1,ntrac(idist)
                 write(99,'(2f15.7,$)') tracer(ip1+i)%x,tracer(ip1+i)%y
                 dz=y0tr(idist)-tracer(ip1+i)%y
                 dz=-dz/8
                 x=tracer(ip1+i)%x
                 tracer(ip1+i)%y=y1tr(idist)-dz
                 write(99,'(f15.7,$)') tracer(ip1+i)%y
                 tracer(ip1+i)%y=tracer(ip1+i)%y-ainittrac_CH94*sin(pi*(x-2.5e0_8))
                 write(99,'(f15.7)') tracer(ip1+i)%y
              enddo
              close(99)
              write(irefwr,*) tracer(ntrac(idist))%x,tracer(ntrac(idist))%y
           endif
      endif
      if (ndist>1 .and. .not.ratio_method) then
         if (print_node) write(irefwr,*) 'PERROR(inittrac): code needs fixing for ndist>1 when not using the ratio method'
         call instop
      endif
      ip1=0
      do idist=1,ndist
         volume_dist(idist) = (x1tr(idist)-x0tr(idist))*(y1tr(idist)-y0tr(idist))
!        scale density in tracers with volume of tracers; this is global quantity 
         if (compress.and.stretch_tracers) then
            volume_per_tracer(idist)=dxtrac*dytrac
         else
            ! pre compressible definition
            volume_per_tracer(idist) = volume_dist(idist)/(ntracg(idist)) ! numprocs = 1 in serial
         endif
         if (T_buoyancy_through_particles) then
            ! debug case where temperature buoyancy is carried through the particles
            do  i=1,ntrac(idist)
                xm=tracer(ip1+i)%x
                ym=tracer(ip1+i)%y
                tracer(ip1+i)%density=Ra*(1-ym+0.1*sin(pi*ym)*cos(pi*xm))*volume_per_tracer(idist)
            enddo
         else
            tracer(ip1+1:ip1+ntrac(idist))%density = -Rb_local*volume_per_tracer(idist)
            if (its_CH94) tracer(ip1+1:ip1+ntrac(idist))%density=tracer(ip1+1:ip1+ntrac(idist))%density*C_peridotite
         endif
         if (print_node) then
            write(irefwr,*) 
            write(irefwr,*) 'PINFO(inittrac_cart): initialized tracers'
            write(irefwr,*) 'Distribution number : ',idist
            write(irefwr,*) '  ntrac = ',ntrac(idist)
            write(irefwr,*) '  density = ',tracer(ip1+1)%density, tracer(ip1+1)%density*ntrac(idist)
         endif
  
         if (mpi_partrac.and.numprocs>1) then
            write(tname,'(''tracer_info.'',i3.3,''_'',i3.3)') myid,idist
         else
            write(tname,'(''tracer_info.'',i3.3)') idist
         endif
         open(9,file=tname)
         write(9,*) 'ntrac  = ',ntrac(idist)
         write(9,*) 'xrange = ', minval(tracer(ip1+1:ip1+ntrac(idist))%x), & 
              & maxval(tracer(ip1+1:ip1+ntrac(idist))%x)
         write(9,*) 'yrange = ', minval(tracer(ip1+1:ip1+ntrac(idist))%y), & 
               & maxval(tracer(ip1+1:ip1+ntrac(idist))%y)
         write(9,*) 'density= ',tracer(ip1+1)%density,tracer(ip1+ntrac(idist))%density
         write(9,*) 'volume_per_tracer : ',volume_per_tracer(idist)
         close(9)

        ip1=ip1+ntrac(idist)
    enddo ! idist=1,ndist  

    ntot=sum(ntrac(1:ndist))

    if (its_CH94) tracer(1:ntot)%meltable=.true.
    
  endif ! ibench_type /= 101


else 

   if (print_node) then
      write(irefwr,*) 'PERROR(inittrac_cart): unknown option'
      write(irefwr,*) 'itracoption = ',itracoption
   endif
   call instop
 
endif  ! itracoption choices


end subroutine inittrac_cart

