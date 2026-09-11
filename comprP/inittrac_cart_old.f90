! INITTRAC_CART
subroutine inittrac_cart_old()
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
character(len=80) :: tname
integer :: nxtrac_here,nyhere
real(kind=8) :: density_ratio,factor_a,sumhi,scale_sumhi,htop,mass_now,massgoal,htot,funccf,stretchcorrect
integer :: nytractot
real(kind=8),dimension(:),allocatable :: htrac,ystretch,zstretch

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
   
else if (itracoption == 9) then
   ! initial condition JGR97 App C
   ntot = 1
   ndist = 1
   ntrac(1) = 1
   tracer(1)%x = 0.5 
   tracer(1)%y = 0.0159221 
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
     !if (ratio_method) then
     !   if (print_node) then
     !     write(irefwr,*) 'PERROR(inittrac_cart): ratio_method'
     !     write(irefwr,*) 'has not been set up yet for ibench_type=4'
     !   endif
     !   call instop
     !endif
!    Special case for Rayleigh-Taylor instability
     if (ainittr(1) <= 0d0) then
       if (print_node) write(irefwr,*) 'PERROR(inittrac_cart): ainittr(1) < 0: ', ainittr(1)
       if (print_node) write(irefwr,*) 'which is inconsistent with ibench_type=4'
       call instop
     endif
     ip1=0
     idist = 1
     nxtrac = int(rlampix/dr)
     dr = rlampix/nxtrac
     dxtrac = dr
     dytrac = dxtrac
     eps_tracer_bound(1:ndist) = 0 ! 0.5*dxtrac
     itrac=0
     if (mpi_partrac.and. (nxtrac/numprocs)*numprocs /= nxtrac) then
        if (print_node) then
          write(irefwr,*) 'PERROR(inittrac): tracers for RT benchmark '
          write(irefwr,*) 'need to be set for parallel computing such '
          write(irefwr,*) 'that nxtrac divides without remainder by '
          write(irefwr,*) 'the number of processors, but here '
          write(irefwr,*) '  nxtrac   = ',nxtrac
          write(irefwr,*) '  numprocs = ',numprocs
        endif
        call instop
     endif
     if (print_node) write(irefwr,*) 'inittrac_cart: nxtrac= ',nxtrac,rlampix,dr
     if (compress.and.stretch_tracers) then
       if (.not.fake_rho_bar.or.step_rho_background.or.exp_rho_background) then
          ! find distances between tracers that decrease as density increases
          nytractot=nint(2.0/dytrac)
          allocate(htrac(nytractot),zstretch(nytractot),ystretch(nytractot))
          if (print_node) write(irefwr,*) 'allocate for nytractot'
          htrac(1)=dytrac
          zstretch(1)=htrac(1)/2
          massgoal=dytrac ! assume rho=1 at top...
          ! make first estimate of layerthicknesses (should decrease with depth)
          htot=htrac(1)
          i=0
          do
             i=i+1
             if (i>nytractot) then
                if (print_node) write(irefwr,*) 'PERROR(inittrac-cart): estimate for nytractot is too small: ',nytractot
                call instop
             endif
             zstretch(i+1)=zstretch(i)+htrac(i)
             htrac(i+1)=htrac(i)
             mass_now=funccf(3,0.0_8,1.0_8-zstretch(i+1),0.0_8)*htrac(i+1)
             ! correct thickness of layer
             htrac(i+1)=htrac(i)*massgoal/mass_now
             ! reassign zstretch
             zstretch(i+1)=zstretch(i)+htrac(i+1)
             htot=htot+htrac(i+1)
             if (htot>1.0_8) exit
          enddo
          nyhere=i+1
          ! now restretch it a bit to make an integer number of tracers fit
          htrac(:)=htrac(:)/htot
          zstretch(1)=0.5*htrac(1)
          ystretch(nyhere)=1.0_8-zstretch(1)
          do i=2,nyhere
             zstretch(i)=zstretch(i-1)+htrac(i)
             ystretch(nyhere-i+1)=1.0_8-zstretch(i)
          enddo
          if (print_node) then
             open(9,file='ystretch.dat')
             write(9,*) 0,ystretch(1)
             do i=2,nyhere
                write(9,*) 0,ystretch(i),1.0_8/(ystretch(i)-ystretch(i-1))
             enddo
             close(9)
          endif
         !write(irefwr,*) 'stop in inittrac_cart'
         !call instop
          
       else
         density_ratio=1+drho_background_dense
         ! first make a regular set of tracer centers in 1D that stretch out along the density profile
         ! we're assuming particle distance for top most particle h_n and that of bottom most particle h_1 are
         ! at a ratio dictated by the density ratio: h_n = density_ratio*h1. Sum h_i = 1. 

         ! if density=1+z*drho_background_dense then the ratio between two adjacent
         ! h's is a: h_i=a*h_i+1. Then h_n=a^n*h_1 or a=density_ratio^1/n. Or a=10^(1/n*log10(density_ratio))
         factor_a=10**(dytrac*log10(density_ratio))
         nytractot=nint(1.0_8/dytrac)
         if (print_node) write(irefwr,*) 'nytractot: ',nytractot,factor_a,density_ratio
         allocate(htrac(nytractot),ystretch(nytractot),zstretch(nytractot))
         ! set trial h1 at top
         htrac(1)=1.0_8
         do i=2,nytractot
            htrac(i)=htrac(i-1)*factor_a
         enddo
         if (print_node) write(irefwr,*) 'htrac1,N: ',htrac(1),htrac(nytractot)
         if (print_node) write(irefwr,*) 'total length of htrac: ',sum(htrac(1:nytractot))
         scale_sumhi=1.0_8/sum(htrac(1:nytractot))
         htrac(1:nytractot)=htrac(1:nytractot)*scale_sumhi
         if (print_node) write(irefwr,*) 'total length of htrac: ',sum(htrac(1:nytractot))
         ! now form the array of center locations 
         ystretch(1)=0.5*htrac(1)
         htop=htrac(1)
         do i=2,nytractot
            ystretch(i)=htop+0.5*htrac(i-1)
            htop=htop+htrac(i)
         enddo
         if (print_node) then
            open(9,file='ystretch.dat')
            do i=1,nytractot
               write(9,*) 0.,ystretch(i)
            enddo
            nyhere=nytractot
            close(9)
         endif
       endif

       ! now fill the tracer array by picking all stretched particles in a column below the interfac
       do ix=myid+1,nxtrac,numprocs  ! stripe tracers in x direction
          xm = (ix-0.5)*dxtrac
          ymax = y1tr(idist)+ainittr(idist)*cos(pi*xm/rlampix)
          iy=0
          do 
            iy=iy+1
            if (ystretch(iy)<=ymax.and.iy<=nyhere) then
               ! particle is still in dense layer
               ym=ystretch(iy)
               itrac = itrac+1
               if (itrac.gt.NTRACMAX) then   
                  if (print_node) then
                    write(irefwr,*) 'PERROR(inittrac_cart): '
                    write(irefwr,*) '    itrac>NTRACMAX '
                    write(irefwr,*) '    itrac,NTRACMAX = ',itrac,NTRACMAX
                  endif
                  call instop
               endif
               tracer(ip1+itrac)%x = xm
               tracer(ip1+itrac)%y = ym
            else
               exit
            endif
          enddo ! loop over y positions
          
          if (ix==1.and.print_node) then 
            open(9,file='firstcolumn.dat')
            do i=1,itrac
               write(9,*) tracer(ip1+i)%x,tracer(ip1+i)%y
            enddo
            close(9)
          endif
       enddo ! loop over x positions
       ntrac(1)=itrac
     else ! not compress
       ! first the set of dense tracers
       do ix=myid+1,nxtrac,numprocs  ! stripe tracers in x direction
          xm = (ix-0.5)*dxtrac
          ymax = y1tr(idist)+ainittr(idist)*cos(pi*xm/rlampix)
          nytrac = int((ymax-0.5*dytrac)/dytrac + 1)
          do iy=nytrac,1,-1
             ym = (iy-0.5)*dytrac
             itrac = itrac+1
             if (itrac.gt.NTRACMAX) then   
                if (print_node) then
                  write(irefwr,*) 'PERROR(inittrac_cart): '
                  write(irefwr,*) '    itrac>NTRACMAX '
                  write(irefwr,*) '    itrac,NTRACMAX = ',itrac,NTRACMAX
                endif
                call instop
             endif
             tracer(ip1+itrac)%x = xm
             tracer(ip1+itrac)%y = ym
!            if (iy==nytrac) write(irefwr,*) xm,ym
          enddo
       enddo
       ntrac(1)=itrac
       if (ratio_method) then
         ! set volume of distribution
         x0tr(2)=x0tr(1)
         x1tr(2)=x1tr(1)
         y0tr(2)=y1tr(1) ! note second set of tracers sit on top of first set
         y1tr(2)=1.0_8
         if (print_node) write(irefwr,*) 'neutral tracers'
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
                  if (print_node) then
                    write(irefwr,*) 'PERROR(inittrac_cart): '
                    write(irefwr,*) '    itrac>NTRACMAX '
                    write(irefwr,*) '    itrac,NTRACMAX = ',itrac,NTRACMAX
                  endif
                  call instop
               endif
               tracer(ip1+itrac)%x = xm
               tracer(ip1+itrac)%y = ym
!              if (iy==nytrac) write(irefwr,*) xm,ym
            enddo
         enddo
         ntrac(2)=itrac
         ndist=2
       endif
     endif ! compress
     if (print_node) then
          write(irefwr,*) 'made RT initial condition'
          ! call instop
     endif

  else if (ibench_type==101) then
 
    if (print_node) write(irefwr,*) 'update for ibench_type==101'
    call instop
    !! "Stokes sphere" benchmark

    if (ratio_method) then

       ! call blobs_init_ratiomethod

    else ! not ratio method; just discretize spheres.
      !! first compute number of tracers
      itrac=0
      do iblob=1,nblob
         nypix=int((2*r_blob(iblob)-dr)/dr+1) ! 
         do iy=nypix/2,-nypix/2,-1 ! start at the top
            yb=y_blob(iblob)+iy*dr ! global coordinates
            theta=acos(iy*dr/r_blob(iblob))
            xbl=r_blob(iblob)*sin(theta)
            nxpix=int(2*xbl/dr) + 1 ! one in center
            itrac=itrac+2*(nxpix/2)+1
         enddo
      enddo
      if (print_node) write(irefwr,*) 'call allocate tracers: ',itrac
      ntrac_allocated=itrac
      call allocate_tracers(ntrac_allocated,1,'inittrac')

      icenterblob=0
      ip=0
      do iblob=1,nblob
         itrac=0
         ! one tracer is at the center of the blob and has local
         ! coordinates (centerpoint) (0,0).
         nypix=int((2*r_blob(iblob)-dr)/dr+1)  ! 
         do iy=nypix/2,-nypix/2,-1 ! start at the top
            yb=y_blob(iblob)+iy*dr ! global coordinates
            theta=acos(iy*dr/r_blob(iblob))
            xbl=r_blob(iblob)*sin(theta)
            nxpix=int(2*xbl/dr) + 1 ! one in center
!                 write(irefwr,*) 'blobs: ',yb,xbl,nxpix,theta,ip+itrac
            do ix=-nxpix/2,nxpix/2
               xb=x_blob(iblob)+ix*dr
               itrac=itrac+1
               if (ip+itrac>ntrac_allocated) then
                  if (print_node) write(irefwr,*) 'error: ip+itrac>ntrac_allocated'
                  if (print_node) write(irefwr,*) iblob,ip+itrac,ntrac_allocated
                  call instop
               endif
               tracer(ip+itrac)%x=xb
               tracer(ip+itrac)%y=yb
               if (ix==0.and.iy==0) icenterblob(iblob)=ip+itrac
            enddo
!                 if (nxpix==1) write(irefwr,*) 'xb,yb: ',xb,yb
         enddo
         ntrac(iblob)=itrac
         volume_per_tracer(iblob)=pi*r_blob(iblob)*r_blob(iblob)/ ntrac(iblob) 
         if (print_node) write(irefwr,*) 'fill tracer density',iblob,ip+1,ip+ntrac(iblob)
         do i=1,ntrac(iblob)
            tracer(ip+i)%density=-Rb_local*delta_rhop(iblob)*volume_per_tracer(iblob)
         enddo
         ip=ip+ntrac(iblob)
      enddo
      if (print_node) then
         open(9,file='blobs.dat')
         do itrac=1,sum(ntrac(1:nblob))
             write(9,'(3e15.7)') tracer(itrac)%x,tracer(itrac)%y, tracer(itrac)%density
         enddo
         open(9,file='centerblobs.dat') 
         do iblob=1,nblob
           write(9,'(3i10,3e15.7)') iblob,icenterblob(iblob), ntrac(1),tracer(icenterblob(iblob))%x-x_blob(iblob), &
             &      tracer(icenterblob(iblob))%y-y_blob(iblob),volume_per_tracer(iblob)
         enddo
         close(9)
         ndist=1
         ntrac(1)=sum(ntrac(1:nblob)) ! collect all blobs as active tracers into first distribution
         ntot=ntrac(1)
       endif
    endif ! not.ratio_method

   else ! not ibench_type=4,101

!    one dense layer without amplitude (for use in 2nd JGR97
!    benchmark or Tackley&King 2003 new benchmark
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
!    nxtrac = rlampix/dr
     nxtrac = int((x1tr(1)-x0tr(1))/dr)  ! rlampix/dr
     dr = (x1tr(1)-x0tr(1))/nxtrac ! rlampix/nxtrac
     dxtrac = dr
     dytrac = dxtrac
     eps_tracer_bound(1:ndist) = 0.5*dxtrac
     ip = 0
     ntot = 0
     do idist=1,ndist
        ymin = y0tr(idist)
        ymax = y1tr(idist)
        nytrac = int((ymax-ymin)/dytrac)
        if (mpi_partrac.and.numprocs>1) then
           nxtrac_here = nxtrac/numprocs
           dx = rlampix/numprocs
           xstart = x0tr(idist)+myid*dx
           xend   = x0tr(idist)+(myid+1)*dx
           if (nxtrac_here * numprocs /= nxtrac) then
              if (print_node) then
               write(irefwr,*) 'PERROR(inittrac): nxtrac_here*numprocs'
               write(irefwr,*) '    should be equal to nxtrac, but: '
               write(irefwr,*) '    nxtrac_here * numprocs = ', nxtrac_here*numprocs
               write(irefwr,*) '    nxtrac : ',nxtrac
              endif ! myid==0
              call instop
           endif
        else
           nxtrac_here = nxtrac
           dx = x1tr(idist)-x0tr(idist) ! rlampix
           xstart = x0tr(idist) ! 0
           xend = x1tr(idist) ! rlampix
        endif
        if (print_node) write(irefwr,'(''idist: '',3i10,4f12.3)') &
           & idist,nxtrac_here,nytrac,x0tr(idist),x1tr(idist),y0tr(idist),y1tr(idist) ! ,0d0,rlampix,ymin,ymax
        ntot = ntot + nxtrac*nytrac
        if (ntot>NTRACMAX) then
           if (print_node) then
              write(irefwr,*) 'PERROR(inittrac_cart): ntot>NTRACMAX'
              write(irefwr,*) 'ntot     = ',ntot
              write(irefwr,*) 'NTRACMAX = ',NTRACMAX
           endif
           call instop
        endif
        itrac=0
        do ix=1,nxtrac_here
           xhere = xstart+(ix-0.5)*dxtrac
           !write(irefwr,*) 'xhere: ',idist,xhere
           do iy=1,nytrac
              itrac=itrac+1
              tracer(itrac+ip)%x = xhere
              tracer(itrac+ip)%y = ymin+(iy-0.5)*dytrac
           enddo
        enddo
        ntrac(idist) = itrac
        ip = ip+itrac
        itrac = 0
      enddo ! idist=1,ndist
      ntot=0
      do idist=1,ndist
         ntot=ntot+ntrac(idist)
      enddo
     
   endif ! tracer distribution options

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
            volume_per_tracer(idist) = volume_dist(idist)/(ntrac(idist)*numprocs) ! numprocs = 1 in serial
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
  
            write(tname,'(''tracer_info.'',i3.3)') idist
            open(9,file=tname)
            write(9,*) 'ntrac  = ',ntrac(idist)
            write(9,*) 'xrange = ', minval(tracer(ip1+1:ip1+ntrac(idist))%x), & 
               & maxval(tracer(ip1+1:ip1+ntrac(idist))%x)
            write(9,*) 'yrange = ', minval(tracer(ip1+1:ip1+ntrac(idist))%y), & 
               & maxval(tracer(ip1+1:ip1+ntrac(idist))%y)
            write(9,*) 'density= ',tracer(ip1+1)%density,tracer(ip1+ntrac(idist))%density
            write(9,*) 'volume_per_tracer : ',volume_per_tracer(idist)
           close(9)
        endif ! myid==0

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

#ifdef MPI
!  write(irefwr,*) 'stopping at end of inittrac_cart_old'
!  call instop
#endif

end subroutine inittrac_cart_old
