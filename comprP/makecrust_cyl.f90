! Makecrust_cyl
! Make crust for gable_plates
subroutine makecrust_cyl(ichoice)
use sepmodulecomio
use coeff
use tracers
use mpetrac
use mparallel
use geometry
use brandenburg
use control
implicit none
integer :: ichoice
logical :: extract
integer :: itrac,ip1,ip2,ntot,idum,iplate
integer :: imelt,izone,lu,iingas,ntracer_ingassed
real(kind=4) :: ran0
real(kind=8) :: basmin,basmax,harzmin,harzmax,ybefore
real(kind=8) :: r,theta,cost,rbefore,x,y,dtheta

!     dummy call, when crust extracton is turned off
if (ichoice == 0) return

ntot=sum(ntrac(1:ndist))

if (ntot>int(1e7)) then
   if (print_node) then
     write(irefwr,*) 'PERROR(makecrust_cyl) update format in makecrust_cyl for write to log files'
     write(irefwr,*) 'for ntot > 1e7: ',ntot
   endif
   call instop
endif

if (mpi_parallel) then
   if (print_node) then
      write(irefwr,*) 'logging of particles needs to be extended to parallel'
   endif
   call instop
endif

!open(109,file='Ba_melted_now')
!open(110,file='Hz_melted_now')

imeltnow=0
! initial crust formation even in entire domain. 
! Only first idist=1 (basalt tracers) when fieldC or tracerC
! Melt harzburgite too (by moving any Hz particles in new crust down) with ratio method
if (print_node) write(irefwr,*) 'makecrust: ichoice, ntrac: ',ichoice,ntrac(1)
if (ichoice == -1) then
   ! initialize meltable on first call
   do itrac=1,ntrac(1) ! loop over basalt tracers
      x=tracer(itrac)%x
      y=tracer(itrac)%y
      r=sqrt(x*x+y*y)
      if (r>r2-zmelt) then
         if (print_node) then
            lu=201
            write(lu,'(i6,i7)') niter,itrac  ! ,0 do not add basalt indicator....
         endif
         cost=y/r
         theta=acos(cost)
         if (x<0.0_8) theta=2*pi-theta
         rbefore=r
         imeltnow(1,1)=imeltnow(1,1)+1
         r=r2-zcrust*ran0(idum)
         tracer(itrac)%x=r*sin(theta)
         tracer(itrac)%y=r*cos(theta)
         tracer(itrac)%meltable=.false.  ! indicates particle has recently been melted
      endif
   enddo ! end loop over basalt tracers
   if (print_node) write(irefwr,'(''basalt melted initially: '',i10)') imeltnow(1,1)
   if (ratio_method.and.ratio_method_overlay .or. duplicate_tracers) then
      ! second tracer distribution is originally the same as the first
      do itrac=ntrac(1)+1,ntot ! loop over Hz tracers
         x=tracer(itrac)%x
         y=tracer(itrac)%y
         r=sqrt(x*x+y*y)
         if (r>=r2-zmelt) then
            tracer(itrac)%meltable=.false.
            ! compact Hz tracers in entire melt layer, not just in the crust
            !if (r>=r2-zcrust) then ! Hz tracer is in the crust; move it down into the residue
            cost=y/r
            theta=acos(cost)
            if (x<0.0_8) theta=2*pi-theta
            rbefore=r
            if (print_node) then
               lu=201
               write(lu,'(i6,i7)') niter,itrac
            endif
            imeltnow(2,1)=imeltnow(2,1)+1
            if (r>=r2-zcrust.and.move_Hz_from_crust) then
               zharz=zcrust+(zmelt-zcrust)*ran0(idum)
               r=r2-zharz
               tracer(itrac)%x=r*sin(theta)
               tracer(itrac)%y=r*cos(theta)
            endif
         endif ! tracer in melt zone
      enddo ! end loop over Hz tracers
      if (print_node) write(irefwr,'(''Harzburgite melted initially: '',i10)') imeltnow(2,1)
   endif
   return
endif


if (ichoice == 1) then ! melt only in the melt zones at a given time step

  ! this sets up an array meltzone(1:nmeltzones) indicating where plates are
  ! diverging 
  call find_melt_zones
  imeltnow=0
  ntracer_ingassed=0
  ! First the basalt tracers. Move from (0,zmelt) to a random position in  (0,zcrust)
  do itrac=1,ntrac(1)
     extract=.false.
     x=tracer(itrac)%x
     y=tracer(itrac)%y
     r=sqrt(x*x+y*y)
     if (r<r2-zmelt) then
        tracer(itrac)%meltable = .true. ! indicates particle can be melted (again)
        tracer(itrac)%ingasable = .true. 
     endif
     if (r>=r2-zmelt .and. (tracer(itrac)%meltable.or.tracer(itrac)%ingasable)) then
         ! particle may be melted or ingassed. Check angular position
         cost=y/r
         theta=acos(cost)
         if (x<0.0_8) theta=2*pi-theta
         if (tracer(itrac)%meltable) then
           do imelt=1,nmeltzone
             izone=meltzone(imelt)
             dtheta=abs(theta-plate_boundaries(izone))
             ! capture particles on the left side of the polar melt region
             if (izone==1 .and. x<0.0_8) then
                dtheta=abs(2*pi-theta-plate_boundaries(izone))
             endif
             if (dtheta <= xmelt) then 
                ! log this tracer's melting
                lu=LU_MELT_START+izone
                if (print_node) write(lu,'(i6,i9)') niter,itrac 
                rbefore=r
                r=r2-zcrust*ran0(idum)
                tracer(itrac)%x= r*sin(theta)
                tracer(itrac)%y= r*cos(theta)
                tracer(itrac)%meltable=.false.  ! indicates particle has recently been melted
                imeltnow(1,izone)=imeltnow(1,izone)+1
             endif
           enddo ! loop over melt zones
         endif ! meltable
         if (tracer(itrac)%ingasable) then
           do iingas=1,ningaszone
             izone=ingaszone(iingas)
             dtheta=abs(theta-plate_boundaries(izone))
             ! capture particles on the left side of the polar melt region
             if (izone==1 .and. x<0.0_8) then
                dtheta=abs(2*pi-theta-plate_boundaries(izone))
             endif
             if (dtheta <= xmelt) then 
                ! log this tracer's ingassing
                lu=LU_INGAS_START+izone
                if (print_node) write(lu,'(i6,i9)') niter,itrac
                tracer(itrac)%ingasable=.false.
                ntracer_ingassed=ntracer_ingassed+1
             endif
           enddo
         endif ! ingasable
     endif ! meltable/ingasable and particle at melt zone depth
  enddo ! loop over basalt tracers
  if (print_node.and.gable_output_choice>0) then
    write(irefwr,'(''degassing          : '',$)') 
    do iplate=1,nplates
       if (diverging(iplate)) write(irefwr,'(i10,$)') iplate
    enddo
    write(irefwr,*) 
    write(irefwr,'(''tracers melted in B: '',$)') 
    do iplate=1,nplates
       if (diverging(iplate)) write(irefwr,'(i10,$)') imeltnow(1,iplate)
    enddo
    write(irefwr,*) 
  endif
  if (ratio_method.and.ratio_method_overlay .or. duplicate_tracers) then
     ! Also move Hz tracers but now from (0,zcrust) to (zcrust,zmelt)
     do itrac=ntrac(1)+1,ntot
        x=tracer(itrac)%x
        y=tracer(itrac)%y
        r=sqrt(x*x+y*y)
        if (r<r2-zmelt) then
           tracer(itrac)%meltable = .true. ! indicates particle can be melted (again)
           tracer(itrac)%ingasable = .true.
        endif
        extract=.false.
        ! ingassing only occurs for basalt tracers
        if ((tracer(itrac)%meltable).and.r>=r2-zmelt) then
            ! meltable Hz particle could be in melting region
            cost=y/r
            theta=acos(cost)
            if (x<0.0_8) theta=2*pi-theta
            do imelt=1,nmeltzone
              !write(irefwr,'(''   Hz: '',i10,2i5,3f12.3)') itrac,imelt,izone,r,theta,dtheta
              izone=meltzone(imelt)
              dtheta=abs(theta-plate_boundaries(izone))
              ! capture particles on the left side of the polar melt region
              if (izone==1 .and. x<0.0_8) dtheta=abs(2*pi-theta-plate_boundaries(izone))
              if (dtheta < xmelt) then 
                 ! log this tracer's melting or ingassing
                 lu=LU_MELT_START+izone
                 write(lu,'(i6,i9,i2)') niter,itrac
                 ! make sure to tag this one as melted
                 tracer(itrac)%meltable=.false.
                 ! and in the case of it being in the crust move it down
                 if (r>r2-zcrust.and.move_Hz_from_crust) then
                    zharz=zcrust+(zmelt-zcrust)*ran0(idum)
                    r=r2-zharz
                    tracer(itrac)%x= r*sin(theta)
                    tracer(itrac)%y= r*cos(theta)
                 endif
                 imeltnow(2,izone)=imeltnow(2,izone)+1
              endif
            enddo ! loop over melt zones
        endif ! meltable and particle at melt zone depth
     enddo ! loop over Hz tracers
     if (print_node.and.gable_output_choice>0) then
        write(irefwr,'(''tracers melted in H: '',$)') 
        do iplate=1,nplates
           if (diverging(iplate)) write(irefwr,'(i10,$)') imeltnow(2,iplate)
        enddo
        write(irefwr,*) 
     endif
  endif ! ratio_method and ratio_method_overlay .or. duplicate_tracers
  if (print_node) write(irefwr,'(''number of tracers ingassed: '',i10)') ntracer_ingassed
     
  !call instop
endif

do iplate=1,nplates
   call flush(LU_MELT_START+iplate)
   call flush(LU_INGAS_START+iplate)
enddo

!close(109)
!close(110)

end subroutine makecrust_cyl
