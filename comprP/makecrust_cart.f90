! Makecrust
! Make crust for CH94
! zmelt=0.08, zcrust=zmelt*C_peridotite=0.01
subroutine makecrust_cart(ichoice)
use sepmodulecomio
use coeff
use tracers
use mpetrac
use mparallel
use geometry
use brandenburg
use control
implicit none
logical :: extract
integer :: itrac,ip1,ip2,ichoice,ntot,idum
integer :: imelt
real(kind=4) :: ran0
real(kind=8) :: basmin,basmax,harzmin,harzmax,ybefore

!     dummy call, when crust extracton is turned off
if (ichoice == 0) return

ntot=sum(ntrac(1:ndist))

imeltnow=0
! initial crust formation even in entire domain. 
! Only first idist=1 (basalt tracers) when fieldC or tracerC
! Melt harzburgite too (by moving any Hz particles in new crust down) with ratio ! method
if (ichoice == -1) then
   ! initialize meltable on first call
   do itrac=1,ntrac(1)
      if (tracer(itrac)%y>(1.0_8-zmelt)) then
         ybefore=tracer(itrac)%y
         imeltnow(1,1)=imeltnow(1,1)+1
         tracer(itrac)%y= 1.0_8-zcrust*ran0(idum)
         tracer(itrac)%meltable=.false.  ! indicates particle has recently been melted
      endif
   enddo
   if (print_node) then
      write(irefwr,'(''basalt melted initially: '',i10)') imeltnow(1,1)
      write(irefwr,*) 'ratio_method=',ratio_method,ratio_method_overlay,ntrac(1)+1,ntrac(2)
      write(irefwr,*) 'Hz min/max: ',minval(tracer(ntrac(1)+1:ntrac(2))%y),maxval(tracer(ntrac(1)+1:ntrac(2))%y)
   endif
   if (ratio_method.and.ratio_method_overlay) then
      ! second tracer distribution is originally the same as the first
      do itrac=ntrac(1)+1,ntot
         if (tracer(itrac)%y>(1.0_8-zcrust)) then
            ybefore=tracer(itrac)%y
            imeltnow(2,1)=imeltnow(2,1)+1
            zharz=zcrust+(zmelt-zcrust)*ran0(idum)
            tracer(itrac)%y= 1.0_8-zharz
            tracer(itrac)%meltable=.false.  ! indicates particle has recently been melted
         endif
      enddo
   endif

else if (ichoice == 1) then
  ! melt only in the melt zones at a given time step

  imeltnow=0
  ! First the basalt tracers. Move from (0,zmelt) to a random position in  (0,zcrust)
  do itrac=1,ntrac(1)
     extract=.false.
     if (tracer(itrac)%meltable.and.tracer(itrac)%y>1.0_8-zmelt) then
       ! particle may be melted
       do imelt=1,nmeltzone
         if (tracer(itrac)%x>=xmelt_left(imelt) .and. tracer(itrac)%x<=xmelt_right(imelt)) then
            imeltnow(1,imelt)=imeltnow(1,imelt)+1
            extract=.true.
         endif
       enddo
       if (extract) then
           ybefore=tracer(itrac)%y
           tracer(itrac)%y = 1.0_8-zcrust*ran0(idum)
           tracer(itrac)%meltable=.false.  ! indicates particle has recently been melted 
        endif
     endif ! meltable and particle in melt zone
     if (tracer(itrac)%y<1.0_8-zmelt) then
        tracer(itrac)%meltable = .true. ! indicates particle can be melted (again)
     endif

  enddo
  if (print_node) write(irefwr,*) 'tracers melted in B: ',imeltnow(1,1:nmeltzone)
  if (ratio_method.and.ratio_method_overlay) then
     ! Also move Hz tracers but now from (0,zcrust) to (zcrust,zmelt)
     do itrac=ntrac(1)+1,ntot
        extract=.false.
        if (tracer(itrac)%meltable.and.tracer(itrac)%y>1.0_8-zcrust) then
           ! figure out if it is in a meltzone
           do imelt=1,nmeltzone
              if (tracer(itrac)%x>=xmelt_left(imelt).and.tracer(itrac)%x<=xmelt_right(imelt)) then
                 imeltnow(2,imelt)=imeltnow(2,imelt)+1
                 extract=.true.
              endif
           enddo ! meltzones
           if (extract) then
              ! Hz particle may be 'melted'
              ybefore=tracer(itrac)%y
              imeltnow(2,imelt)=imeltnow(2,imelt)+1
              zharz=zcrust+(zmelt-zcrust)*ran0(idum)
              tracer(itrac)%y= 1.0_8-zharz
              tracer(itrac)%meltable=.false.  ! indicates particle has recently been melted
           endif ! extract
        endif  ! meltable tracers
     enddo ! loop over second batch of tracers
     if (tracer(itrac)%y<1.0_8-zmelt) then
        tracer(itrac)%meltable = .true. ! indicates particle can be melted (again)
     endif
  endif ! ratio_method and ratio_method_overlay
     
  !call instop
endif

end subroutine makecrust_cart

