subroutine test_tracer_accuracy()
use sepmodulecomio
use mtime
use mpetrac
use tracers
use coeff
use geometry
use sepran_arrays
use mparallel
use control
implicit none
!include 'SPcommon/ctimen'
real(kind=8) :: x2,y2,tstart_tracer,xp,yp
integer :: ntot2,inpcre(20),ntrac2,itracoption2,ndist2
integer :: ntot,ichoice
real(kind=8) :: format,yfaccn,factor,rincre(20),t,dtout2,tmax2,tstep,time_now2
real(kind=4) :: t01


if (print_node) then
   write(irefwr,*) '********************************************'
   write(irefwr,*) '** Test tracer accuracy '
   write(irefwr,*) '** 1: interpolation of analytical function '
endif
t00 = second()
ibuoy_trac=1
if (print_node) write(irefwr,*) 'call detelemtrac(1)'
ichoice=1
call detelemtrac(ichoice,'compr_start')
! Then, check accuracy using a quadratic test function (see func(10)).
if (print_node) write(irefwr,*) 'call detelemtrac(10)'
ichoice=10
call detelemtrac(ichoice,'compr_start')

!write(irefwr,*) 'stop after detelemtrac10'
!call instop



t01 = second()
if (print_node) then
   write(irefwr,*) '**    detelemtrac(2x) took ', t01-t00,' seconds'
endif


if (itracoption /= 0) then
  if (.not.cyl.and.rlampix.ge.1d0.or.(.not.axi.and..not.eighth)) then
     if (print_node) then
        write(irefwr,*) '**  '
        write(irefwr,*) '** 2: Van Keken et al., JGR, 1997, benchmark'
     endif
!    Keep copies of initial conditions 
     dtout2=dtoutp
     tmax2=tmaxp
     time_now2=time_now
     dtoutp = 0.2203712_8
     tmaxp = 0.2203712_8
     xp = 0.5_8
     yp = 0.015922_8
     x2 = tracer(1)%x
     y2 = tracer(1)%y
     ndist2 = ndist
     ntrac2 = ntrac(1)
     itracoption2=itracoption
     call copyvc(isol1,isolold1)

     if (.not.mpi_partrac) then
!       Actual test as intended
        itracoption=1
        if (cyl) then
           ! make consistent with func(4&5)
           if (r1<0.5) then
              tracer(1)%x=  0.5_8 + 0.4192_8 
           else
              tracer(1)%x = 0.5_8 + r1
           endif
           tracer(1)%y=0.0159221_8 - 0.5_8
           if (quart) then
!             rotate 45 degrees counterclockwise
              xp = (tracer(1)%x-tracer(1)%y)*sqrt(0.5_8)
              yp = (tracer(1)%x+tracer(1)%y)*sqrt(0.5_8)
              tracer(1)%x=xp
              tracer(2)%y=yp
           endif
        else ! .not.cyl
           tracer(1)%x=xp
           tracer(2)%y=yp
        endif
        ntot = 1
        ndist = 1
        ntrac(1) = 1
      else
        ntot=ntrac(1)
        ndist=1
        tracer(1)%x=xp
        tracer(1)%y=yp
      endif
      call detelemtrac(1,user_here,'compr_start')
 
!     Synthetic velocity; create velocity vector
!     Total number of entries
      inpcre(1)=10
!     Number of vectors to be created
      inpcre(2)=1
!     Type of vector  (1=solution vector)
      inpcre(3)=1
!     problem number
      inpcre(4)=1
!     sequence number of array of special structure
      inpcre(5)=0
!     0=real vector
      inpcre(6)=0
!     IFILL: filling is defined by IFILL commands
      inpcre(7)=3
!     fill first degree of freedom
      inpcre(8)=1
!     fill using FUNC with ichoice=4
      inpcre(9)=4
!     fill all nodes
      inpcre(10)=0
      call creatv(kmesh1,kprob1,isol1,inpcre,rincre)
!     fill second degree of freedom
      inpcre(8)=2
!     fill using FUNC with ichoice=5
      inpcre(9)=5
!     fill all nodes
      inpcre(10)=0
      call creatv(kmesh1,kprob1,isol1,inpcre,rincre)
      if (print_node) then
         write(ourplotname,'(''SYNVEL.000'')') 
         factor=0d0
         yfaccn=1d0
         format=10d0
         call plotvc(1,2,isol1,isol1,kmesh1,kprob1,format,yfaccn,factor)
      endif
!     if (print_node) then
         ! write(irefwr,*) 'call advance_tracers'
         if (ibench_type /= 4) then
            ! test is useless in box with aspect ratio 0.91
            call advance_tracers_JGR97()
            write(irefwr,*) '**********************************'
          endif
!     endif
!     For testing parallel efficiency of tracing itself
!tt   if (mpi_partrac) then
!tt      if (print_node) write(6,*) 'PINFO(test_tracer_accuracy): mpi_partrac and itracoption=9'
!tt      if (print_node) write(6,*) 'finishing up'
!tt      call instop
!tt   endif

      ! Now do Nate Sime's 2020 rotational test with time-dependent (u,v)=(-y+0.5,x+0.5)(cos^2(t)+0.5)
      ! which has a period of cycling of 2pi. 
      ! NB. Nate initially used a [-1,-1]x[1,1] square with rotation centered around (0,0)
      ! and initial particle position at (0.5,0)
      

!     Restore values
      tmaxp=tmax2
      dtoutp=dtout2
      time_now=time_now2
      ntot = ntot2
      tracer(1)%x = x2
      tracer(2)%y = y2
      ndist = ndist2
      ntrac(1) = ntrac2
      itracoption=itracoption2
      call copyvc(isolold1,isol1)
      if (tracer(1)%y*10000 == rlampix*10000) then
         if (print_node) write(irefwr,*) 'reset problem'
         call instop
      endif
   endif
endif

!if (petest) then
!   write(irefwr,*) 'stopping in test_tracer_accuray'
!   call instop
!endif

if (itracoption == 9) then
   ! We're done since the above code portion executes the tracer test
   if (print_node) then
     write(irefwr,*) 'PINFO(test_tracer_accuracy): stop after tracer test when itracoption=9'
   endif
   call instop()
endif


end subroutine test_tracer_accuracy
