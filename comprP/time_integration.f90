!   TIME_INTEGRATION
!
!   Do time integration of dynamical equations
!
!   Do for NOUT steps
!      tout = tout+dtout
!      While (t<tout)
!          compute CFL time step 
!          Integrate over new time step
!          output intermediate information
!      end
!      Output plots, write solution to file etc.
!   Done
!
!   PvK 071900
!   PvK 070119
!   Updated for 0619
subroutine time_integration
#ifdef MPI
use mpi
#endif
use sepmodulesol
use sepmodulecomio
use convparam
use coeff
use sepran_arrays
use sepran_interface
use geometry
use brandenburg
use control
use mtime
use tracers
use mpetrac
implicit none

! LOCAL VARIABLES
real t01,t11,t10,cpu,t12
real(kind=8) :: tnew,contln(10),format,yfaccn,t_d
integer :: iall,irhs2(5),jsmoot,numarr,ichoice
integer :: inbetween,ibp,imoved,npoint,ielhlp,ip,iout
real(kind=8) :: rmax,subdif,anorm,rmax2,p,q,veloc(20),veloc_d(20)
real(kind=8) :: vrms,vrms2,volint
logical :: output
character(len=80) :: fimage,command,rname
character(len=20) :: caller,cally 
integer :: iremarkstep,nq_a,ihelp,i1,inidgt,ipsol2,ierr,nusol,iinvol(10)
real(kind=8) :: q_a(20),q_a_d(20),tmin2,tmax2,outval(10),Ra_variable
real(kind=4) :: timeA,timeB
integer :: ih,im,is
real(kind=8) :: resmem,get_resmem


save irhs2,iremarkstep

!pedebug=.true.

niter = 0
iremarkstep=0
cpu=0
if (restart) then
   inout=int(time_now/dtoutp)+1
   toutp=inout*dtoutp
else
   inout = 1
   toutp = time_now + dtoutp
endif
if (pedebug.and.print_node) write(irefwr,*) 'inout = ',inout,'toutp=',toutp,'time_now=',time_now
inbetween = 1 ! not 0, duh
imoved = 0

t00 = second()
output = .false.
call pefilxy(2)

if (print_node) then
  write(irefwr,*)
  write(irefwr,*)
  write(irefwr,'(''Start time integration at time '',2e15.7)') time_now,time_now/tscale_dim
  write(irefwr,'(''                Next output at '',2e15.7,i5)') toutp,toutp/tscale_dim,inout
endif
call copyvc(isol1,isolold1)
if (.not. stokes_is_updated) then
   gable_stokes_choice=0
   call cpu_time(timeA)
   call stokesP(1)
   call cpu_time(timeB)
   if (print_node) write(irefwr,*) 'initial stokes solution took: ',timeB-timeA
   if (gable_plates) then
      call cpu_time(timeA)
      call do_gable_solve(iuser_here,user_here,niter,1)
      call cpu_time(timeB)
      if (print_node) write(irefwr,*) 'gable stokes solution took: ',timeB-timeA
   endif
endif
if (print_node) then
   call nusseltP(q_a)
   write(irefwr,'(''       initial solution q_a        : '', 3e12.2)') q_a(1)/q_a(3),q_a(2)/q_a(3),q_a(3)
endif

! Compute RMS velocity (PvK 070119: integr fails on kproba with sepstr)
call pevrms(vrms,isol1)
if (print_node) write(irefwr,*) 'vrms = ',vrms
if (petest) call instop
if (itracoption==1.and.track_cmb_ingas) call logcmb(0)
call time_int_output()


call cpu_time(cpu_after_start)
do
   niter = niter+1
   ! check on resident memory
   !resmem=get_resmem()
   !if (resmem>0.5088) red_flag=.true.
   
!    calculate time step; pefilcof stores velocity in user_here
!    which then is used in pedtcf.
     if (red_flag.and.print_node) write(irefwr,*) 'pefilcof: ',get_resmem()
     call pefilcof(2,iuser_here,user_here)
     if (red_flag.and.print_node) write(irefwr,*) 'pedtcf  : ',get_resmem()
     call pedtcf(kmesh1,iuser_here,user_here,dtcfl)
     tstepp = dtcfl * tfac
     ! see if time step should be limit in early stage (useful for RT instabilities)
     ! t_init_limit_max,dt_init_limit
     if (time_now<t_init_limit_max.and.tstepp>dt_init_limit) then
        tstepp=dt_init_limit
     endif
     tnew = time_now+tstepp
     ! now see if we need to limit time step because we're close to an output step
     if (tnew >= toutp) then
        tstepp = toutp-time_now
        time_now = toutp
        output = .true.
     else if (tnew >= tmaxp) then
         tstepp = tmaxp - tstepp
         time_now = tmaxp
         output = .true.
     else 
         time_now = tnew
         output = .false.
     endif
     if (use_varRa) then
        Ra=min(Rat_max,Ra_variable(time_now))
        !write(6,*) 'time_now, Ra: ',time_now,Ra
     endif

     !petest=.true.
     if (petest.and.print_node) write(irefwr,*) 'dtcfl, tstep, tnew: ',dtcfl,tstepp,tnew,time_now,output,toutp
     !petest=.false.

!    predic  position of markers
     call cpu_time(timeB)
     if (itracoption /= 0) then
        if (tfac>0.6d0.and.predict_with_RK4) then
           !write(irefwr,*) 'predict with RK4'
           call move_tracers4(2,isol1,isol1,tstepp,tfac)
           call cpu_time(timeA)
           if (print_node.and.(pedebug.or.petest)) write(irefwr,*) 'move_tracers took ',timeA-timeB
        else
           ! revert to just Euler predictor
           call cpu_time(timeB)
           call tracvel(kmesh1,kprob1,isol1,user_here)
           !call cpu_time(timeA)
           !if (pedebug.or.petest) write(irefwr,*) 'tracvel took ',timeA-timeB
           call predcoort()
           !call cpu_time(timeA)
           !if (pedebug.or.petest) write(irefwr,*) 'predcoort took ',timeA-timeB
           !call cpu_time(timeB)
           !call tracvel(kmesh1,kprob1,isol1,user_here,1)
           !if (pedebug.or.petest) write(irefwr,*) 'tracvel+predcoort took ',timeA-timeB
        endif
     endif
     !write(irefwr,*) 'done with predicting tracers'
     call cpu_time(timeA)
     cpu_tracers=cpu_tracers+timeA-timeB


     if (red_flag.and.print_node) write(irefwr,*) 'ieulh: ',get_resmem(),matrm2
     if (do_temperature_solution) then
        call cpu_time(timeB)
        call ieulh(kmesh1,kprob1,intmat2,isol2,isolold2,user_here,iuser_here,matr2,matrm2,irhsd2)
        call cpu_time(timeA)
        cpu_heat=cpu_heat+timeA-timeB
     endif

     if (ncor == 1) then
        call copyvc(isol1,isolold1)
        gable_stokes_choice=0
        if (red_flag.and.print_node) write(irefwr,*) 'pevrms: ',get_resmem()
        call pevrms(vrms,isol1)
        if (print_node.and.petest) write(irefwr,*) 'before stokes 2 w/ gable_stokes_choice=0: ',vrms
        if (print_node.and.petest) write(irefwr,*) 'stokesP(2)'
        call cpu_time(timeB)
        if (red_flag.and.print_node) write(irefwr,*) 'stokesP: ',get_resmem()
        call stokesP(2)
        if (red_flag.and.print_node) write(irefwr,*) 'pevrms: ',get_resmem()
        call pevrms(vrms,isol1)
        if (print_node.and.petest) write(irefwr,*) 'after stokes2 w/ gable_stokes_choice=0: ',vrms
        if (isnan(vrms)) call instop
        
        if (gable_plates) call do_gable_solve(iuser_here,user_here,niter,1)
        call cpu_time(timeA)
        cpu_stokes=cpu_stokes+timeA-timeB
        if (do_temperature_solution) then
           call cpu_time(timeB)
           if (red_flag.and.print_node) write(irefwr,*) 'corrh: ',get_resmem(),matrm2
           call corrh(kmesh1,kprob1,intmat2,isol2,isolold2,user_here,iuser_here,matr2,matrm2,irhsd2)
           call cpu_time(timeA)
           cpu_heat=cpu_heat+timeA-timeB
        endif
        if (itracoption>=1) then
           call cpu_time(timeB)
           call move_tracers4(2,isolold1(1),isol1,tstepp,tfac)
           !write(irefwr,*) 'done with correcting tracers'
           call cpu_time(timeA)
           cpu_tracers=cpu_tracers+timeA-timeB
        endif
     endif

     if (itracoption>=1) then
        ! tracer%x,y:= tracer%xnew,ynew
        if (itracoption == 1 .or. itracoption == 3) call copcoor
        if (itracoption == 2) call remarker
     endif

     ! update Stokes
     if (print_node.and.petest) write(irefwr,*) 'stokesP(1)'
     call cpu_time(timeB)
     gable_stokes_choice=0
     if (red_flag.and.print_node) write(irefwr,*) 'stokesP: ',get_resmem()
     call stokesP(1)
     if (petest) then
       call pevrms(vrms,isol1)
       if (print_node) write(irefwr,*) 'before stokes1 w/ gable_stokes_choice=0: ',vrms
     endif
     if (gable_plates) call do_gable_solve(iuser_here,user_here,niter,1)
     call cpu_time(timeA)
     cpu_stokes=cpu_stokes+timeA-timeB

     call cpu_time(timeB)
     if (itracoption==1.and.(its_CH94.or.gable_plates)) then
        call makecrust(1)
     endif
     if (itracoption==1.and.track_cmb_ingas) call logcmb(1)
     call cpu_time(timeA)
     cpu_tracers=cpu_tracers+timeA-timeB
     
     ! intermediate output
     if (red_flag.and.print_node) write(irefwr,*) 'time_int_output: ',get_resmem()
     call time_int_output()
     if (petest) call instop


     if (output) then
        ! big output step
        if (itracoption == 1 .or. itracoption == 2) then
            if (print_node) call tracerout(inout)
        endif
        if (red_flag.and.print_node) write(irefwr,*) 'peplotit: ',get_resmem()
        if (print_node) call peplotit
        toutp=time_now+dtoutp
        inout=inout+1
        output=.false.
        if (print_node) then
           call writbs_netcdf('T_restart.nf',isol2)
           call writbs_netcdf('UV_restart.nf',isol1)
        endif

        !call MPI_Barrier(MPI_COMM_WORLD,ierr)
        !if (print_node) write(irefwr,*) 'stopping at line 105'
        !call instop

     endif
     !write(irefwr,*) 'stop after time_int_output'
     !call instop

     resmem=get_resmem()
     if (resmem>resmem_max) then
        if (print_node) then
           write(irefwr,*) 'Exiting time loop because resident memory exceeds set limit: ',resmem,resmem_max
        endif
        exit
     endif

     if (time_now >= tmaxp) then
        if (print_node) then
           write(irefwr,*) 'Exiting time loop because t >= tmax: ',time_now>=tmaxp,time_now,tmaxp
        endif
        exit
     endif
enddo
!write(6,'(''number of vectors created: '',i5)') nextfreesol-1

!call bmout()

end subroutine time_integration

real(kind=8) function Ra_variable(t)
use coeff
implicit none
real(kind=8) :: t
if (Di<1e-7) then
   Ra_variable=1e4*exp(100*t)
else
   Ra_variable=1e4*exp(50*t)
endif
end function Ra_variable
