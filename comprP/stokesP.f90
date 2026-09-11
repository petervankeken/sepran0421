subroutine stokesP(ichoice)
use sepmoduleoldrouts
use sepmodulecomio
use sepmodulecpack
use sepran_arrays
use sepran_interface
use brandenburg
use geometry
use control
use dtm_elem
use coeff
use tracers
use mpetrac
implicit none
integer,intent(in) :: ichoice
integer :: ipuser=0,ip=0,ipbuoy2=0,iread=0,ndim=2,npoint=0,iprint_solution=-1,ihelp=0
integer :: num_build=0,num_solve=0,build_in(20)=0,solve_in(20)=0,iout=0,ielhlp=0,isubiter=0
real(kind=4) :: t11,t12,timeA,timeB,cpu_here,tbuild0=0.0,tbuild1=0.0
real(kind=8) :: solve_rin(20),vrms,rmax,rhmax
character(len=80) :: cally,rname

!pedebug=.true.
!pedebug=.false.

call cpu_time(timeA)
! make sure b.c. are enforced if necessary
if (cyl) call velocity_bc(0)

if (abs(Rb_local)>1e-7) then
   ibuoy_trac=ichoice
else
   ibuoy_trac=0
endif
   
if (gable_stokes_choice==0) then
   no_buoyancy=.false.
else
   no_buoyancy=.true.
endif


!write(6,*) 'ibuoy_trac: ',ibuoy_trac,no_buoyancy,Rb_local

if (.not.no_buoyancy .and. abs(Rb_local) > 1e-3) then
    ! Prepare lookup tables in the case of chemical buoyancy
    ! NB: replaced in sepran 1020 by tracrhsd to prep info for sepelload2 instead of using el3003
    !if (itracoption==2) then
    !   !    Markerchain method: prepare pixel lookup grid
    !   write(irefwr,*) 'PERROR(stokesP): needs updating for marker chain'
    !   iout=0
    !   !if (ibuoy_trac.eq.1) then
    !      call mardiv(rname,iout)
    !   !else if (ibuoy_trac.eq.2) then
    !   !  call mardiv(coornewm,rname,iout)
    !   !endif
    !endif ! itracoption==2
!
    if (itracoption==1.and.fieldC) then
       ! compute compositional vector C and store in idens
       call tracdens()
       if (ivb) then
            if (print_node) write(irefwr,*) 'PERROR(stokesP) :: not suited yet for ivb'
            call instop
            ! set up pixel grid in case of compositionally dependent viscosity
            ! under development....
            if (ibuoy_trac==1) then
               !     call tractopix(kmesh1,kprob1,coormark,densmark)
            else
               !     call tractopix(kmesh1,kprob1,coornewm,densmark)
            endif
       endif
    else if (itracoption>=1) then
       ! Find barycentric coordinates and element numbers in tracers
       iout=0
       rname='markers.dat'
       !write(irefwr,*) 'mardiv'
       if (itracoption==2) call mardiv(rname,iout)
       !write(irefwr,*) 'mardiv done'
       if (ibuoy_trac > 0) then
           ! prep lookup information for more efficient rhsd assembly
           ! write(irefwr,*) '1 dtout: ',dtout
           !write(irefwr,*) 'call tracrhsd'
           call tracrhsd()
           !write(irefwr,*) 'done with tracrhsd'
           !call instop
       endif
   endif  

endif ! .not.no_buoyancy and |Rb_local|>0


!pedebug=.true.
if (gable_stokes_choice==0) then
   ! default case
   !call cpu_time(t11)
   ! done in bmstokes
   !call pefilcof(1,iuser_here,user_here)
   !call cpu_time(t12)
   if (pedebug.and.print_node) write(irefwr,*) 'pefilcof took : ',t12,t11
   if (pedebug.and.print_node) write(irefwr,*) 'irhsd1: ',irhsd1
   isubiter=0
   call bmstokes(kmesh1,kprob1,intmat1,iuser_here,user_here,matr1,matrm1, &
        & isol1,isolold1,idens,irhsd1)
   !write(6,'(''kmesh1 etc: '',6i10)') kmesh1,kprob1,intmat1,matr1,isol1,irhsd1

   if (compress.and..not.tala.and.ala_subiter) then
      do
         isubiter=isubiter+1
         dif1=anorm(0,3,0,kmesh1,kprob1,isol1,isolold1(1),ihelp)
         if (print_node) write(irefwr,'(''     ALA subiter: '',i5,e15.7)') isubiter,dif1
         if (dif1>=eps_ala_subiter.and.isubiter<N_ala_subiter_max) then
            call copyvc(isol1,isolold1)
            call bmstokes(kmesh1,kprob1,intmat1,iuser_here,user_here,matr1,matrm1, &
                 & isol1,isolold1,idens,irhsd1)
         else
            exit
         endif
       enddo
   endif
      
   !rhmax = anorm(1,3,0,kmesh1,kprob1,irhsd1,isol1,ielhlp)
   !write(6,*) 'rhmax = ',rhmax
   if (pedebug.and.print_node) write(irefwr,*) 'irhsd1: ',irhsd1
     !rmax = anorm(1,3,0,kmesh1,kprob1,isol1,isol1,ielhlp)
     !rmax = anorm(1,5,1,kmesh1,kprob1,isol1,isol1,ielhlp)
     !write(irefwr,*) 'l2 u vel for regular stokes: ',rmax
     !rmax = anorm(1,5,2,kmesh1,kprob1,isol1,isol1,ielhlp)
     !write(irefwr,*) 'l2 v vel for regular stokes: ',rmax
     !write(ourplotname,'(''PLOTS/VEL.'',i3.3)') 0
     !call plotvc(1,2,isol1,isol1,kmesh1,kprob1,0d0,1d0,0d0)
     !write(irefwr,*) 'gable_stokes_choice==0'
     !write(irefwr,*) 'max vel for regular stokes: ',rmax,rhmax
     !if (pedebug) call pevrms(vrms,isol1)
     !write(irefwr,*) 'vrms=',vrms
   call cpu_time(t11)
   if (subtract_rotation) then
      !write(irefwr,*) 'need to add subtract_average_rotation'
      !rmax = anorm(1,1,0,kmesh1,kprob1,isol1,isol1,ielhlp)
      !write(irefwr,*) 'max vel for regular stokes: ',rmax
      cally='stokesP gable_stokes_choice=0'
      call subtract_average_rotation(cally)
      !rmax = anorm(1,1,0,kmesh1,kprob1,isol1,isol1,ielhlp)
      !write(irefwr,*) 'max vel for regular stokes: ',rmax
   endif
     rmax = anorm(1,3,0,kmesh1,kprob1,isol1,isol1,ielhlp)
     !write(irefwr,*) 'max vel for regular stokes: ',rmax,rhmax
   call cpu_time(t12)
   if (pedebug.and.print_node) write(irefwr,*) 'subtract rotation took: ',t12-t11
   if (pedebug) then
      call pevrms(vrms,isol1)
      if (print_node) write(irefwr,*) 'vrms=',vrms
   endif
   !if (petest) call instop
else
   ! special case for force-balanced plates
   !call cpu_time(t11)
   !Done in bmstokes
   !call pefilcof(4,iuser_here,user_here)
   !call cpu_time(t12)
   if (pedebug.and.print_node) write(irefwr,*) 'pefilcof took : ',t12,t11
   if (pedebug.and.print_node) write(irefwr,*) 'irhsd1: ',irhsd1
   call bmstokes(kmesh1,kprob1,intmat1,iuser_here,user_here,matr1,matrm1,&
       & plate_sol(gable_stokes_choice),isolold1,idens,irhsd1)
     rmax = anorm(1,1,1,kmesh1,kprob1,isol1,isol1,ielhlp)
     !write(irefwr,*) 'max u vel for regular stokes: ',rmax
     rmax = anorm(1,1,2,kmesh1,kprob1,isol1,isol1,ielhlp)
     !write(irefwr,*) 'max v vel for regular stokes: ',rmax
     !write(ourplotname,'(''PLOTS/VEL.'',i3.3)') 0
     !call plotvc(1,2,isol1,isol1,kmesh1,kprob1,0d0,1d0,0d0)
   if (pedebug.and.print_node) write(irefwr,*) 'irhsd1: ',irhsd1
   if (pedebug) call pevrms(vrms,plate_sol(gable_stokes_choice))
endif
call cpu_time(timeB)
cpu_stokes=timeB-timeA
nstokes_solves=nstokes_solves+1
if (pedebug.and.print_node) write(irefwr,*) 'stokesP itself took ',timeB-timeA
if (print_node.and.pedebug) write(irefwr,*) 'gable_stokes_choice, vrms = ',gable_stokes_choice,vrms

end subroutine stokesP
