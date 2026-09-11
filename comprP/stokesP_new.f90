! Solve Stokes equation 
! Build matrix
! If necessary construct rhsd for plate solutions
! Solve 
subroutine stokesP_new(ichoice)
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
integer :: ipuser=0,ip=0,ipbuoy2=0,iread=0,ndim=2,npoint=0,iprint_solution=-1,iplate=0,jplate=0
integer :: num_build=0,num_solve=0,build_in(20)=0,solve_in(20)=0,iout=0,ielhlp=0,sinvec(10)=0
real(kind=4) :: t11,t12,timeA,timeB,cpu_here,gable_time
real(kind=8) :: solve_rin(20)=0,vrms,rmax,rhmax,rsinvec(10)=0,teststress,strchain_integrate
character(len=80) :: cally,rname
logical :: first=.true.,jp_orig_remove_rotation=.true.

!pedebug=.true.
!pedebug=.false.

!write(irefwr,*) 'in StokesP'
call cpu_time(timeA)
! make sure b.c. are enforced if necessary

if (cyl) then
   gable_stokes_choice=0
   call velocity_bc(0)
   if (gable_plates) then
     do iplate=1,nplates
        gable_stokes_choice=iplate
        call velocity_bc(0)
     enddo
   endif
endif

if (abs(Rb_local)>1e-7) then
   ibuoy_trac=ichoice
else
   ibuoy_trac=0
endif
   
! First set up matrix and vector with buoyancy
pedebug=.false.
gable_stokes_choice=0
no_buoyancy=.false.
! Test 042924. Code fails in sepprofile on a segfault on second iteration. Maybe decomposition is not stored correctly
! Build every step now for imatrix9=1 or imatrix9=2
 if (first.or.itypv>1.or.imatrix9<3) then
     if (pedebug.and.print_node) write(irefwr,*) 'Stokes: build A and f',gable_stokes_choice,itypv
     build_in=0
     build_in(1)=10 ! maximum number of entries in this array
     build_in(2)= 1 ! build matrix and rhsd vector
     build_in(10)= 3 ! number of old vectors
     first=.false.
 else 
     if (pedebug.and.print_node) write(irefwr,*) 'Stokes: build just f',gable_stokes_choice,itypv
     build_in=0
     build_in(1)=10
     build_in(2)=2 ! build f but not A
     build_in(10)= 3 ! number of old vectors
 endif

call pefilcof(1,iuser_here,user_here)

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
           ! write(irefwr,*) 'call tracrhsd'
           call tracrhsd()
           !write(irefwr,*) 'done with tracrhsd'
           !call instop
       endif
   endif  

endif ! .not.no_buoyancy and |Rb_local|>0

!pedebug=.true.
call sepbuildold(build_in,matr1,intmat1,kmesh1,kprob1,irhsd1,matrm1,isol1,isolold1,iuser_here,user_here)
call cpu_time(timeB)
if (pedebug) write(irefwr,*) 'Done building f (and S)',irhsd1,timeB-timeA

call cpu_time(timeA)
if (gable_plates) then
   call cpu_time(timeA)
   ! set up plate_rhsd(1:nplates) with essential b.c.
   do iplate=1,nplates
      call maver(matr1,plate_sol(iplate),plate_rhsd(iplate),intmat1,kprob1,6)
   enddo
   call cpu_time(timeB)
   !write(irefwr,*) 'nplate*maver took ',timeB-timeA
endif
call cpu_time(timeB)
gable_time=timeB-timeA

solve_in=0
solve_rin=0.0_8
solve_in(1)=3
solve_in(3)=0 ! direct solution method

if (ipetsc9==2) then
   solve_in(1)=10
   solve_in(3)=19 ! mumps
   if (itypv<2) then
      if (pedebug.and.print_node) write(irefwr,*) 'keep MUMPS matrix after factorization',itypv,gable_plates,gable_stokes_choice
      if (mumps_matrix_factorized) then
         if (pedebug.and.print_node) write(irefwr,*) 'reuse factorized MUMPS matrix'
         solve_in(10)=2  ! signal to MUMPS to keep matrix after factorization
      else
         solve_in(10)=1  ! signal to MUMPS to keep matrix after factorization
      endif
      mumps_matrix_kept=.true.
      mumps_matrix_factorized=.true.
   else if (itypv>=2.and.(gable_plates.and.gable_stokes_choice==0)) then
      ! First Stokes solve for gable_plates
      solve_in(10)=1  ! signal to MUMPS to keep matrix after factorization. 
      mumps_matrix_kept=.true.
      mumps_matrix_factorized=.true.
      if (pedebug.and.print_node) write(irefwr,*) 'keep MUMPS matrix after factorization',itypv,gable_plates,gable_stokes_choice
  endif
  solve_in(2)=0
  if (itype_stokes==900.and.mcontv==0) then
     solve_in(2)=1 ! positive definite (1) or not (0)
     ! write(irefwr,*) 'positive definite'
  endif
endif

iread=-1 ! indicate information is taken from solve_in, not stdin
if (pedebug.and.print_node) write(irefwr,*) 'irhsd1: ',irhsd1
if (pedebug.and.print_node) write(irefwr,'(''solve_in: '',20i5)') solve_in(1:14)
if (pedebug.and.print_node) write(irefwr,*) 'Stokes: solve Ax=f with: ',isolmethod9
call cpu_time(timeB)
if (pedebug.and.print_node) write(irefwr,*) 'prepare solvel: ',timeB-timeA

call cpu_time(timeA)
!write(6,'(''solvel with solve_in: '',20i10)') solve_in(1:20)
call solvel(solve_in,solve_rin,matr1,isol1,irhsd1,intmat1,kmesh1,kprob1,iread)
if (cyl) call subtract_average_rotation(kmesh1,kprob1,isol1,iuser_here,user_here)
if (debug) then 
   call pevrms(vrms,isol1)
   if (print_node) write(irefwr,*) 'initial vrms: ',vrms
endif
call cpu_time(timeB)
!pedebug=.true.
if (pedebug.and.print_node) write(irefwr,*) 'solvel took: ',timeB-timeA
pedebug=.false.


call cpu_time(timeA)
if (gable_plates) then

   ! compute stresses due to convective solution
   call strchain_compute(user_here,iuser_here,0,0)
   teststress = strchain_integrate(0)

   do iplate=1,nplates
      TP(iplate) = strchain_integrate(iplate)
   enddo
   !pedebug=.true.
   if (pedebug.and.print_node) then
     !write(irefwr,'(''Stress before '',10e15.7:)') TP(1:nplates)
     call cpu_time(timeB)
     write(irefwr,*) 'build TP: ',timeB-timeA
   endif
   pedebug=.false.

   ! find test solutions for each of the plates
   if (ipetsc9==2) solve_in(10)=2 ! signal to reuse factorized MUMPS matrix
   call cpu_time(timeA)
   do iplate=1,nplates
      call solvel(solve_in,solve_rin,matr1,plate_sol(iplate),plate_rhsd(iplate),intmat1,kmesh1,kprob1,iread)
      if (debug.and.print_node) write(irefwr,*) 'plate vrms: ',iplate,vrms
   enddo

   !write(irefwr,*) 'jp_orig_remove_rotation: ',jp_orig_remove_rotation
   if (jp_orig_remove_rotation) then
      ! Construct the net rotation vector from segments of the test solutions:
      if (nplates /= 8) then
          if (print_node) write(irefwr,*) 'PERROR(do_gable_solve): adjust jlimit for nplates other than 8'
          call instop
      endif
      call jlimit(user_here,iuser_here)
      ! Subtract the net rotation from the test solutions
      rsinvec(2) = 1.0d0
      sinvec(1) = 3
      sinvec(2) = 27
      sinvec(3) = 0
      do iplate=1,nplates
         rsinvec(1) = -1.0d0
         call manvec(sinvec,rsinvec,irotation,plate_sol(iplate),plate_sol(iplate),kmesh1,kprob1)
      enddo
   else
      do iplate=1,nplates
         call subtract_average_rotation(kmesh1,kprob1,plate_sol(iplate),iuser_here,user_here)
      enddo
   endif
   do iplate=1,nplates
      call pevrms(vrms,plate_sol(iplate))
      if (debug.and.print_node) write(irefwr,*) 'vrms ',iplate,'= ',vrms
   enddo

   call cpu_time(timeB)
   if (print_node.and.debug) write(irefwr,*) 'Plate solutions took : ',timeB-timeA
   ! compute stresses and form GP
   do iplate=1,nplates
      gable_stokes_choice = iplate
      call strchain_compute(user_here,iuser_here,iplate,0)
      do jplate=1,nplates
         GP(jplate,iplate) = -strchain_integrate(jplate)
      enddo
   enddo
   !debug=.true.
   if (debug.and.print_node) then
      write(irefwr,*) 'G: '
      do iplate=1,nplates
         write(irefwr,'(8e15.7)') (GP(jplate,iplate),jplate=1,nplates)
      enddo
   endif
   ! Now, solve the force balance equations
   ! gs is the Numerical Recipes Gaussian Solver
   if (pedebug.and.print_node) write(irefwr,*) 'Solving Gv = t with Gaussian Solver'
   if (periodic) then
      call gs(nplates-1,nplates,GP,TP,plate_vel)
   else
      call gs(nplates,nplates,GP,TP,plate_vel)
   endif
   ! Superimpose the Convective solution and the Corrected Plate test Solutions
   rsinvec(2) = 1.0d0
   sinvec(1) = 3
   sinvec(2) = 27
   sinvec(3) = 0

   if (periodic) then
      rsinvec(1) = plate_vel(1)
      call manvec(sinvec,rsinvec,plate_sol(1),isol1,isol1,kmesh1,kprob1)
      call manvec(sinvec,rsinvec,plate_sol(nplates),isol1,isol1,kmesh1,kprob1)
      do iplate=2,nplates-1
         rsinvec(1) = plate_vel(iplate)
         call manvec(sinvec,rsinvec,plate_sol(iplate),isol1,isol1,kmesh1,kprob1)
      enddo

   else
      do iplate=1,nplates
         !write(irefwr,*) 'super impose: ',i
         rsinvec(1) = plate_vel(iplate)
         call manvec(sinvec,rsinvec,plate_sol(iplate),isol1,isol1,kmesh1,kprob1)
      enddo
   endif

endif
call cpu_time(timeB)
if (debug.and.print_node) write(irefwr,*) 'total gable work took',gable_time+timeB-timeA

! Finally, remove null space from cylindrical models if necessary
!if (subtract_rotation) then
!   cally='stokesP'
!   call subtract_average_rotation(cally)
!endif
call cpu_time(t12)
if (pedebug.and.print_node) write(irefwr,*) 'subtract rotation took: ',t12-t11
!pedebug=.true.
if (pedebug) then
   call pevrms(vrms,isol1)
   if (print_node.and.print_node) write(irefwr,*) 'vrms=',vrms
   !call instop
endif
pedebug=.false.

call cpu_time(timeB)
if (pedebug.and.print_node) write(irefwr,*) 'stokesP itself took ',timeB-timeA
if (print_node.and.pedebug) write(irefwr,*) 'gable_stokes_choice, vrms = ',gable_stokes_choice,vrms

end subroutine stokesP_new
