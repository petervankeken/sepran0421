subroutine bmstokes(kmesh1,kprob1,intmat1,iuser_here,user_here,matr1,matrm1,isol1,isolold1,idens,irhsd1)
use sepmoduleoldrouts
use geometry
use control
use brandenburg
use coeff
use sepmodulecomio
implicit none
integer :: kmesh1,kprob1,intmat1,iuser_here(*),matr1,isol1,isolold1,idens,irhsd1,matrm1
real(kind=8) :: user_here(*)
integer :: ipuser,ip,ipbuoy2,iread,ndim,npoint,iprint_solution=-1
integer :: num_build=0,num_solve=0,build_in(20),solve_in(20)
real(kind=4) :: t11,t12,timeA,timeB,cpu_time_here,tbuild00,tbuild0,tbuild1
real(kind=8) :: solve_rin(20),vrms,get_resmem
logical :: first=.true.



call cpu_time(t11)
build_in=0
!if (print_node) pedebug=.true.
cpu_time_here=0
call cpu_time(timeA)
if (print_node) pedebug=.true.
!pedebug=.true.
pedebug=.false.
if (gable_stokes_choice==0) then
  ! no Gable plates. Check to see if matrix needs to be rebuild
  ! if (.not.mumps_matrix_factorized.or.itypv>1) then
  if (first.or.itypv>1) then
     !if (pedebug) write(irefwr,*) 'Stokes: build A and f',gable_stokes_choice,itypv
     build_in=0
     build_in(1)=10 ! maximum number of entries in this array
     build_in(2)= 1 ! build matrix and rhsd vector
     build_in(10)= 3 ! number of old vectors
     first=.false.
   else
     !if (pedebug) write(irefwr,*) 'Stokes: build just f',gable_stokes_choice,itypv
     build_in=0
     build_in(1)=10
     build_in(2)=2 ! build f but not A
     build_in(10)= 3 ! number of old vectors
   endif
   call pefilcof(1,iuser_here,user_here)
   !if (print_node) pedebug=.true.
   if (print_node.and.red_flag) write(irefwr,*) '   build: ',get_resmem()
   call build(build_in,matr1,intmat1,kmesh1,kprob1,irhsd1,matrm1,isol1,isolold1,iuser_here,user_here)
   call cpu_time(timeB)
   if (pedebug) write(irefwr,*) 'Done building f (and S)',irhsd1,timeB-timeA
else
   ! gable_stokes_choice>1: rebuild only f
   !if (pedebug) write(irefwr,*) 'Stokes: build just f',gable_stokes_choice
   build_in=0
   build_in(1)=10
   build_in(2)=2
   build_in(10)= 3 ! number of old vectors
   call pefilcof(4,iuser_here,user_here)
   call cpu_time(tbuild0)
   call build(build_in,matr1,intmat1,kmesh1,kprob1,irhsd1,matrm1,isol1,isolold1,iuser_here,user_here)
   call cpu_time(tbuild1)
   tbuild00=tbuild00+tbuild1-tbuild0
   if (pedebug) write(irefwr,*) 'irhsd1: ',irhsd1
   if (pedebug) write(irefwr,*) 'done building just f for gable: ',irhsd1
endif
if (gable_stokes_choice==8.and.debug) then
   write(irefwr,*) 'building rhsd for plate solution took: ',tbuild00
   tbuild00=0
endif
call cpu_time(timeB)
!if (print_node) pedebug=.true.
if (pedebug) write(irefwr,*) 'build time = ',timeB-timeA
cpu_time_here=timeB-timeA+cpu_time_here

call cpu_time(timeA)
num_build=num_build+1
if (pedebug) write(irefwr,*) 'number of builds: ',num_build,matr1
if (pedebug) write(irefwr,*) 'irhsd1: ',irhsd1
!if (print_node) pedebug=.false.
if (pedebug) call sepgetrhsdinfo(irhsd1)

if (printmatrix) then
   write(irefwr,*) 'original sepran matrix'
   call prinmt(intmat1,matr1,kprob1)
endif

!call prinrv(irhsd1,kmesh1,kprob1,4,1,'irhsd')

!if (print_node) pedebug=.true.
! Solve system of equations (PG8.3)
solve_in=0
solve_rin=0.0_8
if (isolmethod9 < 0) then
  if (pedebug) write(irefwr,*) 'iterative method: ',itype_stokes
  solve_in(1)=14 ! maximum number of entries
  if (ipetsc9==1) then
     solve_in(3)= 20 ! 1=bicgstab; 20=Petsc solver
  else
     solve_in(3) = 1
  endif
  solve_in(4)= ipreco9 ! ILU preconditioner
  solve_in(5)= maxiter9! maxiter
  solve_in(6)= iprint_solution ! print level in sepran 
  solve_in(9)= 1 ! start with given vector
  if (itypv<2) solve_in(10)=1 ! signal sepran to keep decomposition of matrix
  solve_in(11)=ireler9 ! IRELER 1=take relative error compared to initial residual (10: default in sepran)
  solve_in(12)=0 ! compute and print 2 norm of residual
  solve_rin(1)=cgeps9 ! accuracy - but this is ignored with PETSC
  solve_rin(11)=ksp_abs9 ! absolute accuracy with KSP solvers
  solve_rin(12)=ksp_rel9    ! relative accuracy in KSP
  solve_rin(13)=1e6      ! div_tol in KSP
else 
  if (pedebug) write(irefwr,*) 'direct method: ',itype_stokes
  solve_in=0
  solve_in(1)=10
  solve_in(3)=0 ! direct solution method
  if (itypv<2) solve_in(10)=1 ! signal sepran to keep matrix decomposition
  if (ipetsc9==2) then
     solve_in(3)=19 ! mumps
     if (itypv<2) then
       if (pedebug) write(irefwr,*) 'keep MUMPS matrix after factorization',itypv,gable_plates,gable_stokes_choice
       if (mumps_matrix_factorized) then
          if (pedebug) write(irefwr,*) 'reuse factorized MUMPS matrix'
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
       if (pedebug) write(irefwr,*) 'keep MUMPS matrix after factorization',itypv,gable_plates,gable_stokes_choice
     else  if ((itypv<2 .or. (gable_plates.and. gable_stokes_choice>0)) .and. mumps_matrix_factorized) then
       if (pedebug) write(irefwr,*) 'reuse MUMPS matrix after factorization',irhsd1
       solve_in(10)=2  ! reuse matrix when viscosity structure is the same
     else
       if (pedebug) then
          write(irefwr,*) 'do not keep or reuse MUMPS matrix'
          write(irefwr,*) 'itypv=',itypv
       endif
       solve_in(10)=0  ! don't keep matrix because viscosity is changing
       mumps_matrix_kept=.false.
       mumps_matrix_factorized=.false.
     endif
     if (pedebug) write(irefwr,*) 'MUMPS: ',itypv,gable_plates,gable_stokes_choice,solve_in(10), &
        & mumps_matrix_kept,mumps_matrix_factorized
  endif
  solve_in(2)=0
  if (itype_stokes==900.and.mcontv==0) then
     solve_in(2)=1 ! positive definite (1) or not (0)
     ! write(irefwr,*) 'positive definite'
  endif
endif
iread=-1 ! indicate information is taken from solve_in, not stdin
if (pedebug) write(irefwr,*) 'irhsd1: ',irhsd1
if (pedebug) write(irefwr,'(''solve_in: '',20i5)') solve_in(1:14)
if (pedebug) write(irefwr,*) 'Stokes: solve Ax=f with: ',isolmethod9
call cpu_time(timeB)
if (pedebug) write(irefwr,*) 'prepare solvel: ',timeB-timeA
cpu_time_here=timeB-timeA+cpu_time_here

call cpu_time(timeA)
   if (print_node.and.red_flag) write(irefwr,*) '   solvel: ',get_resmem()
call solvel(solve_in,solve_rin,matr1,isol1,irhsd1,intmat1,kmesh1,kprob1,iread)
call cpu_time(timeB)
if (pedebug) write(irefwr,*) 'solvel took: ',timeB-timeA
cpu_time_here=timeB-timeA+cpu_time_here
num_solve=num_solve+1
if (pedebug) write(irefwr,*) 'number of solve: ',num_solve
call cpu_time(t12)
if (petiming .and. print_node) write(irefwr,*) 'bmstokes took: ',t12-t11
if (pedebug) write(irefwr,*) 'bmstokes took: ',t12-t11,cpu_time_here

!call prinrv(isol1,kmesh1,kprob1,1,0,'sol1')
!call prinrv(irhsd1,kmesh1,kprob1,1,0,'rhsd1')
!if (gable_stokes_choice==1) call instop


if (pedebug) write(irefwr,*) 'mumps_matrix: ',mumps_matrix_kept,mumps_matrix_factorized

!call prinrv(isol1,kmesh1,kprob1,4,1,'isol1')

!call pevrms(vrms,plate_sol(gable_stokes_choice))
!if (print_node) write(irefwr,*) 'gable_stokes_choice, vrms = ',gable_stokes_choice,vrms


end subroutine bmstokes
