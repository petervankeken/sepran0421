subroutine steady_heatP()
use sepmoduleoldrouts
use sepmodulecomio
use sepmodulecpack
use sepran_arrays
use sepran_interface
use control
use dtm_elem
use geometry
implicit none
integer :: ipuser,ip,ipbuoy2,iread,ndim,npoint,ipvelx,ipvely,i
integer :: idegfd,ichbvl,jstap,ichoice,iprint_solution=-1
integer :: build_in(20),solve_in(20)
real(kind=8) :: solve_rin(20)

call sepgetmeshinfo(ndim,npoint)

pedebug=.false.
call pefilcof(2,iuser_here,user_here)
!pedebug=.true.
if (pedebug.and.print_node) then
   do i=1,150
      write(irefwr,'(2i15,e15.7)') i,iuser_here(i),user_here(i)
   enddo
endif
pedebug=.false.

! moved to setup_temperature_bc
!idegfd=1
!jstap=0
!ichoice=0
!ichbvl=1002 ! problem 2
!call bvalue(ichoice,ichbvl,kmesh1,kprob1,isol2,T_bot,ibottom,ibottom,idegfd,jstap)
!idegfd=1
!jstap=0
!ichoice=0
!ichbvl=1002 ! problem 2
!call bvalue(ichoice,ichbvl,kmesh1,kprob1,isol2,T_top,itop,itop,idegfd,jstap)
!if (pedebug) write(irefwr,*) 'done with bvalue'

build_in=0
build_in(1)= 2  ! maximum number of entries in this array
build_in(2)= 1 ! build matrix and rhsd vector
! not used if build_in(1)==2
build_in(13)=1 ! iseqin (sequence number of old solution in islold)
build_in(15)=subdivide800 ! isubdivide
!build_in(9)= 2 !iprob
call build(build_in,matr2,intmat2,kmesh1,kprob1,irhsd2,matrm2,isol2,isol2,iuser_here,user_here)

if (pedebug.and.print_node) write(irefwr,*) 'done with build'

if (printmatrix .and. print_node) then
   write(irefwr,*) 'original sepran matrix'
   call prinmt(intmat2,matr2,kprob1)
endif

!if (ipetsc8==2) then
!   write(irefwr,*) 'PERROR(steady_heatP): needs to be tested for solution with MUMPS'
!   call instop
!endif

! Solve system of equations (PG8.3)
solve_in=0
solve_rin=0.0_8
if (isolmethod8 < 0) then
   solve_in(1)=14 ! maximum number of entries
   if (ipetsc8==1) then
      solve_in(3)= 20 ! 1=bicgstab; 20=Petsc solver
   else
      solve_in(3) = 1
   endif
 solve_in(4)=ipreco8 ! ILU preconditioner
 solve_in(5)=maxiter8! maxiter
 iprint_solution=2
 solve_in(6)= iprint8 ! print level in sepran 
 solve_in(9)= 1 ! start with given vector
 solve_in(11)=ireler8 ! IRELER 1=take relative error compared to initial residual (10: default in sepran)
 solve_in(12)=0 ! compute and print 2 norm of residual
 solve_rin(1)=cgeps8 ! accuracy - but this is ignored with PETSC
 solve_rin(11)=ksp_abs8 ! absolute accuracy with KSP solvers
 solve_rin(12)=ksp_rel8    ! relative accuracy in KSP
 solve_rin(13)=1e6      ! div_tol in KSP

else if (isolmethod8 ==0 ) then
 if (ipetsc8==0) then
    ! standard situation: sepran profile method
    solve_in(1)=3
    solve_in(3)=0 ! direct solution method
 else if (ipetsc8==2) then
    solve_in(1)=10
    solve_in(3)=19 ! mumps
    solve_in(10)=0 ! do not keep factorized mtrix
    mumps_matrix_kept8=.false.
    mumps_matrix_factorized=.false.
 else 
    if (print_node) write(irefwr,*) 'PERROR(steady_heatP): unknown option for ipetsc8: ',ipetsc8
    call instop
 endif
endif
iread=-1 ! indicate information is taken from solve_in, not stdin

if (pedebug.and.print_node) write(irefwr,'(''solve_in: '',14i6)') solve_in(1:14)
call solvel(solve_in,solve_rin,matr2,isol2,irhsd2,intmat2,kmesh1,kprob1,iread)

end subroutine steady_heatP
