! Calculate stresses and solve the force balance equations
! Originally - uses Numerical Recipes Gaussian Solver
! Feb. 2007 - uses LAPACK LU decomposition solver 
subroutine do_gable_solve(iuser,user,niterG,ichoice_stokes)
use sepmodulecomio
use convparam
use mtime
use sepran_arrays
use geometry
use brandenburg
use coeff
use control
implicit none
integer :: iuser(*)
real(kind=8) :: user(*)
integer,intent(in) :: niterG,ichoice_stokes
real :: timeA,timeB
real(kind=8) :: rsinvec(3),offdiag,aftertest
real(kind=8) :: teststress,strchain_integrate,vrms
real(kind=8) :: tau_start(10),tau_end(10),tau_drop(10)
integer :: sinvec(3),iplate,jplate

real(kind=8) :: x_c,omega
! Declaration for LAPACK
integer :: INFO,M,N,NRHS,RANK
real(kind=8) :: RCOND
integer,parameter :: LDA=100,LDB=100
integer :: IWORK(100)
real(kind=8) :: A2(LDA,10),B2(LDB,1),S(LDA),vp_sum
! LWORK is the dimension of WORK
integer, parameter :: LWORK=100000
real(kind=8) :: WORK(LWORK)
logical :: jp_orig_remove_rotation=.true.

! Use machine-precision for RCOND
parameter(RCOND = -1.0d0)

! GP - pass through a common block
!real(kind=8) :: GP(8,8),TP(8)
!common /cgbound/ GP

integer :: i,j

timeA=second()
!pedebug=.true.

! Note:  gable_stokes_choice > 0 will set the body force to zero for the Stokes problem, and
! activate the respective partial velocity boundary condition.  The respective
! solution will be stored in plate_sol(i)

! Stresses due to the previous convective solution
call strchain_compute(user,iuser,0,0)
teststress = strchain_integrate(0)

do i=1,nplates
   TP(i) = strchain_integrate(i)
enddo


!if (periodic) then
!   ! Combine plate#1 and plate#n
!   write(irefwr,*)'combining plates 1&n'
!   TP(1) = TP(1) + TP(nplates)
!!   write(irefwr,'(''TP: '',14x,10e15.7:)') (TP(i),i=1,nplates-1)
!else
!!   write(irefwr,'(''TP: '',14x,10e15.7:)') (TP(i),i=1,nplates)
!endif
!pedebug=.true.
if (pedebug.and.print_node) then
  write(irefwr,'(''Stress before '',10e15.7:)') TP(1:nplates)
  call cpu_time(timeB)
  write(irefwr,*) 'build TP: ',timeB-timeA
endif
pedebug=.false.


! compute the kinematic plate test solutions:
do i=1,nplates
   gable_stokes_choice = i
   call stokesP(ichoice_stokes)
   !write(irefwr,*) 'plate_sol: ',i,plate_sol(i)
enddo
timeB=second()
!pedebug=.true.
if (pedebug.and.print_node) write(irefwr,*) 'find test solutions: ',timeB-timeA
!     write(irefwr,*) 'time in stokes: ',timeB-timeA
!do i=1,nplates
!   write(irefwr,'(''GP : '',i5,8e12.3)') i,(GP(j,i),j=1,nplates)
!enddo


! write(irefwr,*)'Creating net-rotation vector'

!jp_orig_remove_rotation=.false.
!write(irefwr,*) 'jp_orig_remove_rotation: ',jp_orig_remove_rotation
if (jp_orig_remove_rotation) then
   ! Construct the net rotation vector from segments of the test solutions:
   if (nplates /= 8) then
       if (print_node) write(irefwr,*) 'PERROR(do_gable_solve): adjust jlimit for nplates other than 8'
       call instop
   endif
   call jlimit(user,iuser)
   ! Subtract the net rotation from the test solutions
   rsinvec(2) = 1.0d0
   sinvec(1) = 3
   sinvec(2) = 27
   sinvec(3) = 0
   do i=1,nplates
      rsinvec(1) = -1.0d0
      call manvec(sinvec,rsinvec,irotation,plate_sol(i),plate_sol(i),kmesh1,kprob1)
   enddo
else
   do i=1,nplates
      call subtract_average_rotation(kmesh1,kprob1,plate_sol(i),iuser_here,user_here)       
   enddo
endif
if (pedebug.and.print_node) then
    call cpu_time(timeB)
    if (print_node) write(irefwr,*) 'remove rotation: ',timeB-timeA
endif
if (debug) then
   do i=1,nplates
      call pevrms(vrms,plate_sol(i))
      if (debug.and.print_node) write(irefwr,*) 'vrms ',i,'= ',vrms
   enddo
endif
 

! Now compute the stresses and form GP
do i=1,nplates
   gable_stokes_choice = i
   call strchain_compute(user,iuser,i,0)
   do j=1,nplates
      GP(j,i) = -strchain_integrate(j)
   enddo
enddo
if (pedebug) then
   call cpu_time(timeB)
   if (print_node) write(irefwr,*) 'form GP: ',timeB-timeA
endif
   
!debug=.true.
if (debug.and.print_node) then
   write(irefwr,*) 'G: '
   do iplate=1,nplates
      write(irefwr,'(8e15.7)') (GP(jplate,iplate),jplate=1,nplates)
   enddo
endif

! Now, solve the force balance equations
! gs is the Numerical Recipes Gaussian Solver
if (pedebug.and.print_node) write(irefwr,*)'Solving Gv = t with Gaussian Solver'
if (periodic) then
   call gs(nplates-1,nplates,GP,TP,plate_vel)
else
   call gs(nplates,nplates,GP,TP,plate_vel)
endif
call cpu_time(timeB)
if (pedebug.and.print_node) write(irefwr,*) 'solve GP v = TP',timeB-timeA

vp_sum=0
do i=1,nplates
   vp_sum=vp_sum+plate_vel(i)
enddo
vp_sum=vp_sum*velocity_scale*100*year_dim
if (print_node.and.gable_output_choice>0) write(irefwr,'(''Plate Velocity = '',10f12.7:)')  &
   & (plate_vel(i)*velocity_scale*100*year_dim,i=1,nplates),vp_sum
!endif

! Superimpose the Convective solution and the Corrected Plate test Solutions
rsinvec(2) = 1.0d0
sinvec(1) = 3
sinvec(2) = 27
sinvec(3) = 0

if (periodic) then
   rsinvec(1) = plate_vel(1)
   call manvec(sinvec,rsinvec,plate_sol(1),isol1,isol1,kmesh1,kprob1)
   call manvec(sinvec,rsinvec,plate_sol(nplates),isol1,isol1,kmesh1,kprob1)
   do i=2,nplates-1
      rsinvec(1) = plate_vel(i)
      call manvec(sinvec,rsinvec,plate_sol(i),isol1,isol1,kmesh1,kprob1)
   enddo
 
else
   do i=1,nplates
      !write(irefwr,*) 'super impose: ',i
      rsinvec(1) = plate_vel(i)
      call manvec(sinvec,rsinvec,plate_sol(i),isol1,isol1,kmesh1,kprob1)
   enddo
endif

!Does this make the rotation worse? 
!if (subtract_rotation) then
!   call subtract_average_rotation(kmesh1,kprob1,isol1,iuser,user)
!endif
!pedebug=.true.
if (pedebug.and.print_node) then
   call pevrms(vrms,isol1)
   write(irefwr,*) 'final vrms: ',vrms
endif

! Check and see if the force balance had its desired effect:
! this quantity is not as important as the 
call strchain_compute(user,iuser,0,1)
teststress = strchain_integrate(0)
if (pedebug.and.print_node) write(irefwr,*)"stress after superposition: ",teststress
if (teststress>1e-8) then
   if (print_node) write(irefwr,*) 'PWARN(do_gable_stokes): teststress is large: ',teststress
endif

!       aftertest = strchain_integrate(1)
!       if (periodic) then
!           aftertest = aftertest + strchain_integrate(nplates)
!       endif

!      do i=1,nplates
!          aftertest = strchain_integrate(i)
!          write(irefwr,*)"Plate ",i," stress after sup. = ",aftertest
!          tau_end(i) = aftertest
!      enddo
!write(irefwr,'(''Stress after '',10f15.7)')  (strchain_integrate(i),i=1,nplates)

do i=1,nplates       
   tau_drop(i) = tau_end(i)
enddo
! Keep track of the stress drop after superposition: it won't be zero 
! if DGELSD is used with an overdetermined system
!write(73,'(10(e14.6,2x))') time_now,tau_drop(1),tau_drop(2),tau_drop(3),tau_drop(4),tau_drop(5),tau_drop(6),tau_drop(7),tau_drop(8)
!call flush(73)
pedebug=.false.


end subroutine do_gable_solve

! gs .... solve the system Ax = b ; A = GP and is passed thru cjp.inc    ********
subroutine gs(idimi,N,GP,b,x)
implicit none


real(kind=8) :: x(N),b(N),GP(N,N)
integer :: i,j,k,n,m,idimi,rowops,ridx,cidx
integer :: fac
real(kind=8) :: rfac

! GP - pass through a common block - its just easier this way.
!real(kind=8) :: GP(8,8)
!common /cgbound/ GP

!write(irefwr,*) 'GP x = plate_vel: '
!do i=1,8
!   write(irefwr,'(8e14.6,2x,e14.6)') GP(i,1:8),b(i)
!enddo

! Initialize Solution Vector
do i=1,idimi
   x(i) = 0.0d0
enddo

! Perform Gaussian Elimination to put Ax=b in RE Form
! rowops = fac(idimi - 1)
ridx = 2
cidx = 1
! Down the rows, starting with row #2
do k=1,(idimi-1)
   do i=ridx,idimi
      rfac = -GP(i,cidx)/GP(ridx-1,cidx)
      b(i) = b(i) + rfac*b(ridx-1)
      ! Across the columns
      do j=1,idimi
         GP(i,j) = GP(i,j) + rfac*GP(ridx-1,j)
      enddo
   enddo
   ridx = ridx + 1
   cidx = cidx + 1
enddo

ridx = idimi - 1
cidx = idimi
! Up the rows, starting with idimi-1
do k=1,(idimi-1)
   do i=ridx,1,-1
      rfac = -GP(i,cidx)/GP(ridx+1,cidx)
      b(i) = b(i) + rfac*b(ridx+1)
      ! Across the columns
      do j=1,idimi
         GP(i,j) = GP(i,j) + rfac*GP(ridx+1,j)
      enddo
   enddo
   ridx = ridx - 1
   cidx = cidx - 1
enddo

do i=1,idimi
   x(i) = b(i)/GP(i,i)
enddo

! Display Matrices
! write(irefwr,*)'A = '
!      do i=1,idimi
!         write(irefwr,'(15(f15.7))GP(i,1),GP(i,2),GP(i,3),GP(i,4),GP(i,5),
!     v               GP(i,6),GP(i,7),GP(i,8),GP(i,9),GP(i,10)
!      enddo
!
!      write(irefwr,*)'b = '
!      do i=1,idimi
!         write(irefwr,*) b(i)
!      enddo
!
!      write(irefwr,*)'x = '
!      do i=1,idimi
!         write(irefwr,*) x(i)
!      enddo
!

end subroutine gs

! Calculate conditions at the ridge
subroutine jlimit(user,iuser)
use sepran_arrays
use brandenburg
implicit none

real(kind=8) :: user(*)
integer :: iuser(*)
integer :: iinvec(10),ipl,plate_soll(5,8)
real(kind=8) :: rinvec(10)

integer icom
logical vcop
common /fv/ icom,vcop


iinvec(1) = 5
iinvec(2) = 32
! iinvec(3) = DOF for manipulation
iinvec(3) = 1
iinvec(4) = 1
! iinvec(5) = Number of vectors
iinvec(5) = 8
rinvec(1)=4.0_8 ! signal option=4 to funvec

!      icom = 1
!      call manvec(iinvec,rinvec,ign(1,igchoice),ign(1,igchoice)v,ign(1,igchoice),kmesh1,kprob1)
!      iinvec(3) = 2
!      icom = 2
!      call manvec(iinvec,rinvec,ign(1,igchoice),ign(1,igchoice)
!     v,ign(1,igchoice),kmesh1,kprob1)

! manvec still expects a 5xiinvec(5) array
plate_soll(1,1:8) = plate_sol(1:8)

! plate 1
icom = 1
ipl = 1
!write(irefwr,*) 'jlimit1'
call manvec(iinvec,rinvec,plate_sol,plate_sol,irotation,kmesh1,kprob1)

iinvec(3) = 2
icom = 2
!write(irefwr,*) 'jlimit2'
call manvec(iinvec,rinvec,plate_sol,plate_sol,irotation,kmesh1,kprob1)
end subroutine jlimit

