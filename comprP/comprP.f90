! Solve convection equations with petsc
! This approach follows closely what we use in the production codes (e.g., comprS)

! PG=Programmer's Guide (ta.twi.tudelft.nl/NW/sepran/pg.pdf)
! SP=Standard Problems Guide (ta.twi.tudelft.nl/NW/sepran/sp.pdf)

! PvK August 2017
! PvK October 2017: expanded for parallel computing
! PvK June 2019: complete rewrite using sepran 0619
program comprP
!use sepmodulemain
use sepmodulekmesh
use sepmodulekprob
use sepmodulecomio  ! provides lu's for reading & writing
use sepmodulecpack  ! parallel computing interface
use sepran_arrays ! access to kmesh, kprob, isol, etc. and interfaces to subroutines
use sepran_interface
use control
use brandenburg
use dtm_elem  ! controls building of pressure mass matrix
#ifdef MPI
use mpi
#endif
implicit none
integer :: iuserh(100),allocate_status,iinvol(10),isoldum=1
real(kind=8),allocatable :: userh(:)
real(kind=8) :: outval(10)

write(6,*) 'comprP_start'
call comprP_start()
!pedebug=.true.
!call sepgetmeshinfo(ndim,npoint)

!iinvol=0
!iuserh=0
!iuserh(1)=100
!iuserh(2)=1
!iuserh(6)=7
!iuserh(10)=2001
!iuserh(11)=6
!allocate(userh(5+npoint),stat=allocate_status)
!if (allocate_status /= 0) then
!   write(irefwr,*) 'allocate userh failed: ',allocate_status
!   call instop
!endif
!userh(6:5+npoint)=1.0_8
!userh(1)=1.0_8*(5+npoint)
!!write(irefwr,*) 'integr from comprP'
!call integr(iinvol,outval,kmesh,kprob,isoldum,iuserh,userh)
!write(irefwr,*) outval(1)

if (steady) then
   write(6,*) 'steady'
   call steady_iteration
else 
   call time_integration
   if (print_node) then
      close(LU_VRMS) 
      close(LU_NU) 
      close(LU_NU_BOT) 
      if (gable_plates) close(LU_ROTATION)
   endif
endif

call pvk_deallocate
!call finish(0)
call instop

end program comprP 

subroutine pvk_deallocate
use sepran_arrays
use sepran_interface
use coeff
implicit none

if (allocated(funcx)) deallocate(funcx)
if (allocated(funcy)) deallocate(funcy)
if (allocated(user_here)) deallocate(user_here)


end subroutine pvk_deallocate
