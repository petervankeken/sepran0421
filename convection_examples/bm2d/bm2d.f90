! updated bm2d with two problem in kprob setup
! rewritten from scratch 
! PvK July 2024
program bm2d
use sepmodulekmesh
use sepran_arrays
use control
use coeff
use geometry
implicit none
!integer :: kmesh=0,kprob=0,isol(2)=0,intmat(2)=0,istream=0
integer :: inpcre(10)=0,ncntln=0,istream_input(5)=0
real(kind=8) :: rincre(10)=0.0_8,contln=0,format=10.0_8,DONE=1.0_8
!integer, parameter :: NUSERMAX=10000
!integer :: iuser_here(NUSERMAX)
!real(kind=8), dimension(:), allocatable :: user_here
integer :: i,num_user_here,allocate_status

print_node=.true.

call bmstart()

if (steady) then
   call steady_iteration
else
   if (print_node) write(irefwr,*) 'not yet suited for time dependence'
   call instop
endif


call finish()

end program bm2d

