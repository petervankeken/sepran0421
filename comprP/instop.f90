subroutine instop
#ifdef MPI
use mpi
#endif
use sepmodulecomio
use mparallel
use control
implicit none
character(len=120) :: command
integer :: ierr
call finish(0)
#ifdef MPI
if (mpi_parallel) then
   ! Sepran has finished; ignore irefwr
   close(lu_out)
   close(7)
   close(14)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   irefwr=6
!  write(irefwr,*) 'finishing up process ',myid
   call MPI_FINALIZE(ierr)
   write(command,'(''rm -f input.'',i3.3)') myid
   call system(command)
else 
!  write(irefwr,*) 'finishing up serial run'
endif
#endif
if (print_node) call system("rm -rf fort??????")
stop
end
