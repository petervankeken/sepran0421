real function second()
#ifdef MPI
use mpi
#endif
implicit none
real time

!#ifdef MPI
!second = MPI_WTIME()
!return
!#endif

call cpu_time(time)
second = time

end function second
