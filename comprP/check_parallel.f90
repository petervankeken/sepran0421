!   CHECK_PARALLEL
!   See whether the code runs in parallel with MPI or not.
!   Currently only used for parallel tracing and is completely 
!   independent of Sepran's domain decomposition for solving the PDEs in parallel
!
!   When using with mpi_parallel=T and mpi_partrac=F the code will just be
!   executed with multiple independent instances (useful for testing, of course).
!
!   020121 PvK Updated for sepran 1020
subroutine check_parallel
use sepmodulecomio
use sepmodulecpack ! provides firstmumps; we need to override the default to avoid errors with MPI_Finalize
use control
use mparallel
#ifdef MPI
use mpi
#endif
#ifdef IFORT
!use IFPORT
#endif
implicit none
logical :: check
integer :: ierr,lu_out_sepran
character(len=180) :: fname,command,hostname

firstmumps=.true.
irefwr=6

#ifdef MPI
! Check on existence of files indicating which type of parallelism is going to be used
! parallel tracing?
! PvK 0421: move to read_namelist
!inquire(file='mpi_partrac',exist=mpi_partrac)
! parallel execution outside of sepran?
!inquire(file='mpi_parallel',exist=mpi_parallel)

! parallel execution guided by sepran?
inquire ( file = 'sepran_par.input', exist = sep_parallel )

if (sep_parallel.and.(mpi_partrac.or.mpi_parallel)) then
   write(irefwr,*) 'PERROR(check_parallel): sepran_par.input exists suggesting sepran parallelism'
   write(irefwr,*) 'This is incompatible with mpi_parallel or mpi_partrac: ',mpi_parallel,mpi_partrac
   stop ! this is called before sepran has started
endif

if (mpi_partrac) mpi_parallel=.true.
if (mpi_parallel) then
   call MPI_INIT(ierr)  !! Start MPI
   if (ierr /= 0) then
      write(irefwr,*) 'PERROR(compr_start): MPI started with error: ',ierr
      stop
   endif
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)  !! Figure out my rank
   if (ierr /= 0) then
      write(irefwr,*) 'PERROR(compr_start): MPI_COMM_RANK reports error: ',ierr
      stop
   endif
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)  !! Figure out how many procs
   if (ierr /= 0) then
      write(irefwr,*) 'PERROR(compr_start): MPI_COMM_SIZE reports error: ',ierr
      stop
   endif
  call hostnm(hostname)
! write(irefwr,*) 'PINFO(compr_start): process ',myid,' of ', numprocs,' is alive with PID: ',getpid(),' on ',hostname
  ! make sure only process 0 provides output to stdout in most cases
  if (myid==0) then
     print_node=.true.
  else 
     print_node=.false.
  endif
! In parallel we to read from multiple instances of the same file rather than from stdin
  if (numprocs>999) then
     if (print_node) then
        write(irefwr,*) 'PERROR(check_parallel): numprocs > 999: ',numprocs
        write(irefwr,*) 'Fix length of index of input.XXX in check_parallel()'
     endif
     call instop
  endif
  lu_in=7
  lu_out=13
  if (myid==0) then
     ! sepran.env7 should have 'Unit number for reading of SEPRAN input' set to 7
     call system('rm -f input.* output.* sep_output.*')
     command='cp $SPHOME/sepran.env7 sepran.env'
     call system(command)
  endif
  ! make sure to have task 0 complete the above before continuing...
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  write(command,'(''cp comprP.in input.'',i3.3)') myid
  call system(command)
  write(fname,'(''input.'',i3.3)') myid
  open(lu_in,file=fname)

  ! reroute sepran output to 'output.000' etc.
  write(fname,'(''output.'',i3.3)') myid
  lu_out_sepran=14
  open(lu_out_sepran,file=fname)


  ! set up output files for all (we'll replace lu_out with 6 after sepran is initialized)
  write(fname,'(''sep_output.'',i3.3)') myid
  lu_out=13
  open(lu_out,file=fname)

  if (myid==0) call system('ls input.*')
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
else
#endif
   if (mpi_parallel.or.mpi_partrac) then
      write(irefwr,*) 'PERROR(check_parallel): mpi_parallel or mpi_partrac: ',mpi_parallel,mpi_partrac
      write(irefwr,*) 'but the code is compiled without -DMPI -fpp'
      call instop
   endif 
   myid=0
   numprocs=1
   write(irefwr,*) 'PINFO(compr_start): serial run'
   lu_in = 5
   print_node=.true.
   call system("rm -f sepran.env")
#ifdef MPI
endif
#endif
!if (print_node) write(irefwr,*) 'stopping at end of check_parallel'
!call instop
end subroutine check_parallel
