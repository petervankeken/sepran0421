subroutine set_up_chempost_logs
use sepmodulecomio
use brandenburg
use control
use mparallel
implicit none

integer :: i
character(len=80) :: fname

if (mpi_parallel) then
   if (print_node) write(irefwr,*) 'PERROR(set_up_chempost_logs): not yet suited for mpi_parallel'
   call instop
endif

if (nplates>10) then
   if (print_node) write(irefwr,*) 'PERROR(set_up_chempost_logs): nplates>10'
   call instop
endif

if (myid==0) then
  call system('mkdir -p chemistry')
  do i=1,nplates
     write(fname,'(''chemistry/outgas.'',i2.2)') i
     open(LU_MELT_START+i,file=fname)
  enddo
  do i=1,nplates
     write(fname,'(''chemistry/ingas.'',i2.2)') i
     open(LU_INGAS_START+i,file=fname)
  enddo
  if (track_cmb_ingas) then
    open(LU_CMB_ENTRY,file='chemistry/cmb_entry.log')
    open(LU_CMB_EXIT,file='chemistry/cmb_exit.log')
  endif
  open(LU_TAU_DROP,file='tau_drop.dat')
  open(LU_TSTEP,file='chemistry/tstep.log')
  !open(301,file='atracer.dat')
endif

end subroutine set_up_chempost_logs
