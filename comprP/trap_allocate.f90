subroutine trap_allocate(alloc_status,subname,arrayname)
use sepmodulecomio
use control
implicit none
integer :: alloc_status
character(len=*) :: subname,arrayname

if (alloc_status/=0) then
   if (print_node) then
      write(irefwr,*) 'PERROR: in subroutine ',subname
      write(irefwr,*) 'error allocating array ',arrayname
      write(irefwr,*) 'alloc_status = ',alloc_status
   endif
   call instop
endif


end subroutine trap_allocate
