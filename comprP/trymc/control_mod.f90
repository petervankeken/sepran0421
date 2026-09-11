module control
character(len=150) :: checkmem_command
contains
   subroutine get_checkmem_command
   !use control
!#ifdef IFORT
!use IFPORT
!#endif
   implicit none
   integer :: iret
   write(checkmem_command,'(''echo 0 $(awk '',a1,''/Rss/ {print "+", $2}'',a1,'' /proc/'',I0,''/smaps) | bc > resmem'')') & 
          & "'","'",getpid()
   write(6,'(a120)') checkmem_command

   end subroutine get_checkmem_command

   real(kind=8) function get_resmem()
!#ifdef IFORT
!use IFPORT
!#endif
   !use control
   implicit none
   integer :: iret,isize
   iret=system(checkmem_command)
   open(9,file='resmem')
   read(9,*) isize
   close(9)
   get_resmem=isize*1.0_8/1024/1024
   end function get_resmem

end module control
