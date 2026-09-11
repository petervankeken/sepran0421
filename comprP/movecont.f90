! MOVECONT
! 
! Rotate the continent formation zone by 
subroutine movecont(toutstep)
use sepmodulecomio
use control
use geometry
implicit none
       
real(kind=8) :: toutstep
real(kind=8) :: dtheta
integer :: i

dtheta = 2*pi*toutstep/revol_cont
! Move second continent only
!      thmincont(2) = thmincont(2) + dtheta
!      thmaxcont(2) = thmaxcont(2) + dtheta
!      *** make sure theta of the continent doesn't
!      *** exceed 2*pi. 
if (thmaxcont(2).ge.2*pi) then
   if (revol_cont.lt.0) then
      ! *** continent revolves in counter clockwise fashion
      if (print_node) write(irefwr,*) 'PERROR(movecont): continent should move'
      if (print_node) write(irefwr,*) '    to larger theta values'
      call instop
   endif                
   dtheta = thmaxcont(2)- thmincont(2)
   thmaxcont(2) = thmaxcont(2)+ 0.1 - 2*pi
   thmincont(2) = thmincont(2)+ 0.1 - 2*pi
endif
if (print_node) write(irefwr,'(''Movecont '',i5,'' to: '',2f8.3)') 2,thmincont(2),thmaxcont(2)
 

end subroutine movecont
