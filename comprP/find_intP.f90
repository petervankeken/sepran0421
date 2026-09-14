subroutine find_intP(plasticity_integral)
use sepmodulecomio
use sepmodulekmesh
use sepran_arrays
use sepmodulevecs
use coeff
use geometry
use control
implicit none
real(kind=8),intent(out) :: plasticity_integral
integer :: allocate_status,i
real(kind=8) :: outval(3)
integer :: iinvol(5)


if (.not.allocated(user_here)) then
    allocate(user_here(npoint+5),stat=allocate_status)
    if (allocate_status/=0) then
       if (print_node) then
          write(irefwr,*) 'PERROR(find_intP): allocation of array user_here failed: ',allocate_status
       endif
       call instop
       user_here(1)=1.0_8*(npoint+5)
    endif
endif

iinvol=0

! put plasticity vector into user_here
call find_intP01(npoint,ks(iplasticity)%sol,user_here(6))
!iuser_here(2)=1
iuser_here(3)=0
iuser_here(4)=0
iuser_here(5)=0
iuser_here(6)=7
iuser_here(7)=intrule900
iuser_here(8)=icoor900
iuser_here(9)=0
iuser_here(10)=2001
iuser_here(11)=6

call integr(iinvol,outval,kmesh1,kprob1,iplasticity,iuser_here,user_here)
!deallocate(user_here)

plasticity_integral=outval(1)/volume

end subroutine find_intP


subroutine find_intP01(npoint,plasticity,user)
implicit none
integer,intent(in) :: npoint
real(kind=8),intent(in) :: plasticity(*)
real(kind=8),intent(out) :: user(*)

user(1:npoint)=plasticity(1:npoint)

end subroutine find_intP01
