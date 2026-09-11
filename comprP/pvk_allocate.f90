subroutine pvk_allocate(array,array_length,name_array)
use sepmodulecomio
use sepmodulecpack
use control
implicit none
real(kind=8),allocatable :: array(:)
integer,intent(in) :: array_length
character(len=*),intent(in) :: name_array
integer :: allocate_status

allocate(array(array_length),stat=allocate_status)
if (allocate_status /= 0 .and. print_node) then
   write(irefwr,*) 'allocate of array ',name_array(1:len_trim(name_array)),' failed'
   write(irefwr,*) 'length = ',length
   write(irefwr,*) 'allocate_status = ',allocate_status 
   call instop
endif

end subroutine pvk_allocate

