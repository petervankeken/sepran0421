subroutine handle_err(errcode)
use sepmodulecomio
use control
implicit none
include 'netcdf.inc'
integer errcode
if (print_node) write(irefwr,*) 'netcdf error: ',nf_strerror(errcode)
call instop
end subroutine handle_err
