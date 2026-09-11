subroutine getsepraninfomatr(matr)
use sepmodulecomio
use sepmodulematr
use sepmoduleintmat
use sepmodulematstr
use sepmodulekprob
use control
implicit none
integer, intent(in) :: matr

if (print_node) write(irefwr,*) 'iprob, jmethod: ',realmt(matr)%iprob,realmt(matr)%jmethod,realmt(matr)%decomp

end subroutine getsepraninfomatr
