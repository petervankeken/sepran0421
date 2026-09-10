! PvK July 2024
subroutine steady_heatP
use sepran_arrays
use sepmoduleoldrouts
implicit none
integer :: build_in(20)=0,solve_in(20)=0,iread
real(kind=8) :: solve_rin(20)

call fillcoef800()

build_in=0
build_in(1)=9 ! maximum number of entries
build_in(2)=1 ! build matrix and rhsd vector
build_in(9)=2 ! IPROB
call build(build_in,matr(2),intmat(2),kmesh,kprob,irhsd(2),matrm(2),isol(2),isol(2),iuser_here,user_here)

solve_in=0
solve_rin=0.0_8
solve_in(1)=3 ! max number of entries
solve_in(3)=0 ! direct solution method
iread=-1 ! indicate information is provide by array solve_in
call solvel(solve_in,solve_rin,matr(2),isol(2),irhsd(2),intmat(2),kmesh,kprob,iread)

end subroutine steady_heatP
