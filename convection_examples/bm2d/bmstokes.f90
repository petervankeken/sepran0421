! PvK July 2024
subroutine bmstokes()
use sepmoduleoldrouts
use sepmodulecomio
use sepran_arrays
use control
use coeff
implicit none
integer :: build_in(20),solve_in(20),inpsol(20),iread
real(kind=8) :: solve_rin(20),vrms
integer :: ichoice
logical :: first=.true.
save first

ichoice=1
call fillcoef900(ichoice)
build_in=0
if (first.or.itypv>1) then
   build_in(1) =10 ! maximum number of entries
   build_in(2) = 1 ! build matrix and rhsd vector
   build_in(10)= 3 ! number of old vectors
   first=.false.
else
   build_in(1) =10 ! maximum number of entries
   build_in(2) = 2 ! build just rhsd vector
   build_in(10)= 3 ! number of old vectors
endif
call build(build_in,matr(1),intmat(1),kmesh,kprob,irhsd(1),matrm(1),isol(1),islold(1),iuser_here,user_here)

inpsol=0
solve_rin=0.0_8
if (isolmethod9==-1) then
   call solve(1,matr(1),isol(1),irhsd(1),intmat(1),kprob)
else if (isolmethod9==0) then
   inpsol(1)=3
   inpsol(2)=1 ! profile method
   inpsol(3)=0
else
   if (print_node) write(irefwr,*) 'PERROR(bmstokes): unsuited yet for isolmethod<>0'
   call instop
endif
iread=-1
call solvel(inpsol,solve_rin,matr(1),isol(1),irhsd(1),intmat(1),kmesh,kprob,iread)
pedebug=.false.
if (pedebug) then
   call pevrms(vrms)
   write(6,*) 'in bmstokes: ',vrms
endif
pedebug=.false.
   
end subroutine bmstokes
