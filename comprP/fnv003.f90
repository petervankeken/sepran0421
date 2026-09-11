! FNV003
! Viscosity description with modelv=103
! Dummy version for isoviscous rheology
real(kind=8) function fnv003(x1,x2,x3,v1,v2,v3,secsqr,numold,maxunk,uold)
use sepmodulecomio
use coeff
use control
implicit none
integer :: numold,maxunk
real(kind=8) :: x1,x2,x3,v1,v2,v3,secsqr,uold(numold,maxunk)
integer :: jtypv,ival,ivalfind1
real(kind=8) :: temp,temp_a,pefvis

fnv003=1d0
if (numold<2) then
   if (print_node) then
       write(irefwr,*) 'PERROR(fnv003): numold < 2'
       write(irefwr,*) 'Second vector should contain temperature here'
   endif
   call instop
endif
jtypv=abs(itypv)
if (jtypv==0) then
   fnv003=1
else if (jtypv==1) then
   ival=ivalfind1(x1,x2)
   fnv003=viscl(ival)
else
   temp=uold(2,1)
   ival = ivalfind1(x1,x2)
   temp_a = 0.0_8
   fnv003=pefvis(x1,x2,temp,temp_a,ival,secsqr)
   !if (abs(x2)<1e-3) write(irefwr,*) 'fnv003: ',fnv003
endif

end function fnv003
