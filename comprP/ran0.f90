real(kind=4) function ran0(idum)
integer :: idum,IA,IM,IQ,IR,MASK
parameter(IA=16807,IM=2147483647,AM=1e0/IM, IQ=127773,IR=2836,MASK=123459876)
integer :: k

idum = ieor(idum,MASK)
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
ran0 = AM*idum
idum = ieor(idum,MASK)
end function ran0

