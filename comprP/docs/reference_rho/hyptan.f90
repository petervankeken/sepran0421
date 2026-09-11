program hyptan
implicit none
integer :: i
real(kind=8) :: z0=400,z,w=10
do i=1,2885
   z=i*1
   write(6,*) z,tanh((z-z0)/w)
enddo
end program hyptan
