program adia
implicit none
real(kind=8) :: Ts,Di,Tbars,deltaTr,deltaTbar,deltaTp,z,Ttot,Tp,Tbar,y
integer :: i,j
character(len=80) :: fname

Ts=273
deltaTr=3000
Tbars=1600
Di=0.5
do j=1,4 
   Di=0.25*j
   write(fname,'(''Tbar.'',i3.3)') nint(Di*1e2)
   open(9,file=fname)
   deltaTbar=Tbars*(exp(Di)-1)
   deltaTp=deltaTr-deltaTbar
   write(6,*) deltaTbar,deltaTp
   do i=1,101
      y=(i-1)*0.01
      z=1-y
      Ttot=Ts+deltaTr*(1-y)
      Tbar=Tbars*exp(Di*z)
      Tp=Ttot-Tbar+Tbars-Ts
      write(9,'(4f15.7)') y,Ttot,Tbar,Tp
   enddo
   close(9)
enddo
      


end program adia
