! PvK July 2024
real(kind=8) function func(ichoice,x,y,z)
use convparam
use geometry
implicit none
integer :: ichoice
real(kind=8) :: x,y,z

func=0.0_8
if (ichoice==1) then
   ! to make compatible with skye
   func=1.0_8-y+0.1*sin(pi*y)*cos(pi*x/wavel_perturb)
else if (ichoice==2) then
   ! horizontal component of velocity
   func = -40*pi*sin(pi*x)*cos(pi*y)
else if (ichoice==3) then
   ! vertical component of velocity
   func = 40*pi*sin(pi*y)*cos(pi*x)
else 
   write(6,*) 'PERROR(func): ichoice <> 1-3: ',ichoice
   call instop
endif


end function func


