function funcbc( ifunc, x, y, z)
use sepmodulecomio
use mtime
use coeff
use convparam
use geometry
use brandenburg
use control
implicit none
!include 'SPcommon/ctimen' ! for time t

real(kind=8) :: funcbc,x,y,z,vel,omega,v_left,v_right
real(kind=8) :: t_here
integer :: ifunc
real(kind=8), parameter :: Delta_x=1.0_8 ! CH94 equations 10+11
real(kind=8) :: x_c,vel_right,vel_left,x_l,r,th
real(kind=8),external :: ramp

omega = u_nought*pi/5d0
funcbc=0.0_8

if (ifunc==3) then
   if (.not.gable_plates.or.gable_stokes_choice==0) then
      funcbc=0.0_8
      return
   endif
   if (cyl) then
      r=sqrt(x*x+y*y)
      th=acos(y/r)
      if (x<=0.0_8) then
         th=2.0_8*pi-th
      endif
      funcbc=ramp(th,gable_stokes_choice)
      if (periodic) then
         if (gable_stokes_choice==1) then
            funcbc=ramp(th,1)+ramp(th,nplates)
         else
            funcbc=ramp(th,gable_stokes_choice)
         endif
      else
         funcbc=ramp(th,gable_stokes_choice)
      endif ! periodic
   else ! .not.cyl
      funcbc=ramp(x,gable_stokes_choice)
   endif ! cyl
   ! write(109,'(3f15.7)') x,y,funcbc
   ! write(irefwr,'(''bc: '',i5,5f15.7)') gable_stokes_choice,x,y,r,th,funcbc
     

else if (ifunc == 10) then
   ! PvK implementation of CH94 boundary condition
   x_c = x_nought+Delta_x*cos(omega*time_now) ! CH94 equation 11
   vel_left=u_nought + omega*Delta_x*0.5*sin(omega*time_now)
   vel_right=-u_nought + omega*Delta_x*0.5*sin(omega*time_now)
   ! overwrite to get downwelling in center
!  x_c = 2
!  vel_left = u_nought
!  vel_right = -u_nought
   if (x<1e-9) then
      funcbc=0
   else if (x>xcmax-1e-9) then
      funcbc=0
   else if (x<=x_c-dx_two_elements) then 
      funcbc = vel_left
!     if (x<dx_two_elements) then
!        x_l=x/(2*dx_two_elements)
!        funcbc = vel_left * (6*x_l**5-15*x_l**4+10*x_l**3)
!     endif
   else if (x>=x_c+dx_two_elements) then
      funcbc = vel_right
!     if (x>xcmax-dx_two_elements) then
!        x_l=(xcmax-x)/(2*dx_two_elements)
!        funcbc = vel_right * (6*x_l**5-15*x_l**4+10*x_l**3)
!     endif
   else 
      x_l = (x - (x_c-dx_two_elements))/(2*dx_two_elements)
      funcbc = vel_left + (6*x_l**5-15*x_l**4+10*x_l**3)*(vel_right-vel_left)  !  smootherstep=6x^5-15x^4+10x^3
   endif
!  write(irefwr,'(3f12.3)') x,funcbc

else if (ifunc == 12) then
   ! ??
   funcbc=100d0

else if (ifunc == 11) then
   if (print_node) write(irefwr,*) 'PERROR(funcbc) :: ichoice=11 is obsolete'
   call instop
endif

end function funcbc
