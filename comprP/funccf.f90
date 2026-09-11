real(kind=8) function funccf(ichoice,x,y,z)
use sepmodulecomio
use control
use convparam
use geometry
use coeff
use eos
implicit none
integer :: ichoice
real(kind=8) :: x,y,z,r,zr,zdim,rho_nd
real(kind=8) :: temp,drho_ndz

funccf=0.0_8
if (ichoice==10) then
   funccf=0
   return
endif

if (cyl) then
   r = sqrt(x*x+y*y)
   zr = r2-r
else
   zr = 1-y
endif
zdim=zr*height_dim

if (ichoice==1) then
   temp = 1-y + 0.1*sin(pi*y)*cos(pi*x)
   funccf = 1e4*temp

else if (ichoice==3) then

   ! density
   funccf=get_rhobar_dim(zdim)/rho_dim

else if (ichoice==4) then
   ! conductivity
   rho_nd=get_rhobar_dim(zdim)/rho_dim
 
   if (icondtype==1) then
      ! as in PEPI01
      funccf=(rho_nd/rho_av)**pcond
   else if (icondtype==2) then
      ! get rid of dependence on <rhobar>
      funccf=rho_nd**pcond
   else
      funccf=1.0_8
   endif
!   write(irefwr,*) zdim,zr,rho_nd

else if (ichoice==5) then
   ! thermal expansivity
   funccf=exp(-zr)
   funccf=get_alpha_dim(zdim,exp(-zr))/alpha_dim

else if (ichoice == 7) then
   ! compositional compressibility beta
   if (ibetatype == 11) then
!     compressibility beta as used in CH94
      funccf = (CH94_s*exp(-CH94_s*zr))/(1d0-exp(-CH94_s))
   else if (ibetatype == 12) then
!     compressibility beta as used in CH94 doubled
      funccf = (CH94_s*exp(-CH94_s*zr))/(1d0-exp(-CH94_s))
   else if (ibetatype == 2) then
!     compressibility beta as used in CH94 with beta=1 at top
      funccf = exp(-CH94_s*zr)
   else
      funccf = 1d0
   endif

else if (ichoice==13) then
   ! d rho / dz 
   funccf=drho_ndz(zr)

else
   if (print_node) write(irefwr,*) 'PERROR(funccf): incorrect option ichoice = ',ichoice
   call instop
endif

end function funccf

!real(kind=8) function rho_n(z)
!use control
!use coeff
!use convparam
!implicit none
!real(kind=8) :: z
!if (compress) then
! if (fake_rho_bar) then
!   if (step_rho_background) then
!     rho_n=1.0_8
!     if (z>0.5_8) rho_n=1.0_8+drho_background_dense
!   else  if (exp_rho_background) then
!     rho_n=exp(0.5*z)
!   else
!     rho_n=1+drho_background_dense*z
!!   endif
! else
!     ! Jarvis & McKenzie, 1980
!     rho_n = exp(Di*z/Grueneisen)
!  endif
!else if (ibench_type==2) then
!     rho_n = exp(Di*z/Grueneisen) ! PEPI01
!else
!   rho_n = 1d0
!endif
!!write(irefwr,*) 'fake_rho_bar: ',compress,fake_rho_bar,step_rho_background,z,rho_n,drho_background_dense
!
!end function rho_n

real(kind=8) function drho_ndz(z)
use coeff
use convparam
implicit none
real(kind=8) :: z
if (compress) then
   ! Jarvis & McKenzie, 1980
   drho_ndz = Di/Grueneisen*exp(Di*z/Grueneisen)
else
   drho_ndz = 1d0
endif

end function drho_ndz

!real(kind=8) function rho_a(z)
!use coeff
!use convparam
!implicit none
!real(kind=8) :: z
!
!if (compress) then
!   rho_a=(-DiG*palpha*z+1)**(-1.0_8/palpha)
!else
!   rho_a=1.0_8
!endif
!
!end function rho_a
