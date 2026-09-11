!   PEFBUOY1
!   calculate the buoyancy force in the Stokes equation.
!
!   Note: Sepran multiplies this with density rho to find
!   the complete buoyancy terms.
!
!   input:
!   ichoice     - select first/second component
!   x,y        - cartesian coordinates of the evaluation point
!   temp       - temperature value in the evaluation point
!
!   output: value of the buoyancy force component
real(kind=8) function pefbuoy1(ichoice,x,y,temp,pressure)
use sepmodulecomio
use coeff
use convparam
use control
implicit none
integer :: ichoice
real(kind=8) :: x,y,temp,pressure,adia
real(kind=8) :: aload,prespi,bigGamma,alpha_eff,funccf,z,rho,get_prespi
real(kind=8) :: bigGamma_norm,temp_norm,prespi_norm
integer :: iph
real(kind=8) :: r,cost,sint,temp_adiabatic,adiabat

! spherical coordinates
r = sqrt(x*x+y*y)
sint=x/r
cost=y/r

if (compress) then
   if (print_node) then
     write(irefwr,*) 'PERROR(pefbuoy1): should not be called'
     write(irefwr,*) 'for compressible convection'
   endif
   call instop
endif

! thermal buoyancy force
aload = Ra*temp

! Depth dependent expansivity; multiply with thermal part only
if (ialphatype>=1) then
   alpha_eff = funccf(5,x,y,z)
else
   alpha_eff = 1
endif
aload = aload * alpha_eff

if (compress.and..not.TALA) then
!  add effect of pressure term; has opposite sign of Ra*T
   aload = aload - DiG*pressure
   if (print_node) write(irefwr,*) 'aload: ',aload,pressure
endif


if (ichoice.eq.1) then
   pefbuoy1 = aload * sint
else if (ichoice.eq.2) then
   pefbuoy1 = aload * cost
else 
  if (print_node) write(irefwr,*) 'PERROR(pefbuoy1): unknown option: ',ichoice
  call instop
endif

end function pefbuoy1
