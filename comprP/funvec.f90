!   FUNCVEC
!
!   PvK 082808
subroutine funvec(rinvec,reavec,nvec,coor,outvec)
use sepmodulecomio
use geometry
use brandenburg
use eos
use convparam
use coeff
use control
implicit none
real(kind=8) :: rinvec(*),reavec(*),coor(*),outvec(*)
integer :: nvec
integer :: ioption,iph
real(kind=8) :: x,y,y1,r,T0,temp,bigGamma,prespi,get_prespi,th,cost,zdim,zbot_dim
real(kind=8) :: adiabatic_difference
real(kind=8), external :: ramp
integer icom
logical vcop
common /fv/ icom,vcop
   
ioption = nint(rinvec(1))
if (ioption< 1 .or. ioption>5) then
   if (print_node) write(irefwr,*) 'PERROR(funvec): 0< ioption < 5: ',ioption
   call instop
endif

   

if (cyl) then
   x = coor(1)
   y = coor(2)
   r = sqrt(x*x+y*y)
   y1 = r-r1
   cost=y/r
   th=acos(cost)
   if (coor(1)<0.0_8) then
      th=2.0*pi-th
   endif
else  
   y1 = coor(2)
endif
zdim=(1-y1)*height_dim
zbot_dim=height_dim

if (ioption.eq.1) then
!  Add adiabat Tbar to the solution vector
   if (iadiabat.eq.2) then
      adiabatic_difference=get_Tbar_dim(zdim)-get_Tbar_dim(zbot_dim)
      if (.not.delta1K) adiabatic_difference=adiabatic_difference/deltaT_dim
      outvec(1) = reavec(1)+adiabatic_difference
   else if (iadiabat.eq.3) then
      adiabatic_difference=get_Tbar_dim(zdim)
      if (.not.delta1K) adiabatic_difference=adiabatic_difference/deltaT_dim
      outvec(1) = reavec(1)+adiabatic_difference
   endif
endif

if (ioption.eq.2) then
   T0 = rinvec(2)
   if (Di.le.0d0.or.T0.le.0) then
      if (print_node) write(irefwr,*) 'PERROR(funvec): Di <= 0 or T0<=0'
      call instop
   endif
   y1 = coor(2)
   Temp = reavec(1)
   outvec(1) = log( (Temp+T0)/T0 ) - (1-y1)*Di
endif

if (ioption.eq.3) then
!  compute bigGamma
   temp = reavec(1)
   if (nph.gt.2) then
      if (print_node) then
         write(irefwr,*) 'PERROR(funvec): nph>2'
         write(irefwr,*) 'nph = ',nph
      endif
      call instop
   endif
   iph = 1
   bigGamma = 0
   do iph=1,nph
     !prespi = (1-y1) - phz0(iph) -gamma(iph)*(temp-pht0(iph))
     prespi=get_prespi(1-y1,temp,iph)
     bigGamma = bigGamma + 0.5 + 0.5*tanh(prespi/phdz(iph))
   enddo
   outvec(1) = bigGamma
endif

if (ioption == 4) then
   ! JP's functionality for jlimit
   if (icom==1) then
     outvec(1) = reavec(1)*ramp(th,5) + reavec(2)*ramp(th,6) + reavec(3)*ramp(th,7) + reavec(4)*ramp(th,8) + &
         & reavec(5)*ramp(th,1) + reavec(6)*ramp(th,2) + reavec(7)*ramp(th,3) + reavec(8)*ramp(th,4)
   endif
   if (icom==2) then
      outvec(1) = reavec(1)*ramp(th,5) + reavec(2)*ramp(th,6) + reavec(3)*ramp(th,7) + reavec(4)*ramp(th,8) + &
         & reavec(5)*ramp(th,1) + reavec(6)*ramp(th,2) + reavec(7)*ramp(th,3) + reavec(8)*ramp(th,4)
    endif
endif

end subroutine funvec

