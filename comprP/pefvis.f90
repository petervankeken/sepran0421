! PEFVIS
!  
! Version for sigma,T dependent rheology in layered case
! Find value of ival through /layer/zint(10).
!  
! Since ival is modified inside this version of pefvis
! it should be treated as a dummy argument without changing it.
!  
! PvK 032910: Modified for compressible convection. Add adia
! which should have local adiabatic temperature
!  
! PvK 082511: addition of delta1K (where temperature is in K)
real(kind=8) function pefvis(x,y,temp,adia,jval,secinv)
use sepmodulecomio
use convparam
use coeff
use control
use geometry
use brandenburg
implicit none
real(kind=8) :: x,y,temp,secinv,func,adia,yh,zh
integer :: ival,jval
real(kind=8) :: z,E,V,t0,z0,A,A0,etanon,etalayer,etatemp,dr_smooth,r_smooth,smoothstep
real(kind=8) :: beta,secsqr,pow,radius,realtemp,y1,r
real(kind=8),parameter :: A_LZ=13.8155_8/3500.0_8  ! from Leng&Zhong, 2010 (page 210)
integer ::  ifirst,ivalfind,icounter
logical :: printvis=.false.
data ifirst,icounter/0,0/
save ifirst,icounter

pefvis=1.0_8
if (itypv == 0 .or. itypv==10) then
   ! viscosity is constant or blob dependent
   ! blob effect will be added in coefvis
   pefvis=1.0_8
   return
endif

realtemp=temp
ival=1
if (cyl) then
   r=sqrt(x*x+y*y)
else
   r=y
endif
if (ivl.or.cyl.or.axi) ival = ivalfind(x,y)
smoothstep_check: if (ivl.and..not.ivl_smoothstep) then
   A0 = viscl(ival)
else if (ivl) then
   if (nlay>2) then
      if (print_node) write(irefwr,*) 'PERROR(pefvis): ivl & smoothstep needs work for nlay=',nlay
      call instop
   endif
   dr_smooth=0.03
   r_smooth=0.0
   if (.not.cyl) then
      if (print_node) write(irefwr,*) 'PERROR(pefvis): ivl & smoothstep needs work for cyl'
      call instop
   endif
   if (r>r_zint(1)+0.5*dr_smooth) then
      A0=viscl(1)
   else if (r<r_zint(1)-0.5*dr_smooth) then
      A0=viscl(2)
   else
      r_smooth=1.0_8-(r-(r_zint(1)-0.5*dr_smooth))/dr_smooth
      A0=viscl(1)+(viscl(2)-viscl(1))*smoothstep(r_smooth)
   endif
    !if (x<1e-3) write(irefwr,'(''pefvis ivl smooth: '',4f15.7)') r,r_smooth,pefvis,r_zint(1)
endif smoothstep_check
if (abs(itypv) == 1) then
   pefvis=A0
   return
endif

if (delta1K) realtemp=realtemp/deltaT_dim
zh=1-y
if (cyl) then
   r=sqrt(x*x+y*y)
   zh=r2-r
endif

if (itypv == 2 .or. itypv==3) then
   !! Default Blankenbach et al. case 2a
   pefvis=A0*exp(-b_eta*realtemp+c_eta*zh)
   if (gable_plates .and. ival==1) pefvis=A0
   !if (abs(x)<1e-3) write(irefwr,'(''eta: '',i5,7f12.6)') ival,x,y,r,realtemp,b_eta,A0,pefvis
endif

if (pefvis <= 0.0e0_8) then
   if (print_node) then
      write(irefwr,*) 'PERROR(pefvis): viscosity <= 0'
      write(irefwr,*) 'itypv, ivl: ',itypv,ivl
      write(irefwr,*) 'x,y,cyl,axi ',x,y,cyl,axi
      write(irefwr,*) 't,ival,rt   ',temp,ival,realtemp
      write(irefwr,*) 'secinv      ',secinv,secsqr
      write(irefwr,*) 'viscl(ival) ',viscl(ival)
      write(irefwr,*) 'A0,b,c    : ',A0,-b_eta,c_eta
      write(irefwr,*) 'pefvis      ',pefvis
   endif
   call instop
endif
    
end function pefvis

real(kind=8) function smoothstep(x)
implicit none
real(kind=8),intent(in) :: x

if (x<0) then
   smoothstep=0.0_8
else if (x>1) then
   smoothstep=1.0_8
else
   smoothstep=3*x*x-2*x*x*x
endif

end function smoothstep
