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
!use brandenburg
!use tracers
implicit none
real(kind=8) :: x,y,temp,secinv,func,adia,yh,zh
integer :: ival,jval
real(kind=8) :: z,E,V,t0,z0,A,A0=1.0_8,etanon,etalayer,etatemp,dr_smooth,r_smooth,smoothstep
real(kind=8) :: beta,secsqr,pow,radius,realtemp,y1,r
real(kind=8),parameter :: A_LZ=13.8155_8/3500.0_8  ! from Leng&Zhong, 2010 (page 210)
integer ::  ifirst,ivalfind,icounter
logical :: linear_only=.false.
data ifirst,icounter/0,1000/
save ifirst,icounter

pefvis=1.0_8
if (itypv == 0 ) then ! .or. itypv==10) then
   ! viscosity is constant !pvk or blob dependent
   !pvk! blob effect will be added in coefvis
   pefvis=1.0_8
   return
endif

A0=1.0_8
realtemp=temp
ival=1
if (cyl) then
   r=sqrt(x*x+y*y)
else
   r=y
endif

!pvk! find prefactor A0
!pvkif (weak_zones) then
!pvk    ! hack
!pvk    A0=1.0_8
!pvk!   if (cyl) then
!pvk!      if (print_node) write(6,*) 'PERROR(pefvis): weak_zones and cyl'
!pvk!      call instop
!pvk!   endif
!pvk    if (y>=0.95) then 
!pvk       A0=eta_plate ! lithosphere
!pvk        if (x<=0.05.or.x>=rlampix-0.05) A0=eta_wz
!pvk    endif
!pvk    write(6,'(''wz: '',3e15.7)') x,y,A0
!pvkelse 
!pvk   if (ivl.or.cyl.or.axi) ival = ivalfind(x,y)
!pvk   smoothstep_check: if (ivl.and..not.ivl_smoothstep) then
!pvk      A0 = viscl(ival)
!pvk      pefvis=A0
!pvk       !if (x<1e-3) write(irefwr,'(''pefvis ivl : '',4f15.7)') r,r_smooth,pefvis,r_zint(1)
!pvk   else if (ivl) then
!pvk      if (nlay>2) then
!pvk         if (print_node) write(irefwr,*) 'PERROR(pefvis): ivl & smoothstep needs work for nlay=',nlay
!pvk         call instop
!pvk      endif
!pvk      dr_smooth=0.03
!pvk      r_smooth=0.0
!pvk      if (.not.cyl) then
!pvk         if (print_node) write(irefwr,*) 'PERROR(pefvis): ivl & smoothstep needs work for cyl'
!pvk         call instop
!pvk      endif
!pvk      if (r>r_zint(1)+0.5*dr_smooth) then
!pvk         A0=viscl(1)
!pvk      else if (r<r_zint(1)-0.5*dr_smooth) then
!pvk         A0=viscl(2)
!pvk      else
!pvk         r_smooth=1.0_8-(r-(r_zint(1)-0.5*dr_smooth))/dr_smooth
!pvk         A0=viscl(1)+(viscl(2)-viscl(1))*smoothstep(r_smooth)
!pvk      endif
!pvk       if (x<1e-3) write(irefwr,'(''pefvis ivl smooth: '',4f15.7)') r,r_smooth,pefvis,r_zint(1)
!pvk   endif smoothstep_check
!pvkendif
!pvk!endif ! weak_zones
!pvkif (abs(itypv) == 1) then
!pvk   pefvis=A0
!pvk!   write(6,*) 'pefvis = A0 = ',A0
!pvk   return
!pvkendif
!pvk
!pvkif (delta1K) realtemp=realtemp/deltaT_dim
zh=1-y
if (cyl) then
   r=sqrt(x*x+y*y)
   zh=r2-r
endif

if (itypv == 2 .or. itypv==3) then
   !! Default Blankenbach et al. case 2a
   pefvis=A0*exp(-b_eta*realtemp+c_eta*zh)
!pvk   if (gable_plates .and. ival==1) pefvis=A0
!  if (abs(x)<1e-3) write(irefwr,'(''eta: '',i5,7f12.6)') ival,x,y,r,realtemp,b_eta,A0,pefvis
   return
endif

!pvkif ((itypv == 4 .or. itypv==5) .and. .not.tosi15) then
!pvk   if (print_node) then
!pvk      write(irefwr,*) 'PERROR(pefvis): tosi15 should be T'
!pvk      write(irefwr,*) 'when itypv > 3 for this implementation'
!pvk   endif
!pvk   call instop
!write(6,*) 'itypv, tosi_case15: ',itypv,tosi15_case
if ((itypv == 4.or.itypv==5) .and. tosi15) then
   etalinh=exp(-b_eta*realtemp+c_eta*zh)
   ! if ((icounter/1000)*1000==icounter.and.printvis) then
   !        write(6,'(''etalinh1: '',5e15.7)') x,y,etalinh,b_eta,realtemp
   !  endif
       if (tosi15_case == 1 .or. tosi15_case == 3) then
           pefvis = A0*etalinh
           return
       endif
!      write(6,*) 'secinv: ',secinv
       if (secinv>1e-8_8) then
           ! note sqrt(2) term is necessary to correct for
           ! difference of definition of strainrate and invariant
           ! between Tosi15 and sepran.
           etaplasth = eta_star + sqrt(2.0)*sigma_y/sqrt(secinv)
       else
           etaplasth = 1e8  ! arbitrarily high
       endif
    
    
       if (ival==1) then
          linear_only=.false.
          pefvis = 2.0_8/(1.0_8/etalinh + 1.0_8/etaplasth)
       else 
          linear_only=.true.
          pefvis = 2.0_8*etalinh
       endif
       printvis=.false.
       if ((icounter/1000)*1000==icounter.and.printvis) then
          write(6,'(''etalinh2: '',6e15.7,l5)') x,y,etalinh,etaplasth,pefvis,realtemp,linear_only
       endif
!pvk
!pvk    if (tackley .and. tosi15) then
!pvk       ! write(6,*) 'tackley!'
!pvk       ! special case of Tosi benchmark with use of Tackley
!pvk       ! rheological :: parameters for MN.
!pvk
!pvk       etalinh=exp(etalin_z_dep*zh)*exp(27.631/(realtemp+1.d0-0.12))*exp(-27.631/(0.64+1.d0-0.12))*5.706e-6*visc_factor
!pvk
!pvk
!pvk       !!!!etalinh=A0*etalinh
!pvk
!pvk       etaplasth = eta_star + (sigma_y+sigma_b*zh)/ (sqrt(2.d0)*sqrt(secinv))
!pvk       if ((icounter/1000)*1000==icounter .and. printvis) then
!pvk           write(6,'(''eta_star etc: '',8e15.7)') eta_star,sigma_y,sigma_b,secinv,etaplasth,zh,z_yieldlimit,A0
!pvk       endif
!pvk
!pvk       if (zh<z_yieldlimit) then
!pvk          pefvis = min(etalinh,etaplasth)
!pvk       else
!pvk          etaplasth=1.0e3 ! for output purposes
!pvk          pefvis = etalinh
!pvk       endif
!pvk       ! return NB do not return to test for negative viscosity 
!pvk
!pvk    else if (tackley_newvis) then !Tosi style viscosity average
!pvk
!pvk       etalinh=exp(4.6d0*zh)*exp(27.631/(realtemp+1.d0-0.12))*exp(-27.631/(0.64+1.d0-0.12))*1.8e-3*0.4855*visc_factor
!pvk
!pvk       etalinh=A0*etalinh
!pvk
!pvk       etaplasth = eta_star + (sigma_y+sigma_b*zh)/(sqrt(2.d0)*sqrt(secinv))
!pvk
!pvk       if (zh<z_yieldlimit) then
!pvk          pefvis = 1d0/(1d0/etalinh + 1d0/etaplasth)
!pvk       else
!pvk          etaplasth=1.0e3 ! for output purposes
!pvk          pefvis=etalinh
!pvk       endif
!pvk       ! return NB do not return to test for negative viscosity
!pvk    else if (tackley2000) then
!pvk         !viscosity from Tackley 2000
!pvk         etalinh=exp(23.03d0*(1d0/(realtemp+1d0)-0.5d0)) *exp(-23.03d0*(1d0/(0.64+1d0)-0.5d0))*5.21e-3*visc_factor
!pvk         !This is set to 1 at realtemp=Ts_T
!pvk         etalinh=A0*etalinh
!pvk         etaplasth = eta_star + (sigma_y+sigma_b*(1.0-y))/(sqrt(2d0)*sqrt(secinv))
!pvk
!pvk         if (zh<z_yieldlimit) then
!pvk           pefvis = min(etalinh,etaplasth)
!pvk         else
!pvk           etaplasth=1.0e3 ! for output purposes
!pvk           pefvis=etalinh
!pvk         endif
!pvk       ! return NB do not return to test for negative viscosity
!pvk
!pvk      if ((icounter/1000)*1000==icounter.and.printvis) then
!pvk         write(6,'(''etalinh3: '',5e17.6,3L5)') x,y,etalinh,etaplasth,pefvis,tackley,tackley_newvis,tackley2000
!pvk      endif
!pvk    endif   ! extra Tackley options
!pvk
   endif
!pvk
icounter=icounter+1

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
