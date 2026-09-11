subroutine check_stress_tensor(iuser,user)
use sepmoduleoldrouts
use sepmodulecomio
use sepmodulekmesh
use sepmodulevecs
use coeff
use control
use sepran_arrays
implicit none
integer :: iuser(*)
real(kind=8) :: user(*)
integer :: igradv,istress,mystress
integer :: iinder(200)
real(kind=8) :: stress11my,maxpressure
real(kind=8) :: gradv11,gradv12,gradv22,div11,stress11,stress12,stress22
real(kind=8) :: stress21,gradv21,div22
integer :: ihelp,nparm,iwork(4,100),i,iinvec(10)
real(kind=8) :: work(100),rinvec(10),visdip(2),visdip1(2),getrmsdiv

call sepactsolbf1(isol1)
call pefilcof(1,iuser,user)
call coefvis
   
iinder(1)=8
iinder(2)=1
iinder(3)=0
iinder(5)=0
iinder(6)=0
iinder(7)=0
iinder(8)=2

if (print_node) write(irefwr,*) 'compute grad v'
iinder(4)=2
call deriv(iinder,igradv,kmesh1,kprob1,isol1,iuser,user)

iinder(4)=4
if (print_node) write(irefwr,*) 'compute div v'
call deriv(iinder,idivv,kmesh1,kprob1,isol1,iuser,user)

do i=1,100
   iwork(1,i)=0
   iwork(2,i)=0
   iwork(3,i)=0
   iwork(4,i)=0
    work(i) = 0d0
enddo
iwork(1,1) = 0
iwork(1,2) = 1
iwork(1,3) = intrule900 
iwork(1,4) = icoor900
iwork(1,5) = mcontv
if (itype_stokes==900) then
!  penaltyfunction parameter
   eps = 1d-6
else
   eps = 0
endif
!  6: epsilon
iwork(1,6) = 0
 work(6) = eps
!  7: rho
if (compress) then
   iwork(1,7) = 3
else 
   iwork(1,7) = 0
    work(7)   = 1d0
endif
!  12: eta assumed constant 1 here
iwork(1,12) = 0
 work(12)   = 1d0

do i=6,iuser(1)
   iuser(i)=0
enddo
do i=6,nint(user(1))
   user(i)=0d0
enddo
nparm=12
!if (print_node) write(irefwr,*) 'bmout line 561'
call fil103(1,1,iuser,user,kprob1,nparm,iwork,work,kmesh1)
if (print_node) write(irefwr,*) 'compute stress tensor'
! icheld = 6 for the stress tensor
iinder(4)=6
call deriv(iinder,istress,kmesh1,kprob1,isol1,iuser,user)

! if (print_node) write(irefwr,*) 'compute viscous dissipation'
iinder(4)=11
call deriv(iinder,iphi,kmesh1,kprob1,isol1,iuser,user)

gradv11 = anorm(1,1,1,kmesh1,kprob1,igradv,igradv,ihelp)
gradv22 = anorm(1,1,4,kmesh1,kprob1,igradv,igradv,ihelp)
gradv12 = anorm(1,1,2,kmesh1,kprob1,igradv,igradv,ihelp)
gradv21 = anorm(1,1,3,kmesh1,kprob1,igradv,igradv,ihelp)
div11   = anorm(1,1,1,kmesh1,kprob1,idivv,idivv,ihelp)
if (eos_type==0) then
   ! compute sqrt (int u^2) in a roundabout way (ichnorm=7 in anorm or equivalent in manvec is not available)
   div22 = getrmsdiv(kmesh1,kprob1,isol1,ks(idivv)%sol,iuser,user)
endif


stress11= anorm(1,1,1,kmesh1,kprob1,istress,istress,ihelp)
stress22= anorm(1,1,2,kmesh1,kprob1,istress,istress,ihelp)
stress12= anorm(1,1,4,kmesh1,kprob1,istress,istress,ihelp)
stress21= anorm(1,1,4,kmesh1,kprob1,istress,istress,ihelp)
visdip = anorm(1,1,1,kmesh1,kprob1,iphi,iphi,ihelp)

if (print_node) then 
   write(irefwr,*) 'NORMS: '
   write(irefwr,*) 'grad v = ',gradv11,gradv22,gradv12
   write(irefwr,*) 'div  v = ',div11,div22
   write(irefwr,*) 'stress = ',stress11,stress22,stress12
   write(irefwr,*) 's11:     ',2*gradv11-2d0/3*div11,stress11
   write(irefwr,*) 's22:     ',2*gradv22-2d0/3*div11,stress22
   write(irefwr,*) 'visdip:  ',visdip
   if (itype_stokes==903) then
      maxpressure= anorm(1,1,3,kmesh1,kprob1,isol1,isol1,ihelp)
      write(irefwr,*) 'max pressure:',maxpressure
   else if (itype_stokes==900) then
      maxpressure= anorm(1,1,1,kmesh1,kprob1,idivv,idivv,ihelp)
      write(irefwr,*) 'max div u:',maxpressure
   endif
call checkstress001(kmesh1,kprob1,isol1,idivv,igradv,istress,iphi,ivisc)
endif



end subroutine check_stress_tensor

subroutine checkstress001(kmesh,kprob,isol,idivv,igradv,istress,iphi,ivisc)
use sepmodulecomio
use sepmodulekmesh
use sepmodulekprob
implicit none
integer :: kmesh,kprob,isol,idivv,igradv,istress
integer :: iphi,ivisc
integer :: ndef_div,ndef_grad,ndef_stress,i
integer :: ipkprf,ipvisdip,ndef_visdip,ipvisc

indprf = indprfi
!if (indprf.ne.0) call ini070(indprf)
!ndef_div = idivv(5)/npoint
!ndef_grad = igradv(5)/npoint
!ndef_stress = istress(5)/npoint
!ndef_visdip = iphi(5)/npoint
!write(irefwr,*) 'ndef: ',npoint,ndef_div,ndef_grad,ndef_stress

write(irefwr,*) idivv,igradv,istress,iphi,ivisc
write(irefwr,*) npoint,nunkp,indprf
write(irefwr,*) ks(ivisc)%sol(1),ks(iphi)%sol(1)
call checkstress002(coor,npoint,ks(idivv)%sol,ks(igradv)%sol,ks(istress)%sol,nunkp, &
    &          kprobfi,indprf,ks(iphi)%sol,ks(ivisc)%sol)
end subroutine checkstress001

subroutine checkstress002(coor,npoint,div,grad,stress,nunkp,kprobf,indprf,visdip,visc)
use sepmodulecomio
use convparam
use coeff
implicit none
integer :: npoint,indprf,kprobf(*),nunkp
real(kind=8) :: coor(2,npoint),div(1,npoint),grad(4,npoint)
real(kind=8) :: stress(6,npoint),sep_stress11,dstress11,stress11
real(kind=8) :: visdip(*),visc(*)  ! ,uvp(*)
integer :: i,j1
real(kind=8) :: stress11tot,sep_stress11tot,divtot,dstress11tot
real(kind=8) :: presstot,press,visdiptot,dvisdip,visdip1
real(kind=8) :: sep_visdiptot,sep_visdip,x,y,pefvis
real(kind=8) :: s11max(2),visdipmax(2)

stress11tot=0
divtot=0
dstress11tot=0
sep_stress11tot=0
presstot=0
visdiptot=0
sep_visdiptot=0
dvisdip=0
write(irefwr,*) 'compress = ',compress,npoint
s11max=0
visdipmax=0
do i=1,npoint
   if (.not.compress) then
      visdip1 = 4*grad(1,i)*grad(1,i)+4*grad(4,i)*grad(4,i)
      visdip1 = visdip1+2*(grad(2,i)+grad(3,i))*(grad(2,i)+grad(3,i))
      visdip1 = 0.5*visc(i)*visdip1
   else 
      visdip1 = (4d0/3*grad(1,i)-2d0/3*grad(4,i))**2
      visdip1 = visdip1 + 2*(grad(2,i)+grad(3,i))**2
      visdip1 = visdip1 + (4d0/3*grad(4,i)-2d0/3*grad(1,i))**2
      visdip1 = 0.5*visc(i)*visdip1
   endif
   sep_visdip = visdip(i)
   visdipmax(1)=visdipmax(1)+visdip(i)
   visdipmax(2)=visdipmax(2)+visdip1
   dvisdip = dvisdip+(sep_visdip-visdip1)
   !write(irefwr,*) i,visdip1,dvisdip
   visdiptot = visdiptot+visdip1
   sep_visdiptot = sep_visdiptot+sep_visdip

   stress11 = 2*grad(1,i) - 2d0/3d0*div(1,i)
   sep_stress11 = stress(1,i)
   s11max(1)=s11max(1)+abs(sep_stress11)
   s11max(2)=s11max(2)+abs(stress11)
   dstress11 = (sep_stress11-stress11)
   stress11tot = stress11tot+abs(stress11)
   sep_stress11tot = sep_stress11tot+abs(sep_stress11)
   dstress11tot = dstress11tot+abs(dstress11)
   divtot = divtot+abs(div(1,i))
   !if (indprf.eq.0) then
   !   j1 = (i-1)*nunkp+3
   !else 
   !   j1 = kprobf(i)+3
   !endif
   !press = uvp(j1)
   !presstot = presstot + abs(press)
enddo
stress11tot = stress11tot/npoint
sep_stress11tot = sep_stress11tot/npoint
dstress11tot = dstress11tot/npoint
divtot = divtot/npoint
presstot = presstot/npoint
write(irefwr,'(''stress11: '',5f12.3)') stress11tot,sep_stress11tot,dstress11tot,divtot,presstot
!write(irefwr,*) 'dimensional pressure average: ',presstot*pressure_scale
visdiptot = visdiptot/npoint
sep_visdiptot = sep_visdiptot/npoint
dvisdip=dvisdip/npoint
write(irefwr,'(''visdipt : '',5f12.3)') visdiptot,sep_visdiptot,dvisdip,dvisdip/visdiptot*100
write(irefwr,*) 'visdipmax: ',visdipmax(1:2)
write(irefwr,*) 's11pmax  : ',s11max(1:2)

end  subroutine checkstress002

real(kind=8) function getrmsdiv(kmesh,kprob,isol,div,iuser,user)
use sepmodulevecs
use sepmodulekmesh
use sepmodulekprob
use geometry
use coeff
implicit none
integer :: iuser(*),kmesh,kprob,isol
real(kind=8) :: div(*),user(*)
real(kind=8) :: volint,x,y,rho,funccf
integer :: ihelp,i

write(irefwr,*) 'user: ',iuser(1),user(1)
if (eos_type /= 0) then
   write(irefwr,*) 'PERROR(getrmsdiv): not yet suited for eos_type /= 0'
   call instop
endif

write(irefwr,*) 'user: ',user(1),iuser(1)
iuser(2)=1
iuser(3:5)=0
iuser(6)=7
iuser(10)=2001
iuser(11)=6
call sepactsolbf1(isol)

call getrmsdiv01(npoint,ndim,coor,ks(isol)%sol,div,nunkpi,nphysi,indprfi,kprobfi,indprpi,kprobpi,user(6))
getrmsdiv=sqrt(volint(0,1,1,kmesh,kprob,isol,iuser,user,ihelp))/volume

end function getrmsdiv

subroutine getrmsdiv01(npoint,ndim,coor,usol,div,nunkp,nphys,indprf,kprobf,indprp,kprobp,user)
implicit none
integer :: npoint,ndim,nunkp,indprf,kprobf(*)
integer :: indprp,nphys,kprobp(npoint,nphys)
real(kind=8) :: usol(*),div(*),user(*),coor(ndim,npoint),funccf,functional,drhodz,rho,w
real(kind=8) :: x,y
integer :: j1,j2,i

!write(irefwr,*) 'indprf/p: ',indprf,indprp,nunkp,nphys,ndim
do i=1,npoint
   x=coor(1,i)
   y=coor(2,i)
   if (indprf == 0 .and. indprp == 0) then
      j1 = (i-1)*nunkp+1
      j2 = (i-1)*nunkp+2
   else if (indprp /= 0) then
      j1 = kprobp(i,1)
      j2 = kprobp(i,2)
   else
      j1 = kprobf(i) + 1
      j2 = j1 + 1
   endif
   w=usol(j2)
   rho=funccf(3,x,y,y)
   drhodz=funccf(13,x,y,y)
   ! from chain rule on div(rho u)
   functional=rho*div(i) - w*drhodz
   !if ((i/100)*100==i) then
   !  write(irefwr,'(i5,6f12.3)') i,x,y,rho,drhodz,w,div(i)
   !endif
   user(i) = functional*functional
enddo

end subroutine getrmsdiv01
