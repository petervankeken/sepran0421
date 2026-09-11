! COEFVIS
! Make array that contains viscosity in the nodal points
! PvK 070819 adapted for 0619
subroutine coefvis
use sepmodulecomio
use sepmodulevecs
use sepmodulekmesh
use sepran_arrays
use coeff
use control
implicit none
interface
   subroutine coefv2(ndim,npoint,coor,temp,adia,eta,secinv,rho)
   integer,intent(in) :: ndim,npoint
   real(kind=8),intent(in) :: coor(ndim,npoint),temp(npoint),adia(npoint),secinv(npoint),rho(npoint)
   real(kind=8),intent(out) :: eta(npoint) 
   end subroutine coefv2
end interface
integer :: iuserlc(100),inpcre(10)
real(kind=8) :: userlc(100),rincre(10)
integer :: ix=0,ielhlp=0,jdegfd=0

!if (print_node) pedebug=.true.

if (pedebug) write(irefwr,*) 'in coefvis: ',ivisc

iuserlc=0
userlc=0
iuserlc(1)=100
iuserlc(2)=1
 userlc(1)=100.0_8
if (ivisc==0) then
   inpcre(1)=6
   inpcre(2)=2
   inpcre(3)=1
   inpcre(4)=1
   inpcre(5)=0
   inpcre(6)=0
   call creatv(kmesh1,kprob1,ivisc,inpcre,rincre)
   if (pedebug) write(irefwr,*) 'created ivisc'
endif
if (ivn.or.tackley) then
   call deriva(2,10,ix,jdegfd,1,isecinv,kmesh1,kprob1,isol1,isol1,iuser_here,user_here,ielhlp)
else 
   isecinv=1 ! needs to be something different from 0..
endif

if (idens==0) idens=1

!if (pedebug) write(irefwr,*) 'isol2: ',isol2
call sepactsolbf1(isol2)
call sepgetmeshinfo(ndim,npoint)
if (pedebug) write(irefwr,*) 'isol2: ',isol2,iadia,ivisc,isecinv,idens
call coefv2(ndim,npoint,coor,ks(isol2)%sol,ks(iadia)%sol,ks(ivisc)%sol,ks(isecinv)%sol,ks(idens)%sol)

end subroutine coefvis

subroutine coefv2(ndim,npoint,coor,temp,adia,eta,secinv,rho)
use sepmodulecomio
use geometry
use control
use coeff
implicit none
integer,intent(in) :: ndim,npoint
real(kind=8),intent(in) :: coor(ndim,npoint),temp(npoint),adia(npoint),secinv(npoint),rho(npoint)
real(kind=8),intent(out) :: eta(npoint) 
real(kind=8) :: etamin,etamax,temp_a,temp_l,secsqr,pefvis,x,y
integer :: i,ival,ivalfind

!if (print_node) pedebug=.true.
etamin=5e9
etamax=0


!write(irefwr,*) 'ivl: ',ivl,nlay,r_zint(1:nlay)
if (itypv==0) then
   eta(1:npoint)=1.0_8
else
   ival=1
   secsqr=0
   do i=1,npoint
      x = coor(1,i)
      y = coor(2,i)
!     write(irefwr,*) 'x,y: ',x,y
      if (ivl) ival = ivalfind(x,y)
!     write(irefwr,*) 'ival: ',ival
      temp_l=temp(i)
!     write(irefwr,*) 'temp_l',temp_l
      temp_a=adia(i)
!     write(irefwr,*) 'temp_a ',temp_a
      if (ivn.or.tackley) secsqr=secinv(i)*secinv(i)
      !if (pedebug) write(irefwr,'(''before pefvis: '',5e15.7)') x,y,temp_l,temp_a,secsqr
      eta(i) = pefvis(x,y,temp_l,temp_a,ival,secsqr)
   enddo
endif
!if (pedebug) write(irefwr,*) 'min/max viscosity = ',minval(eta(1:npoint)),maxval(eta(1:npoint))
      
end subroutine coefv2

