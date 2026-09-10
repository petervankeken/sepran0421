! PvK July 2024
! COEFVIS
! Set up ivisc() with updated viscosity
subroutine coefvis
use sepmodulecomio
use sepmodulevecs
use sepmodulekmesh
use sepmodulecomio
use sepran_arrays
use sepmoduleoldrouts
use coeff
use control
implicit none
interface
   subroutine coefv2(ndim,npoint,coor,temp,adia,eta,etalin,etaplast,secinv,rho)
   integer,intent(in) :: ndim,npoint
   real(kind=8),intent(in) :: coor(ndim,npoint),temp(npoint),adia(npoint),secinv(npoint),rho(npoint)
   real(kind=8),intent(out) :: eta(npoint),etalin(npoint),etaplast(npoint)
   end subroutine coefv2
end interface
integer :: iuserlc(100),inpcre(10),ix,ielhlp=0,jdegfd=0
real(kind=8) :: userlc(100),rincre(10)

iuserlc=0
 userlc=0.0_8
inpcre=0
rincre=0
iuserlc(1)=100
iuserlc(2)=1 ! signal there is only one element group
 userlc(1)=100*1.0_8

if (ivisc==0) then
   if (print_node) write(irefwr,*) 'PERROR(coefvis): ivisc=0'
   call instop
endif

if (ivn.or.tackley.or.tosi15) then
   ! find strain rate in nodal points
   iuserlc(1)=100
   iuserlc(2)=1
   iuserlc(3:5)=0
   iuserlc(6)=7
   iuserlc(7)=0 ! type of NS equation
   iuserlc(8)=1 ! type of constitutive equation
   iuserlc(9)=0 ! type of interpolation
   iuserlc(10)=icoor900 ! type of coordinate system
   iuserlc(11)=mcontv ! compressibility formulation
   iuserlc(12)=-7 ! eps
   iuserlc(13)=-8 ! rho
   iuserlc(18)=-9 ! eta
   iuserlc(15:100)=0
    userlc(2:5)=0.0_8
    if (itype_stokes==903) then
       userlc(6)=0.0_8
    else
       userlc(6)=penalty_parameter
    endif
    userlc(8)=1.0_8
    userlc(9)=1.0_8
   call deriva(2,10,ix,jdegfd,1,isecinv,kmesh,kprob,isol(1),isol(1),iuserlc,userlc,ielhlp)
   pedebug=.false.
   if (pedebug) then
     write(irefwr,*) 'secinv sum: ',anorm(1,1,1,kmesh,kprob,isecinv,isecinv,ielhlp)
     write(irefwr,*) 'secinv max: ',anorm(1,3,1,kmesh,kprob,isecinv,isecinv,ielhlp)
     write(irefwr,*) 'uvel   max: ',anorm(1,3,1,kmesh,kprob,isol(1),isol(1),ielhlp)
     write(irefwr,*) 'wvel   max: ',anorm(1,3,2,kmesh,kprob,isol(1),isol(1),ielhlp)
   endif
   pedebug=.false.
endif

call sepactsolbf1(isol(2))
!write(6,*) 'before coefv2: ',iadia,ivisc,isecinv,ivisclin,iviscplast,idens
!write(6,*) 'ndim: ',ndim
call coefv2(ndim,npoint,coor,ks(isol(2))%sol,ks(iadia)%sol,ks(ivisc)%sol,ks(ivisclin)%sol,ks(iviscplast)%sol, &
          & ks(isecinv)%sol,ks(idens)%sol)



end subroutine coefvis

subroutine coefv2(ndim,npoint,coor,temp,adia,eta,etalin,etaplast,secinv,rho)
use sepmodulecomio
use geometry
use control
use coeff
implicit none
integer,intent(in) :: ndim,npoint
real(kind=8),intent(in) :: coor(ndim,npoint),temp(npoint),adia(npoint),secinv(npoint),rho(npoint)
real(kind=8),intent(out) :: eta(npoint),etalin(npoint),etaplast(npoint)
real(kind=8) :: temp_a,temp_l,secsqr,pefvis,x,y
integer :: i,ival,ivalfind,ipl,ipp
real(kind=8) :: vislmin,vislmax,vispmin,vispmax

!write(6,*) 'coefv2, tmin/tmax : ',minval(temp(1:npoint)),maxval(temp(1:npoint))

vislmin=1e9
vislmax=-1e9
vispmin=1e9
vispmax=-1e9
ipl=0
ipp=0
eta=0.0_8
etaplast=0.0_8
etalin=0.0_8

!write(irefwr,*) 'ivl: ',ivl,nlay,r_zint(1:nlay)
if (itypv==0) then
   eta(1:npoint)=1.0_8
else if (.not.tackley.and..not.tosi15) then
   ival=1
   secsqr=0
   do i=1,npoint
      x = coor(1,i)
      y = coor(2,i)
      if (ivl) ival = ivalfind(x,y)
      temp_l=temp(i)
      temp_a=adia(i)
      secsqr=secinv(i)*secinv(i)
      eta(i) = pefvis(x,y,temp_l,temp_a,ival,secsqr)
   enddo
else
      ! tackley or tosi
   ival=1
   secsqr=0
   do i=1,npoint
      x = coor(1,i)
      y = coor(2,i)
      if (ivl) ival = ivalfind(x,y)
      temp_l=temp(i)
      temp_a=adia(i)
      secsqr=secinv(i)*secinv(i)
      eta(i) = pefvis(x,y,temp_l,temp_a,ival,secsqr)
      etalin(i) = etalinh
      etaplast(i) = etaplasth
   enddo
endif

etamin=minval(eta(1:npoint))
etamax=maxval(eta(1:npoint))
!write(6,'(''PINFO(coefv2): etamin/max = '',2f12.4)') etamin,etamax
!call instop

end subroutine coefv2
