! coefpress
!
! create pressure vector for use in specification of ALA buoyancy term
! PvK 17 Aug 2020
subroutine coefpress()
use sepmodulecomio
use sepmodulekmesh
use sepmodulekprob
use sepmodulevecs
use control
use coeff
use sepran_arrays
implicit none
integer :: iinder(10),ipoint,ihelp
real(kind=8) :: pee,quu,anorm,maxpressure

if (ipress==0) then
   iu1lc(1)=0
    u1lc(1)=0.0_8
   call creavc(0,1001,1,ipress,kmesh1,kprob1,iu1lc,u1lc,iu1lc,u1lc)
endif
call sepactsolbf1(isol1)

if (ishape/=4) then
   if (print_node) write(irefwr,*) 'PERROR(coefpress): may need adjustments for ishape other than 4'
   call instop
endif

if (itype_stokes==903) then
   call coefpress01(kmeshc,ks(ipress)%sol,ks(isol1)%sol,indprpi,kprobpi,npoint,nphysi,nelem)
else 
   if (print_node) write(irefwr,*) 'PERROR(coefpress): until further notice only suitable for TH elements'
   call instop
endif

! some test code that doesn't lead to decent results
!pvk else if (itype_stokes==900) then
!pvk    iinder(1)=8
!pvk    iinder(2)=1
!pvk    iinder(3)=0
!pvk    iinder(4)=4  ! (pressure per ! 4
!pvk    iinder(5)=0
!pvk    iinder(6)=0
!pvk    iinder(7)=0
!pvk    iinder(8)=2
!pvk    call deriv(iinder,idivv,kmesh1,kprob1,isol1,iuser_here,user_here)
!pvk    ! NB fixed assumed penalty function parameter
!pvk    call algebr(3,1,idivv,idivv,ipress,kmesh1,kprob1,1.0e3_8,0.0_8,pee,quu,ipoint)
!pvk    !iinder(4)=20
!pvk    !iuser_here(6) =8
!pvk    !iuser_here(8) =0
!pvk    !iuser_here(9) =1
!pvk    !iuser_here(10)=0
!pvk    !iuser_here(11)=0
!pvk    !iuser_here(12)=mcontv
!pvk    !iuser_here(13)=-7
!pvk    ! user_here(7)=1.0e-6_8
!pvk    !write(irefwr,*) 'before finding pressure for 900'
!pvk    !call deriv(iinder,ipress,kmesh1,kprob1,isol1,iuser_here,user_here)
!pvk endif
!pvk maxpressure= anorm(1,1,1,kmesh1,kprob1,ipress,ipress,ihelp)
!pvk write(irefwr,*) 'max pressure:',maxpressure
!pvk if (idivv>0) then
!pvk   maxpressure= anorm(1,1,1,kmesh1,kprob1,idivv,idivv,ihelp)
!pvk   write(irefwr,*) 'max div u:',maxpressure
!pvk endif

end subroutine coefpress

subroutine coefpress01(kmeshc,press,uvp,indprp,kprobp,npoint,nphys,nelem)
use sepmodulecomio
use control
implicit none
integer :: nelem,npoint,nphys,indprp
integer :: kmeshc(*),kprobp(npoint,nphys)
real(kind=8) :: press(*),uvp(*),pres_local(6)
integer :: ip,ielem,j,j1,nodno

!write(irefwr,*) 'nelem: ',nelem,npoint,nphys,indprp
if (indprp==0) then
   ! not a quadratic Taylor-Hood element
   if (print_node) write(irefwr,*) 'PERROR(coefpress01): indprp==0'
   call instop
endif
if (nphys /= 3) then
   ! Not a quadratic Taylor Hood element
   if (print_node) write(irefwr,*) 'PERROR(coefpress01): nphys /= 3'
   if (print_node) write(irefwr,*) 'nphys = ',nphys
   call instop
endif
ip=0
do ielem=1,nelem ! loop over elements
   do j=1,6 ! loop over nodal points in element
      nodno=kmeshc(ip+j)
      j1=kprobp(nodno,3)
      if (j1 /= 0) then
         ! got a point with pressure
         pres_local(j)=uvp(j1)
      else
         pres_local(j)=0.0_8
      endif
   enddo
   ! fill the rest of the points by averaging
   ! Fill rest of points by simple averaging
   pres_local(2)=0.5*(pres_local(1)+pres_local(3))
   pres_local(4)=0.5*(pres_local(3)+pres_local(5))
   pres_local(6)=0.5*(pres_local(5)+pres_local(1))
   ! Store pressure in global vector
   do j=1,6
      nodno=kmeshc(ip+j)
      press(nodno) = pres_local(j)
   enddo
   ip=ip+6
enddo ! ielem=1,nelem
!write(irefwr,*) 'minimum/maximum pressure: ',minval(press(1:npoint)),maxval(press(1:npoint))

end subroutine coefpress01
