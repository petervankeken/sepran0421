!   PEFILCOF
! 
!   Generic routine for specifying coefficients for the Stokes
!   and heat equation in mantle convection applications.
!   
!   Updated to conform with element 903
!   PvK 071705
subroutine pefilcof(ichoice,iuser,user)
use coeff
use geometry
use control
implicit none
integer :: ichoice,iuser(*)
real(kind=8) :: user(*)
 
if (cyl) then
   if (ichoice==2.or.stokes_via_user) then
      call pefilcof_cart(ichoice,iuser,user)
   else
      call pefilcof_cyl_old(ichoice,iuser,user)
   endif
else 
   call pefilcof_cart(ichoice,iuser,user)
endif
    
end subroutine pefilcof

!   PEFILCOF_CART
!
!   Modified for ichoice=1 Stokes equation. Use modelv=103
!
!   PvK 981108/040900
subroutine pefilcof_cart(ichoice,iuser,user)
use sepran_arrays
use sepran_interface
use coeff
use brandenburg
use control
implicit none
integer :: ichoice
integer :: iuser(*)
real(kind=8) :: user(*)
real(kind=8) :: work(100),f2

integer :: iwork(4,100),itime,modelv,npoint,iincop(6),i
integer :: iinpri(5),icurvs(5),iprhocp,k,nparm,j,inputcr(10)
real(kind=8) :: rinputcr(10),format,d2adiabatdz2,adiabat
character(len=80) :: text
integer :: ishape1,kelmf1,iu1(3)
real(kind=8) :: u1(3)

real(kind=8), parameter :: DONE=1.0_8

if (ichoice==1 .or. ichoice==4) then
!     PvK 031910 rewritten to fill iuser/user directly which is            
!     easier when coefficients that are defined in the nodal points
!     (which is the case with eos_type /= 0)
!     irho, ialpha, icp, ibulkmod contain the reference values
!     ipress, ivisc contain the current values for pressure and viscosity
!
!     irho should have been filled in compr_start
!     if eos_type /=0 then ialpha, icp, ibulkmod should have been filled in compr_start
!     for compress press needs to be filled too
!     idens contains the chemical buoyancy (if used)

!     fill ivisc when abs(itypv)>0
     if (itypv/=0.and..not.gable_plates) then
        call coefvis()
     else if (itypv/=0.and.(gable_plates.and.gable_stokes_choice==0)) then
        call coefvis()
     endif
 
     if (compress.and..not.tala.and.itype_stokes==903) then
        ! find pressure by extracting it from the solution vector
        ! for element 900 we need to first fill coefficients and then use deriv
        call coefpress(kmesh1,kprob1,isol1,ipress,iuser,user)
     endif
     call fillcoef900(ichoice,iuser,user)
     if (.not.tala.and.compress.and.itype_stokes==900) then
        ! now that coefficients are filled for either TALA or out of date pressure we can find pressure from deriv and update
        ! coefficients
        if (print_node) then
        endif
        !all instop
        call coefpress900()
        ! and update coefficients 
        call fillcoef900(ichoice,iuser,user)
     else if (tala.and.compress.and.itype_stokes==900) then
        call setpresstozero()
     endif
      

else if (ichoice==2) then


   call fillcoef800(iuser,user)

else if (ichoice==3) then

   do i=6,100
      iuser(i)=0
   enddo
!  Element group 1
   iuser(6)=10
!  Boundary element group 1
!  iuser(7)=30

!  itime, modelv, numint, icoor, mcont
   iuser(10)=0
   iuser(11)=1
   iuser(12)=intrule800+100*interpol800
   iuser(13)=1
   iuser(14)=0
!  eps rho omega f1 f2 f3 eta 
   iuser(15)=-6
   iuser(16)=-7
   iuser(17)=0
   iuser(18)=0
   iuser(19)=0
   iuser(20)=0
   iuser(21)=-7
    user(6)=1d-6
    user(7)=1d0

endif

end subroutine pefilcof_cart

subroutine coefpress900()
use sepmodulecomio
use sepmodulecpack
use sepmodulekmesh
use sepmodulevecs
use sepran_arrays
use sepran_interface
use control
implicit none
integer ::  iinder(8),iu1(3)


iinder(1)=8
iinder(2)=1
iinder(3)=0
iinder(4)=20  ! pressure in the elements
iinder(5)=0
iinder(6)=0
iinder(7)=0
iinder(8)=2
! ipresse is special vector containing value of pressure in the elements (centroids)
call deriv(iinder,ipresse,kmesh1,kprob1,isol1,iuser_here,user_here)
! shift so that P is all positive (to match 903 with essential b.c.)
call print_and_shift_P(ks(ipresse)%sol,nelem)
if (ipress==0) then
   iu1=0
    u1=0
   call creavc(0,1001,1,ipress,kmesh1,kprob1,iu1,u1,iu1,u1)
endif
call convert_elem_to_coor(npoint,nelem,ndim,npelm,coor,kmeshc,ks(ipresse)%sol,ks(ipress)%sol)

end subroutine coefpress900

subroutine setpresstozero()
use sepmodulekmesh
use sepmodulevecs
use sepran_arrays
implicit none  

call setpressreallytozero(ks(ipress)%sol,npoint)

end subroutine setpresstozero

subroutine setpressreallytozero(press,npoint)
implicit none
integer :: npoint
real(kind=8) :: press(*)
press(1:npoint)=0.0_8
end subroutine setpressreallytozero
  
