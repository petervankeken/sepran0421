subroutine averageT(ichoice,iuser,user,avT)
use sepmodulecomio
use sepran_arrays
use coeff
use convparam
use geometry
use control
implicit none
real(kind=8) :: user(:),avT(:)
integer :: iuser(:)
real(kind=8) :: volint,rinvec(2)
integer :: ihelp,i,ichoice,iinvec(2)
!     include 'pecof800.inc'
!     include 'ccc.inc'
!     include 'dimensional.inc'
!     include 'compression.inc'

iuser(2:iuser(1))=0
user(7:int(user(1)))=0.0_8
iuser(2)=1
iuser(6)=7
iuser(7)=intrule800+100*interpol800
iuser(8)=icoorsystem
iuser(9)=0
! specify f
iuser(10)=-6
 user(6)=1d0
!call sepactsolbf1(isol2)

if (isol3==0) then
   if (print_node) then
      write(irefwr,*) 'PERROR(averageT):: isol3 has not been filled'
      write(irefwr,*) 'Make sure to use compute_total_T'
   endif
   call instop
endif
   

iinvec(1)=2
iinvec(2)=3
if (ichoice==2) then
   ! for some reason we need to copy iphi from problem 1 to a vector with problem 2
   call change_problem_for_phi(iphi,iphi2)
   avT(1) = volint(0,2,1,kmesh1,kprob1,iphi2,iuser,user,ihelp)
   call manvec(iinvec,rinvec,iphi2,iphi2,iphi2,kmesh1,kprob1)
else  if (ichoice==1) then
   avT(1) = volint(0,2,1,kmesh1,kprob1,isol3,iuser,user,ihelp)
   call manvec(iinvec,rinvec,isol3,isol3,isol3,kmesh1,kprob1)
else if (ichoice==3) then
   avT(1) = volint(0,2,1,kmesh1,kprob1,icompwork,iuser,user,ihelp)
   call manvec(iinvec,rinvec,icompwork,icompwork,icompwork,kmesh1,kprob1)
endif

avT(3)=rinvec(2) ! maximum value
avT(5)=rinvec(1) ! minimum value


avT(1) = avT(1)/volume
if (ichoice == 1) then
   if (ibench_type>0.and.Kequivalent) then
      ! subtract surface temperature 
      avT(1)=avT(1)-273.0 ! patch for CY85 Ts_nondimK
   endif
  ! temperature
  avT(2) = avT(1)*deltaT_dim
else if (ichoice == 2) then
  ! viscous dissipation (unit: work = J/s)
  ! This needs work
  !avT(2) = rho_dim*avT(1)/height_dim
  !avT(2) = avT(2)*rkappa_dim*rkappa_dim*rkappa_dim
  if (delta1K) avT(1)=avT(1)/deltaT_dim
else if (ichoice == 3) then
  ! vertical advective heat transport
  ! This needs work
  ! avT(2) = rho_dim*cp_dim*DeltaT_dim*avT(1)
  ! avT(2) = avT(2)*height_dim*height_dim
  if (delta1K) avT(1)=avT(1)/deltaT_dim
else if (ichoice == 4) then
  avT(2) = avT(1)*rho_dim
else if (ichoice == 5) then
  avT(2) = avT(1)*rkappa_dim
else if (ichoice == 6) then
  avT(2) = avT(1)*alpha_dim
endif


end subroutine averageT

subroutine change_problem_for_phi(iphi,iphi2)
use sepmodulekmesh
use sepmodulesol
implicit none
integer :: iphi,iphi2

call change_problem_for_phi01(ks(iphi)%nusol,ks(iphi2)%nusol,npoint,ks(iphi)%sol,ks(iphi2)%sol)


end  subroutine change_problem_for_phi

subroutine change_problem_for_phi01(nusol1,nusol2,npoint,usol1,usol2)
use sepmodulecomio
use control
implicit none
integer :: nusol1,nusol2,npoint
real(kind=8) :: usol1(*),usol2(*)
integer :: i

if (nusol1 /= nusol2) then
   if (print_node) write(irefwr,*) 'PERROR(change_problem_for_phi): nusol1 <> nusol2'
   call instop
endif
do i=1,nusol2
   usol2(i)=usol1(i)
enddo

end subroutine change_problem_for_phi01
