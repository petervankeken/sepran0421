subroutine averageT(ichoice,iuser,user,avT)
use sepmodulecomio
use sepmodulekmesh
use sepmodulevecs
use sepran_arrays
use coeff
use convparam
use geometry
use control
implicit none
real(kind=8) :: user(*),avT(*)
integer :: iuser(*)
real(kind=8) :: volint,rinvec(2)
integer :: ihelp,i,ichoice,iinvec(2),iinvol(5)
!     include 'pecof800.inc'
!     include 'ccc.inc'
!     include 'dimensional.inc'
!     include 'compression.inc'

iuser(2:iuser(1))=0
user(7:int(user(1)))=0.0_8
! input for volint
iuser(2)=1
iuser(6)=7
iuser(7)=intrule800+100*interpol800
iuser(8)=icoorsystem
iuser(9)=0
! specify f
iuser(10)=-6
 user(6)=1d0
!call sepactsolbf1(isol2)

! input for integr
iinvol(1)=5
iinvol(2)=1
iinvol(3)=1
iinvol(4)=0
iinvol(5)=0

if (isol(3)==0) then
   if (print_node) then
      write(irefwr,*) 'PERROR(averageT):: isol(3) has not been filled'
      write(irefwr,*) 'Make sure to use compute_total_T'
   endif
   call instop
endif
   

ihelp=0
iinvec(1)=2
iinvec(2)=3
if (ichoice==1) then
!  write(6,*) 'averageT: ',ichoice
   avT(1) = volint(0,2,1,kmesh,kprob,isol(3),iuser,user,ihelp)
   !call integr(iinvol,avT(1),kmesh,kprob,isol(3),iuser,user)
   call manvec(iinvec,rinvec,isol(3),isol(3),isol(3),kmesh,kprob)
!  write(6,*) 'averageT temp = ',avT(1)
else if (ichoice==2) then
!  write(6,*) 'averageT: ',ichoice
   ! for some reason we need to copy iphi from problem 1 to a vector with problem 2
   ! Average viscous dissipation as computed by sepran from eta * sqrt (II)
!  write(6,*) 'iphi->iphi2'
   call change_problem_for_phi(iphi,iphi2)
   avT(1) = volint(0,2,1,kmesh,kprob,iphi2,iuser,user,ihelp)
!   write(6,*) 'iphi2: '
!  call print_vector(npoint,ks(iphi2)%sol)
   !call integr(iinvol,avT(1),kmesh,kprob,iphi2,iuser,user)
   call manvec(iinvec,rinvec,iphi2,iphi2,iphi2,kmesh,kprob)
!  write(6,*) 'averageT iphi2 = ',avT(1)
else if (ichoice==3) then
!  write(6,*) 'averageT: ',ichoice
!  write(6,*) 'icompwork: ',icompwork
!  call print_vector(npoint,ks(icompwork)%sol)
   avT(1) = volint(0,2,1,kmesh,kprob,icompwork,iuser,user,ihelp)
   !call integr(iinvol,avT(1),kmesh,kprob,icompwork,iuser,user)
   call manvec(iinvec,rinvec,icompwork,icompwork,icompwork,kmesh,kprob)
!  write(6,*) 'averageT work = ',avT(1)
else if (ichoice==7) then
!  write(6,*) 'averageT: ',ichoice
   ! Average viscous dissipation as computed by PvK from grad v
!  write(6,*) 'ivisdip: ',ivisdip,ivisdip2
!  write(6,*) 'ivisdip->ivisdip2'
!  write(6,*) 'visdip: '
!  call print_vector(npoint,ks(ivisdip)%sol)
   call change_problem_for_phi(ivisdip,ivisdip2)
!  write(6,*) 'visdip2: '
!  call print_vector(npoint,ks(ivisdip2)%sol)
!  write(6,*) 'ichoice==7'
   call copyvc(ivisdip,isol(3))
   avT(1) = volint(0,2,1,kmesh,kprob,isol(3),iuser,user,ihelp)
!  write(6,*) 'average visdip after copy to isol3:',avT(1)
   avT(1) = volint(0,2,1,kmesh,kprob,ivisdip,iuser,user,ihelp)
   call manvec(iinvec,rinvec,ivisdip2,ivisdip2,ivisdip2,kmesh,kprob)
!  write(6,*) 'averageT visdip: ',avT(1)
else if (ichoice==10) then
!  write(6,*) 'averageT: ',ichoice
   ! average conductivity
   !avT(1)=volint(0,2,1,kmesh,kprob,icond,iuser,user,ihelp)
   call integr(iinvol,avT(1),kmesh,kprob,icond,iuser,user)
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
else if (ichoice == 2 .or. ichoice == 7) then
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
interface
  subroutine change_problem_for_phi01(nusol1,nusol2,npoint,usol1,usol2)
  integer,intent(in) :: nusol1,nusol2,npoint
  real(kind=8),intent(in) :: usol1(*)
  real(kind=8),intent(out) :: usol2(*)
  end subroutine change_problem_for_phi01
end interface
integer :: iphi,iphi2

call change_problem_for_phi01(ks(iphi)%nusol,ks(iphi2)%nusol,npoint,ks(iphi)%sol,ks(iphi2)%sol)


end  subroutine change_problem_for_phi

subroutine change_problem_for_phi01(nusol1,nusol2,npoint,usol1,usol2)
use sepmodulecomio
use control
implicit none
integer,intent(in) :: nusol1,nusol2,npoint
real(kind=8),intent(in) :: usol1(*)
real(kind=8),intent(out) :: usol2(*)
integer :: i

if (nusol1 /= nusol2) then
   if (print_node) then 
      write(irefwr,*) 'PERROR(change_problem_for_phi): nusol1 <> nusol2'
      write(irefwr,*) 'nusol1, nusol2 = ',nusol1,nusol2
   endif
   call instop
endif
do i=1,nusol2
   usol2(i)=usol1(i)
enddo

end subroutine change_problem_for_phi01

subroutine print_vector(npoint,sepranvector)
use sepmodulecomio
use control
implicit none
integer,intent(in) :: npoint
real(kind=8),intent(in) :: sepranvector(*)
integer :: i

if (print_node) then
  do i=1,npoint,npoint/10
     write(irefwr,*) sepranvector(i)
  enddo
endif

end subroutine print_vector
