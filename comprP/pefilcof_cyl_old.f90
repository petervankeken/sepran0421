subroutine pefilcof_cyl_old(ichoice,iuser,user)
use sepmodulecomio
use sepran_arrays
use sepmodulekmesh
use sepmodulevecs
use coeff
use control
implicit none
integer :: iuser(*),ichoice
real(kind=8) :: user(*),work(100),f2

integer :: iwork(4,100),jtime,modelv,iincop(6),i
integer :: iinpri(5),icurvs(5),iprhocp,k,nparm,j,inputcr(10)
real(kind=8) :: rinputcr(10),format
character(len=80) :: text

real(kind=8) :: DONE
parameter(DONE=1.0_8)

if (ichoice == 1 .or. ichoice==4) then
   ! clear work/iwork
   do i=1,100
      do k=1,4
        iwork(k,i) = 0
      enddo
      work(i) = 0d0
   enddo
   do i=1,iuser(1)-6
      iuser(5+i)=0
   enddo


   !  Coefficients for element 900
   !  First five, integer :: info
   modelv=  103
   jtime = 0
   iwork(1,1) = jtime
   iwork(1,2) = modelv

   iwork(1,3) = intrule900 + 100*interpol900
   iwork(1,4) = icoor900
   iwork(1,5) = mcontv

   if (itype_stokes.eq.900) then
      ! penaltyfunction parameter
      eps = 1e-6_8
   else
      eps = 0
   endif
   f2     = 1.0_8

   ! 6: epsilon
   iwork(1,6) = 0
    work(6) = eps
   !   7: rho
   if (compress) then
      ! Density is variable: specify with funccf(3)
      if (eos_type == 0) then
         iwork(1,7) = 3
      else
         ! density is stored in an old vector
         iwork(1,7) = 2003
         ! use old vector number 2
         iwork(2,7) = 4
         ! degree of freedom 1
         iwork(3,7) = 1
      endif
   else 
      iwork(1,7) = 0
       work(7)   = 1.0_8
   endif
!  8: omega
!  9: f1. First dof of third vector in islold
!  This contains only the thermal (+phase change) components
!  The chemical buoyancy is added in elm900_pvk.
   if (ichoice==1) then
      iwork(1,9) = 2003
      iwork(2,9) = 3
      iwork(3,9) = 1
   endif
!  10: f2. Second dof of third vector in islold.
!  This contains only the thermal (+phase change) components
!  The chemical buoyancy is added in elm900_pvk.
   if (ichoice==1) then
      iwork(1,10) = 2003
      iwork(2,10) = 3
      iwork(3,10) = 2
   endif
!  11: f3
   iwork(1,11) = 0
    work(11) = 0.0_8
!  12: eta 
!  constant viscosity OR specified through fnv003
   iwork(1,12) = 0
    work(12)   = 1.0_8
   call cofcopy(0,isolold1(2:2),isol2,kmesh1,kprob1,DONE)
   call f12copy(isolold1(3:3),isol2,isol1,ipress,kmesh1,kprob1)
   nparm=12
   call fil103(1,1,iuser,user,kprob1,nparm,iwork,work,kmesh1)
else
   if (print_node) write(irefwr,*) 'PERROR(pefilcof_cyl_old): call with ichoice=',ichoice
endif

if (print_node) then
   write(irefwr,*) 'iwork / work: '
   do i=1,12
       write(irefwr,'(3i10,e15.7)') (iwork(j,i),j=1,3),work(i)
   enddo
   call print_coeffs_in_isolold(ks(isolold1(2))%sol,ks(isolold1(3))%sol,npoint)
endif

end subroutine pefilcof_cyl_old

subroutine print_coeffs_in_isolold(uold1,uold2,npoint)
use sepmodulecomio
use control
real(kind=8) :: uold1(*),uold2(*)
integer :: npoint
integer :: i,ip

if (print_node) then
  do i=1,npoint,npoint/10
     ip=2*i
     write(irefwr,'(4e15.7)') uold1(ip-1),uold1(ip),uold2(ip-1),uold2(ip)
  enddo
endif

end subroutine print_coeffs_in_isolold
