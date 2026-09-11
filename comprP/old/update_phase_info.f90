!     UPDATE_PHASE_INFO
!     Update phase change info to account for adiabat/reference temperature
!     For EBA the phase boundary reference temperature is set to zero
!     For ALA it is set to that of the adiabat
!     PvK October 2011
subroutine update_phase_info(iuser,user)
use sepmodulecomio
use convparam
use coeff
use mparallel
use geometry
use convparam
use control
implicit none
real(kind=8) :: rho_max,rho1,rho2,funccf,z,rl,height,dz,T_bot_dim
real(kind=8) :: deltaT_dim_old,adiabat,deltaT_dim_eos,Tso_dim_old
real(kind=8) :: avHl(4),user(*)
integer :: lu,find_eos_iz,iz,i,iph,iuser(*)

if (EBA.and.(hannah.or.sky)) then
   if (print_node) then
      if (nph>0) write(irefwr,'('' Reset phase change T to 0 '')')
   endif
   pht0_d(1:nph) = T0_nondim
   pht0(1:nph) = 0d0
   return
else if (EBA.or.BA) then
   if (delta1K) then
      pht0(1:nph) = pht0_d(1:nph)
   else
      pht0(1:nph) = pht0_d(1:nph)/deltaT_dim
   endif
endif

if (compress) then ! .and.(.not.solve_for_Tperturb)) then
!  Reset pht0 to ambient T, which is that of the adiabat
   if (nph > 0 .and. print_node) then
      write(irefwr,*)
      write(irefwr,'('' Reset phase change T to adiabatic T '')')
   endif
   do i=1,nph
!     Find corresponding entry in look up table
!     z is depth in km; pht0_d is in K; pht0 is in C equivalent
      z = phz0_d(i)*1e-3
      iz = find_eos_iz(z,eos_z_d,eos_np)
      dz = eos_z_d(iz+1)-eos_z_d(iz)
      rl = (z-eos_z_d(iz))/dz
      pht0_d(i) = eos_T_d(iz)*(1d0-rl)+eos_T_d(iz+1)*rl  ! in K
      if (delta1K) then
         pht0(i) = (pht0_d(i)-T0_nondim)
      else
         pht0(i) = (pht0_d(i)-T0_nondim)/deltaT_dim
      endif
   enddo
endif ! compress .and. (.not.solve_for_Tperturb)

if (print_node) then
      write(irefwr,*) 'Updated phase change info, nph=',nph
      do i=1,nph
         write(irefwr,'('' phase change '',i5,'':  phT0_d (in K) = '',f12.1,''; phT0 = '',4f12.3)') i,pht0_d(i),pht0(i), &
               &    phz0(i),R2-phz0(i)
      enddo
endif

end subroutine update_phase_info
