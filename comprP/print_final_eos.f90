subroutine print_final_eos()
use sepmodulecomio
use coeff
use convparam
use mparallel
use control
implicit none
integer :: iph,lu,i
character(len=120) :: fname

if (print_node) then
  open(11,file='eos_T.dat')
  open(12,file='eos_rho.dat')
  do iph=1,nph
    lu=12+iph
    write(fname,'(''eos_Gamma'',i2.2,''.dat'')') iph
    open(lu,file=fname)
  enddo
  do i=1,eos_np
    write(11,'(4e15.7)') eos_z(i),eos_T(i),eos_z_d(i),eos_T_d(i)
    write(12,'(4e15.7)') eos_z(i),eos_rho(i),eos_z_d(i),eos_rho_d(i)
    do iph=1,nph
      lu=12+iph
      write(lu,'(2e15.7)') eos_z_d(i),eos_Gamma(iph,i)
    enddo
  enddo
  close(11)
  close(12)
  do iph=1,nph
    lu=12+iph
    close(lu)
  enddo
endif
! output updated phase change info
do i=1,nph
   if (Ra.eq.0) then
      if (print_node) write(irefwr,'(i4,6e16.7)') i,gamma(i),phz0(i),phdz(i),pht0(i),0d0,Rb(i)
   else
      if (print_node) write(irefwr,'(i4,7e16.7)') i,gamma(i),phz0(i),phdz(i),pht0(i), gamma(i)*Rb(i)/Ra,Rb(i)
   endif
enddo
!if (print_node) then
!   write(irefwr,*) 'deltaT_dim = ',deltaT_dim,' height_dim= ',height_dim,' rho_dim=',rho_dim
!   write(irefwr,*) 'drho: ',drho_rel(1),drho_ph_d(1),drho_ph_rel(1)
!endif
end subroutine print_final_eos
