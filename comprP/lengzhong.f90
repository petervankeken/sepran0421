subroutine lengzhong()
use sepmodulecomio
use mparallel
use convparam
use coeff
use control
implicit none 
integer, dimension(10) :: lo_index_ph,hi_index_ph
integer :: iph,i,npha,N,iz,lu,find_eos_iz
parameter(N=2000)
real(kind=8), dimension(10) :: rho1,rho2
real(kind=8), dimension(N) :: rho,temp_d
real(kind=8) :: alpha_mod,cp_mod,prespi,bigGamma,dGdpi,dGdpi_d,get_prespi
real(kind=8) :: drhoph,drhophtot,rl,adiabat,Di_dim,Di_frac
real(kind=8) :: rhoph,alpha,cp,z,dz,nph_orig,rho_loc
character(len=80) :: fname

! Correct for choice of Di given on input; 
!!!! ! LZ10 uses alpha=4e-5, cp=1000.
!!!! Di_dim = (4e0/3)*alpha_dim*g_dim*(R2_dim-R1_dim)/(cp_dim/1.25)
!!!! Di_frac = Di/Di_dim
!!!! eos_T_d(1) = 1280+273
!!!! temp_d = eos_T_d(1)
!!!! do iz=1,eos_np-1
!!!!    z = eos_z_d(iz)*1e3
!!!!    dz = (eos_z_d(iz+1)-eos_z_d(iz))*1e3
!!!!    alpha = eos_alpha_d(iz)
!!!!    temp_d(iz+1) = temp_d(iz) + dz*Di_frac*alpha*g_dim*temp_d(iz)/cp_dim
!!!! enddo ! iz
!!!! 
!!!! open(11,file='temp_LZ.dat')
!!!! open(12,file='alpha_LZ.dat')
!!!! open(13,file='rho_LZ.dat')
!!!! do i=1,eos_np
!!!!   eos_T_d(i)=temp_d(i)
!!!!   if (delta1K) then
!!!!      eos_T(i) = eos_T_d(i)  - T0_dim
!!!!   else 
!!!!      eos_T(i) = (eos_T_d(i) - T0_dim) /deltaT_dim
!!!!   endif
!!!!   write(11,*) eos_z_d(i),temp_d(i)
!!!!   write(12,*) eos_z_d(i),eos_alpha_d(i)
!!!!   write(13,*) eos_z_d(i),eos_rho_d(i)
!!!! enddo
!!!! close(11)
!!!! close(12)
!!!! close(13)
!!!! ! update various global parameters
!!!! Ts1_dim = eos_T_d(eos_np) - T0_dim
!!!! Tso_dim = eos_T_d(1) - T0_dim
!!!! write(irefwr,*) 'Leng & Zhong 2010: '
!!!! ! update non-dimensional equivalents for adiabatic T at top and bottom
!!!! Ts1_nondim = Ts1_dim
!!!! Tso_nondim = Tso_dim
!!!! deltaT_dim = 3500
!!!! deltaT_dim_pot = deltaT_dim - (Ts1_dim - Tso_dim)
!!!! deltaT_dim_TBL = 1280
!!!! deltaT_dim_LBL = deltaT_dim - deltaT_dim_TBL - (Ts1_dim-Tso_dim)
!!!! Ra = Ra_orig / deltaT_dim  ! following definition LZ10 eq (4))
!!!! open(422,file='Leng_Zhong.dat')
!!!! write(422,*) ' Tso_dim, Ts1_dim = ',Tso_dim,Ts1_dim
!!!! write(422,*) ' deltaT_dim       = ',deltaT_dim
!!!! write(422,*) ' deltaT_dim_pot   = ',deltaT_dim_pot
!!!! write(422,*) ' deltaT_dim_TBL   = ',deltaT_dim_TBL
!!!! write(422,*) ' deltaT_dim_LBL   = ',deltaT_dim_LBL
!!!! write(422,*) ' Ra               = ',Ra,Ra_orig
!!!! close(422)
!!!! 
end subroutine lengzhong
