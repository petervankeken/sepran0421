! PvK July 2024
module coeff
real(kind=8) :: Ra,Di,Grueneisen=1.0,DiG,DiRa,Ra_orig,Rb_local,penalty_parameter
integer :: iqtype,ichoice_init_temp  
integer, parameter :: NLAYMAX=10
real(kind=8) :: q_layer(NLAYMAX),q_layer_d(NLAYMAX)
logical :: compress=.false.,TALA=.false.,EBA,BA,divufromeos,quad903,quad903_use_funcc3,quad800,ivl_smoothstep=.false.
logical :: subtract_rotation,qwithrho,qnondim,stokes_via_user,T_buoyancy_through_particles,fake_rho_bar
logical :: pvk_buoy_special=.false.,bilinearC,use_effective_alpha_Hr

real(kind=8) :: etamin=0.0_8,etamax=0.0_8

logical :: krad=.false.
integer :: krad_choice=0

integer :: intrule900,icoor900,icoorsystem,interpol900
integer :: intrule800,icoor800,itype_stokes,interpol800
integer :: itypv,mcontv,pemcont,eos_type
logical :: printvis=.false.
real(kind=8) :: viscl(NLAYMAX)=1.0_8
logical :: ivl,ivt,ivb,ivn,dostokes
integer,parameter :: EOS_NPMAX=1351
integer :: eos_np
logical :: eos_exists
character(len=80) :: eos_data_base
real(kind=8),dimension(EOS_NPMAX) :: eos_p_d,eos_T_d,eos_alpha_d,eos_cp_d,eos_K_d,eos_rho_d, &
    & eos_g_d,eos_z_d,eos_alpha,eos_cp,eos_K,eos_rho,eos_T,eos_z
real(kind=8),dimension(10,EOS_NPMAX) :: eos_Gamma
real(kind=8) :: eos_pav,eos_Tav,eos_alphaav,eos_cpav,eos_Ta_av, &
     & eos_Kav,eos_rhoav,eos_gav,eos_zav,eos_condav
real(kind=8) :: rho_av,r_rho0,z_rho0,rho_z0,rhomax,drho_background_dense
logical :: step_rho_background,exp_rho_background,stretch_tracers,Rb1isRb,Tbariszero


logical :: tosi15=.false.,tackley=.false.,tackley2000=.false.,tackley_newvis=.false.
integer :: ibench_type,tosi15_case ! ,nmax_stokes_sub_iter
! ibench_type=2  : van Keken 2001 (PEPI01)
! ibench_type=3  : King et al. 2010;  (King2010)
! ibench_type=4  : RT benchmark from van Keken et al. 1997 (JGR97)
! ibench_type=5  : thermochemical benchmark from van Keken et al. 1997 (JGR97)
! ibench_type=7  : cylindrical benchmark from van Keken 2001 (??)
! ibench_type=8  : special test of Leng & Zhong 2008 set up
! ibench_type=101: Stokes sphere test 
! ibench_type=300: Stokes cpu test
real(kind=8) :: b_eta,c_eta,Tshift,CH94_b,CH94_c,z_yieldlimit,etalin_z_dep
real(kind=8) :: eta_star,sigma_y,dnu_T,visc_factor,sigma_b
real(kind=8) ::  plume_viscosity_min,visminmin,visminmax
real(kind=8) :: plume_viscosity_max
real(kind=8) :: etalinh,etaplasth
logical :: output_vislin,output_visplast

logical :: its_CH94,CH94_fs,stokes_subdivide

integer :: ibuoy_trac
logical tracerC,fieldC,ratio_method,old_rhsd,ratio_method_overlay
logical truncateC,no_buoyancy,duplicate,test_C_through_T

logical :: vara
real(kind=8) :: Ra_final,Ra_start

! cpephase.inc
!     *** Ra       - Thermal Rayleigh number
!     *** Rb       - Phase buoyancy number used in heat equation
!     *** Rb1      - Modified Phase Buoyancy numbers for heat equation in
!     compressible convection
!     *** Di       - Dissipation number
!     *** gamma    - CC slope phase change
!     *** gRbRa    - gamma(1)*Rb/Ra
!     *** g2RbRaDi = gamma(1)*gRbRa*Di
!     *** phz0     - depth phase change (at T=0)
!     *** phdz     - depth interval of phase change
!     *** nph      - number of phase changes in model
! phase_rho: rho1 above and rho2 below phase change iph
integer, parameter :: NPHASE_MAX=10
real(kind=8),dimension(NPHASE_MAX) :: Rb,Rb1,Rb_rel,Rr,gamma,gamma1,phz0,phdz,pht0,drho_ph_d,drho_ph_rel,drho_rel
real(kind=8),dimension(NPHASE_MAX) :: gamma_d,phz0_d,phdz_d,pht0_d
real(kind=8),dimension(NPHASE_MAX) :: glRbRa,gl2RbRaDi,sign_gamma
real(kind=8),dimension(NPHASE_MAX,2) :: phase_rho
integer :: nph=0
logical :: absolute1,absolute2,solve_for_Tperturb
logical :: keep_deltaT_dim,do_rhoprime

! strawblob
real(kind=8) :: xc_blob,yc_blob,width_blob,amplitude_blob


end module coeff
