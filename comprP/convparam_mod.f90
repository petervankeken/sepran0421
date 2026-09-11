module convparam
! dimensional.inc
real(kind=8) :: alpha_dim,rho_dim,cp_dim,Ts_dimK=273.0_8,rk_dim,Tbars_dimK=273.0_8,dCY85=0.0_8,dCY85_nondim
real(kind=8) :: rkappa_dim,g_dim,eta0_dim,rhoc_dim
real(kind=8) :: cpc_dim,R2_dim,R1_dim,deltaT_dim,dtout_d,tmax_d
real(kind=8) :: year_dim,tscale_dim,height_dim,QBSE_dim
real(kind=8) :: Ks_dim,deltaT_dimK
real(kind=8) :: velocity_scale,pressure_scale,deltaT_dim_pot
real(kind=8) :: T0_nondim,Tm0_dim,Tm0_nondim,deltaT_dim_LBL
real(kind=8) :: Ts_eos0,deltaT_dim_TBL
real(kind=8) :: etamax_TBL,T_top_orig,T_bot_orig,T_bot_pot
real(kind=8) :: Tjump_tot_nondim,Tjump_pot_nondim,T_bot_adia
real(kind=8) :: Tjump_tot_dim,Tjump_pot_dim,T_top_adia
logical :: use_etamax_TBL,subtract_bigGamma_norm,Tbar0,bigGammabar0
! Make sure we record the various temperatures in K equivalents
real(kind=8) :: Ts_nondimK,Tbars_nondimK
real(kind=8) :: T0_dim=273.0_8,T_bc_top
! deltaT_dim_LBL/TBL : lower and top boundary layer temperature differences (used in plume models)

! surface values of alpha, diffusivity, bulk modulus
parameter(alpha_dim=3e-5,Ks_dim=1e11)
parameter(rkappa_dim=1e-6)
! Average value for conductivity (since we used k=(rho/rhoav)**pdif)
parameter(rk_dim=7.37)
! reference mantle density 
parameter(rho_dim=3300d0)
! dimensional radii of CMB and Earth's surface. 
parameter(R1_dim=3493e3,R2_dim=6371e3,height_dim=R2_dim-R1_dim)
! density and specific heat in the core
parameter(rhoc_dim=11000d0, cpc_dim=500d0)
! Bulk Silicate Earth heating 
parameter(QBSE_dim=4.8e-12)
! gravity
parameter(g_dim=9.8d0)
! specific heat
parameter(cp_dim=1250d0)
! scale for mantle viscosity
parameter(eta0_dim=1e22)
! year in seconds
parameter(year_dim=365.25*3600*24)
! tscale_dim translates Byr to non-dimensional units
parameter(tscale_dim=rkappa_dim*year_dim*1e9/(height_dim*height_dim))
! velocity scale (m/s)
parameter(velocity_scale=rkappa_dim/height_dim)
! pressure scale (Pa)
parameter(pressure_scale=rkappa_dim*eta0_dim/(height_dim*height_dim))

! depth_thermodyn.inc
integer :: icondtype,ialphatype,ibetatype,iadiabat
real(kind=8) :: pcond,palpha,geometry_factor,CH94_d,CH94_s,C_peridotite,ainittrac_CH94
real(kind=8) :: alphamin1,alphamax1,u_nought
namelist /CH94_nl/ CH94_d, CH94_s, u_nought, C_peridotite, ainittrac_CH94

! plume.inc
logical :: hannah_test_mesh
integer :: symcurve,nstep_hannah
! func(): initial condition for plume models
! deltah_dim: boundary layer thickness in km (130)
! kwave_perturb: wave number of perturbation (16)
! ampini_perturb: amplitude of perturbation (0.5)
real(kind=8) :: deltah_dim=130.0_8, kwave_perturb=16,ampini_perturb=0.5,plume_tracer_start
real(kind=8) :: deltah_nondim,rlamtemp
real(kind=8) :: r_plume_max,deltaT_top,deltaT_bot,Tso_plume_dim
real(kind=8) :: r_low_cond
integer :: plume_eos_type
logical :: hannah_stop,plume_with_top_BL,plume_model
logical :: ignore_latent_heat,eos_ph_latent_heat
logical :: eos_prep_only,low_cond_top,zero_coarse_vel
logical :: restoreUMplume,zero_uppermantle_vel
logical :: zero_vel_700km,zero_vel_400km
real(kind=8) :: vis_ord

integer :: plume_bench_type
namelist /plume_nml/ vis_ord,deltah_dim,kwave_perturb, &
     &                      ampini_perturb,r_plume_max, &
     &                      plume_tracer_start,r_low_cond, &
     &                      plume_bench_type,symcurve,plume_eos_type, &
     &                      nstep_hannah,hannah_test_mesh, &
     &                      ignore_latent_heat,eos_ph_latent_heat, &
     &                      eos_prep_only,low_cond_top,zero_coarse_vel, &
     &                      restoreUMplume,zero_uppermantle_vel, &
     &                      zero_vel_700km,zero_vel_400km, &
     &                      use_etamax_TBL

real(kind=8) :: rho_u_0,rho_u_tr,rho_l,rho_ecl0,rho_ecl_cmb
real(kind=8) :: dep1,drho1,dep2,drho2,drho_s1,drho_s2
logical :: sky_filter

! degas.inc
!     *** decay times in per By
real(kind=8),parameter :: rl238=1.55e-1,rl235=9.85e-1,rl232=4.95e-2,rl40=5.54e-1
! time scale in By
real(kind=8),parameter :: tscale=2.637e2
! initial concentrations (in 1e14 atoms/g)
real(kind=8),parameter :: U238_0=1.0674,U235_0= 0.32813,Th232_0=2.4322,K40_0=56.9129,K40_eff=0.1024
integer :: ndegaszone(10)
integer,parameter :: NSPECIES=5,NCHEMCASE=4,NDEGAS_DEPTHS=4,MAXDEGASZONE=4
real(kind=8) :: thmindegas(NSPECIES,MAXDEGASZONE),thmaxdegas(NSPECIES,MAXDEGASZONE)
real(kind=8) :: rdegas(NDEGAS_DEPTHS),zdegas(NDEGAS_DEPTHS)
integer,dimension(NCHEMCASE) :: ndegas,nnow,ncont,nnowc
integer :: numdegas_depth,ifollowdegas
logical :: vardegaszones

real(kind=8) :: U238,U235,Th232,K40,He4,Ar40
real(kind=8) :: dU238,dU235,dTh232,dK40,dHe4,dAr40
real(kind=8) :: He4loss(4),He3loss(4),Ar40loss(4)
real(kind=8) :: U238loss(4),U235loss(4),Th232loss(4),K40loss(4)

end module convparam
