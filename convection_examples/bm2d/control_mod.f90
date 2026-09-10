! PvK July 2024
module control

  logical :: print_node,conductive,read_velocity,stokes_is_updated,restart,stokes_only,do_not_do_Stokes,tracers_are_set_up
  integer :: irestart,krestart
  logical :: printmatrix,pedebug,petiming,delta1K,Kequivalent,petest,check_on_curved_elements,temperature_in_C
  integer :: subdivide800,subdivide900
  logical :: no_temperature_solution,do_temperature_solution,noGMT=.true.,compute_curl_v=.false.
  integer :: needed_neighbors,needed_neighbors_neighbors,missed
  logical :: use_tracers_from_nate=.false.,use_tracers_from_cian=.false.,start_from_steady=.false.,cartplume=.false.
  integer :: ntracers_from_nate=0,ntracers_from_cian=0,nres_GMT=120,nstokes_solves=0
  integer :: N_ala_subiter_max,nmax_stokes_sub_iter
  ! RaT: solve with f2=Ra*T, not f2=Ra*T'
  logical :: ala_subiter,debugTbars=.false.,cian_notation=.false.,RaT=.false.,no_mumps_info=.true.,single_Stokes=.false.
  logical :: chiTisinvrhobar=.false.,CY85_Eq22=.false.,add_dynamic_pressure_term_to_work,cian_ALA=.false.,cpisnotcv=.true.
  logical :: altNu=.false.,add_wT_Nu=.false.,HoLiu87=.true.,add_phi_Nu=.false.,add_wT_Nu_quad=.false.,add_phi_Nu_quad=.false.
  logical :: add_P_to_Nu=.true.,pressure_test=.false.,addptermtowork=.false.,stokes_sub_iter
  real(kind=8) :: eps_ala_subiter,arbitrary_p_scale,eps_stokes_sub_iter
  ! Jared Whitehead parameters for Ra=Ra(t)
  logical :: use_varRa,update_phase_rho_via_gradient=.false.,update_phase_rho=.false.
  logical :: update_phase_Tbar=.false.,update_phase_Tbar_via_gradient=.false.
  logical :: rhobar_times_bigGamma,strawblob,compute_work_Tprime,weak_zones=.false.,wz_reduce_Ra=.false.
  real(kind=8) :: eta_wz=1e0_8,eta_plate=1e0_8
  real(kind=8) :: Rat_max
  ! memory allocation by PvK
  integer :: npvkalloc=0

  character(len=120) :: checkmem_command
  real(kind=8) :: resmem_max=10.0
  logical :: red_flag=.false.


  logical :: mumps_matrix_kept=.false.,mumps_matrix_factorized=.false.,Coutput=.false.
  logical :: mumps_matrix_kept8=.false.,mumps_matrix_factorized8=.false.

  integer :: isolmethod8,maxiter8,iprint8,ireler8,ipreco8
  real(kind=8) :: cgeps8,ksp_abs8,ksp_rel8
  integer :: isolmethod9,maxiter9,iprint9,ireler9,ipreco9
  real(kind=8) :: cgeps9,ksp_abs9,ksp_rel9

  integer :: metupw

! bstore.inc
  character(len=80) :: f2name,Tstartfile,UVstartfile
  logical :: netcdf=.true.,output_velocity_solution,outputsphan,output_viscosity_breakdown

  ! cimage/colimage
  real(kind=8) :: ximin,ximax,yimin,yimax
  logical :: make_image
  integer :: ipmax
  integer :: nsteps,nulstep

  ! cperson
  logical pvk,sky,hannah,amy,zhangyi

  ! extrainput
  integer :: iextra_input
  character(len=80) :: file_extrainput

  ! tloop + petime
  real(kind=8) :: petoutstep,tfac,tstart_d
  integer :: noutput,nbetween,ncor,idia

  ! peiter
  integer :: nitermax,nout
  real(kind=8) :: tstepmax,tvalid,difcor
  real(kind=8) :: difcormax
  logical :: print_header

  ! pecpu
  real(kind=4) :: t00=0.0,t1=0.0,t2=0.0,t3=0.0,cput
  real(kind=4) :: dcpu=0.0,cpu_now=0.0,cpu_then=0.0,cpu_first=0.0,cpu_total=0.0,cpustart=0.0
  real(kind=4) :: cpu_heat=0.0,cpu_stokes=0.0,cpu_tracers=0.0,cpu_after_start=0.0
  real(kind=4), external :: second

  ! pedebug
  logical :: pdebug

  ! pesteady
  real(kind=8) :: eps_convergence,relax,relaxtosi,dif1,dif2,dif
  integer :: nsteady_max,isolution_type,niter=0,min_iter_steady
  logical :: steady,did_not_converge

  ! penoniter
  real(kind=8) :: subeps
  integer :: nsubmax,nsub(10)

  ! plotpvk
  character(len=80) :: ourplotname,vtkname
  integer :: inout
 
  ! zhchem
! integer ZHNTRACMAX
! parameter(ZHNTRACMAX=10)
! logical meltable(ZHNTRACMAX)
! integer nonmeltable_index(ZHNTRACMAX/10)

  ! peverbose
  logical verbose

  ! new comprP functionality
  integer, parameter :: lu_nml=101
  integer :: ipetsc8,imatrix8,ipetsc9,imatrix9

  ! logical unit information
  integer, parameter :: LU_NU=21,LU_VRMS=20,LU_CMB_ENTRY=221,LU_CMB_EXIT=222,LU_MELT_START=200,LU_ROTATION=18
  integer, parameter :: LU_INGAS_START=210,LU_TAU_DROP=73,LU_TSTEP=60,LU_HEATFLOW=19,LU_NU_BOT=17,LU_SURFACEVEL=16
  integer, parameter :: LU_WORK=70,LU_PHI=71,LU_MEM=72,LU_SCRATCH=199
  ! used: 16-21, 60, 70-73, 200-222

  integer :: gable_output_choice=0

end module control
