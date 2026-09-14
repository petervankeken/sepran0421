subroutine read_namelist_input() 
use sepmodulecomio
use geometry
use coeff
use mtime
use mpetrac
use convparam
use brandenburg
use sepmodulecomio
use tracers
use mparallel
use eos
use control
implicit none
integer :: i,jtypv,iiqtype,iiiqtype,kfollowchem,mcont
logical :: time_is_nondimensional
real(kind=8) :: dr_min

namelist /comprP_nml/ axi,cyl,steady,nsteady_max,eps,jtypv, &
     & Ra,Rb_local,Di,deltaT_dim,delta1K,Kequivalent,Ts_eos0, T0_dim, &
     & ipetsc8,imatrix8, ipetsc9,imatrix9,isolmethod8,isolmethod9, &
     & nlay,zint_d,q_layer_d,iiqtype,krestart,itop,ibottom,isideboundary,Tstartfile,UVstartfile,tmax_d,tfac,dtout_d, &
     & nitermax,printmatrix, &
     & nbetween,ncor,mcont,itypv,viscl,R1,R2,iclc,icloc,nbp,iboundpoints,ipmax, &
     & pedebug,petiming,intrule800,intrule900,interpol900,interpol800,output_velocity_solution,itype_stokes, &
     & isolmethod8,maxiter8,iprint8,ireler8,cgeps8,ksp_abs8,ksp_rel8,ipreco8, & 
     & isolmethod9,maxiter9,iprint9,ireler9,cgeps9,ksp_abs9,ksp_rel9,ipreco9,time_is_nondimensional,metupw,iextra_input, & 
     & relax,relaxtosi,wavel_perturb,half,quart,eighth,nr_output,nth_output,subtract_rotation, &
     & gable_plates,nplates,whalf,recomputeG,istress_curve,min_iter_steady,stokes_via_user,b_eta,c_eta,subdivide800,subdivide900, &
     & petest, tracerC, fieldC, duplicate_tracers, zmelt, xmelt, d_cmb, track_cmb_ingas, initialize_crust, &
     & itracoption,stop_after_setup_tracers,ndist,kfollowchem,dr,dr_min,npix_radial,x0tr,x1tr,y0tr,y1tr,ainittr, &
     & check_on_curved_elements, &
     & predict_with_RK4,T_bot,gable_output_choice,ibench_type,iadiabat,divufromeos,eos_type,TALA,stokes_only,do_not_do_Stokes, & 
     & drho_background_dense,step_rho_background,T_buoyancy_through_particles,fake_rho_bar, & 
     & no_temperature_solution,ichoice_init_temp, &
     & pvk_buoy_special,compute_curl_v,exp_rho_background,stretch_tracers,tstartp,use_tracers_from_nate,ntracers_from_nate, &
     & nochain,ainit,y0m,bilinearC,ratio_method,truncateC,t_init_limit_max,dt_init_limit,use_tracers_from_cian,ntracers_from_cian, &
     & start_from_steady,nres_GMT,ialphatype,palpha,icondtype,pcond,ivl_smoothstep,old_gable_plates,cartplume,kwave_perturb, &
     & deltah_dim, ampini_perturb,insulbot,ala_subiter,eps_ala_subiter,N_ala_subiter_max, &
     & nph,gamma_d,phz0_d,phdz_d,pht0_d,drho_ph_d,ignore_latent_heat,do_rhoprime,subtract_bigGamma_norm,solve_for_Tperturb, &
     & debugTbars,cian_notation,RaT,Ts_dimK,Tbars_dimK,drho_rel,dCY85,use_effective_alpha_Hr,use_varRa,Rat_max,no_mumps_info, &
     & chiTisinvrhobar,CY85_Eq22,add_dynamic_pressure_term_to_work,arbitrary_p_scale,update_phase_rho_via_gradient, &
     & update_phase_rho,update_phase_Tbar,update_phase_Tbar_via_gradient, &
     & Rb1isRb,rhobar_times_bigGamma,Tbariszero,compute_work_Tprime, &
     & strawblob,xc_blob,yc_blob,width_blob,amplitude_blob,cian_ALA,altNu,single_Stokes,add_wT_Nu,HoLiu87,add_phi_Nu, & 
     & add_wT_Nu_quad,netcdf,&
     & add_phi_Nu_quad,add_P_to_Nu,pressure_test,krad,krad_choice,penalty_parameter,tosi15,tosi15_case,cpisnotcv,addptermtowork, &
     & sigma_y,eta_star,stokes_sub_iter,eps_stokes_sub_iter,nmax_stokes_sub_iter,file_extrainput,weak_zones,wz_reduce_Ra, & 
     & eta_wz,eta_plate,printvis,stop_after_initialize_temp_vel,stop_after_first_real_Stokes


! this is called before sepran starts
stop_after_initialize_temp_vel=.false.
stop_after_first_real_Stokes=.false.
printvis=.false.
eta_wz=1.0_8
eta_plate=1.0_8
itop=3
ibottom=1
isideboundary=2
weak_zones=.false.
wz_reduce_Ra=.false.
stokes_sub_iter=.false.
eps_stokes_sub_iter=1e-2
nmax_stokes_sub_iter=30
eps_ala_subiter=1e-6
N_ala_subiter_max=100
addptermtowork=.false.
cpisnotcv=.false.
penalty_parameter=1e-6
krad_choice=0 ! choices for adding radiative effect to computation of conductivity
krad=.false.
!altNu=.true.  ! use HoLiu et al. 1987 flux computation in top row
HoLiu87=.true.
add_wT_Nu=.true.
add_wT_Nu_quad=.false.
add_P_to_Nu=.true.
cian_ALA=.false.
compute_work_Tprime=.false.
Tbariszero=.false.
xc_blob=0.5_8
yc_blob=0.8_8
width_blob=0.05_8
amplitude_blob=1.0_8
strawblob=.false.
rhobar_times_bigGamma=.false. ! ignore old definition to be compatible with Cian (May/June 2021)
add_dynamic_pressure_term_to_work=.false.
Rb1isRb=.false.
arbitrary_p_scale=1.0_8
chiTisinvrhobar=.false.
CY85_Eq22=.true.
update_phase_rho_via_gradient=.true.
update_phase_rho=.false.
update_phase_Tbar_via_gradient=.false.
update_phase_Tbar=.false.

no_mumps_info=.true. ! avoid writing mumps output to mumps.info
netcdf=.true.
use_varRa=.false.
Rat_max=1e10
use_effective_alpha_Hr=.false.
T0_dim=273.0_8 !  kept for backward compatibility
Ts_eos0=273.0_8 !        ,,
Ts_dimK=273.0_8
Tbars_dimK=273.0_8
RaT=.false.
cian_notation=.false.
debugTbars=.false.
solve_for_Tperturb=.false.
subtract_bigGamma_norm=.false.
do_rhoprime=.false.
eos_type=0
ignore_latent_heat=.false.
gamma_d = 0_8       ! dimensional Clapeyron slope
phz0_d = 0_8        ! dimensional depth of phase change
phdz_d = 0_8        ! dimensional transition thickness of phase change
pht0_d = 0_8        ! dimensional temperature of phase change at phz0_d 
drho_ph_d = 0_8      ! dimensional density contrast across phase change
drho_rel=-1e5      ! old variable indicating relative density contrast
ala_subiter=.false.
iadiabat=2
ialphatype=0
palpha=-2.0_8
icondtype=0
pcond=1.0_8
nres_GMT=120
start_from_steady=.false.
t_init_limit_max=0.0_8
dt_init_limit=0.0_8
truncateC=.false.
kfollowchem=0
ratio_method=.false.
bilinearC=.false.
npix_radial=400
use_tracers_from_nate=.false.
use_tracers_from_cian=.false.
ntracers_from_nate=0
ntracers_from_cian=0
stretch_tracers=.false.
ichoice_init_temp=1
fake_rho_bar=.false.
exp_rho_background=.false.
no_temperature_solution=.false.
do_temperature_solution=.true.
T_buoyancy_through_particles=.false.
drho_background_dense=1.0
step_rho_background=.false.
gable_output_choice=0
predict_with_RK4=.true.
initialize_crust=.false.
d_cmb=0.01
track_cmb_ingas=.false.
duplicate_tracers=.false.
petest=.false.
check_on_curved_elements=.true.
stokes_only=.false.
do_not_do_Stokes=.false.

subdivide800=0
subdivide900=0
b_eta=0.0_8
c_eta=0.0_8
min_iter_steady=0
stokes_via_user=.true. ! better idea for cylindrical models
gable_plates=.false.
min_iter_steady=1
nplates=8
whalf=0.05_8
recomputeG=.true.
istress_curve=43
subtract_rotation=.false.
wavel_perturb=-1 ! makes default wavelength twice the size of the box
relax=0.0_8
relaxtosi=0.0_8
iextra_input=0 ! assume no input beyond comprP.in / comprP.nml
file_extrainput='pvk.nml'
itype_stokes=900
metupw=0
iclc=4
icloc(1:4)=(/1,2,3,4/)
ibottom=1
itop=3
nbp=0
isolmethod8=1 ! CG
maxiter8=10000 
iprint8=-1
ireler8=0
ipreco8=3
cgeps8=1e-8
ksp_abs8=1e-8
ksp_rel8=1e-8
isolmethod9=1 ! CG
maxiter9=10000 
iprint9=-1
ireler9=0
ipreco9=3
cgeps9=1e-8
ksp_abs9=1e-8
ksp_rel9=1e-8
R1=0_8
R2=1.0_8
icoorsystem=0 ! don't adopt code for axisymmetry yet
output_velocity_solution=.true.
delta1K=.true.  ! scale temperature by 1 K
Kequivalent=.false. ! if true temperature is numerically equivalent to dimensional T in K
deltaT_dim=3000.0_8
Di=0.0_8
Rb_local=0.0_8
petiming=.false.
pedebug=.false.
itypv=0
jtypv=0
mcont=0
nlay=1 ! no layered coefficients
do i=1,nlay
   viscl(i)=1.0_8
enddo
zint_d=0.0_8
iiqtype=0
q_layer_d=0.0_8
ncor=1
printmatrix=.false.
axi=.false.
cyl=.false.
full=.true.
conductive=.false.
steady=.true.
nsteady_max=40
eps=1.0e-4_8
Ra=1.0e4_8
ipetsc8=0 ! do not use petsc
ipetsc9=0 ! do not use petsc
imatrix8=2 ! compact matrix
imatrix9=1 ! compact matrix
q_layer=0 ! no internal heating
krestart=0
UVstartfile='UV_start.nf'
Tstartfile='T_start.nf'
read_velocity=.true.
idia=2
tfac=-0.5
tmax_d=1e-3
dtout_d=1e-3
nitermax=10000

! tracer buoyancy
itracoption=0
fieldC=.false.
tracerC=.false.
stop_after_setup_tracers=.false.

! geometry
itop=3
ibottom=1
T_top=0.0_8
T_bot=1.0_8
volume=1.0_8
cyl=.false.
axi=.false.
insulbot=.false.
nr_output=100
nth_output=200

! coeff
icoor900=0
intrule900=0
interpol900=0
interpol800=0
icoor800=0
intrule800=0
compress=.false.

nph=0

read(lu_nml,NML=comprP_nml)
if (cian_ALA.and.CY85_Eq22) then
   write(6,*) 'PERROR(read_namelist): cian_ALA and CY85_Eq22 cannot both be true'
   stop
endif
if (stokes_sub_iter) then
   write(6,*) 'PERROR(read_namelist): stokes_sub_iter needs to be implemented here'
   stop
endif

if (strawblob) then
   ! make sure to update rhobar and Tbar
   update_phase_rho=.true.
   !not yet
   !update_phase_Tbar=.true.
endif

if (cyl) then
   altNu=.false. ! HoLiu et al. 1987 flux method only works in Cartesian for now.
endif

if (nint(Ts_dimK).ne.273) then
   write(6,*) 'PERROR(read_namelist): use Ts_dimK=273 throughout. If you wish to emulate CY85 Ts=1500 or similar'
   write(6,*) '                       specify dCY85=Ts-273 instead'
   stop
endif
! for backward compatibility
if (nint(T0_dim).ne.273) then
   dCY85=T0_dim-Ts_dimK
endif
if (nint(Ts_eos0).ne.273) then
   Tbars_dimK=Ts_eos0
endif
if (abs(drho_rel(1))<1) then
   ! old approach setting relative rather than absolute density difference at phase change
   drho_ph_d=drho_rel*rho_dim
endif
drho_rel=drho_ph_d/rho_dim

deltah_nondim=deltah_dim/height_dim
rlamtemp=2*pi/kwave_perturb

Tbar0=.true. ! assume Tbar=0 except when solving for T', of course
if (solve_for_Tperturb) Tbar0=.false.
bigGammabar0=.true. ! same for phase function Gamma

if (ialphatype>=1.and.palpha>1e-7) then
   if (print_node) then
      write(irefwr,*) 'PERROR(read_name_list_input): palph is positive'
      write(irefwr,*) ' but you probably intend it to be negative...'
      write(irefwr,*) ' palpha=',palpha
   endif
  stop
endif

if (Kequivalent) delta1K=.true.

if (Kequivalent.or.delta1K) then
   deltaT_dimK=1
else
   deltaT_dimK=deltaT_dim
endif

if (delta1K) then
   if (print_node) write(irefwr,*) 'PWARN(read_name_list_input): not yet suited for delta1K for most applications'
   !if (.not.cartplume) stop
endif
if (periodic) then
   if (print_node) write(irefwr,*) 'PERROR(read_name_list_input): not suited for periodic'
   stop
endif

do_temperature_solution=.not.no_temperature_solution

! Set parameters for EOS etc.
pemcont=mcont
Ra_orig=Ra
if (delta1K) then
   Ra=Ra_orig/deltaT_dim
   ! Ra_deltaT_rescale=deltaT_dim
endif
write(irefwr,*) 'Ra: ',Ra,delta1K,deltaT_dim

if (tfac<0 .or. time_is_nondimensional) then
   dtout_d = dtout_d/tscale_dim
   tmax_d = tmax_d/tscale_dim
   tfac = abs(tfac)
endif

! make sure not to stretch tracers with incompressible flow
if (itracoption==1.and.stretch_tracers) then
   if (mcont==0) then
      if (print_node) write(irefwr,*) 'PWARN(read_namelist_input): resetting stretch_tracers for itracoption=1 to F'
      stretch_tracers=.false.
   endif
endif

if (jtypv>0) then
  itypv = jtypv-(jtypv/100)*100
  tackley = (jtypv.ge.100)
else
  itypv = jtypv
  tackley = .false.
endif

write(irefwr,*) 'mcont, imatrix9: ',mcont,imatrix9
if (mcont==1.and.(imatrix9==53.or.imatrix9==54.or.imatrix9==1)) then
   if (print_node) then
      write(irefwr,*) 'PERROR(read_namelist_input): inconsistent input'
      write(irefwr,*) 'mcont = ',mcont
      write(irefwr,*) 'but imatrix9 = ',imatrix9
   endif
   stop
endif

ifollowchem=kfollowchem-(kfollowchem/100)*100
chemwithrho = (kfollowchem/100)*100.eq.1
if (chemwithrho) then
  if (print_node) write(irefwr,*) 'PERROR(compr:chemwithrho): figure this out first'
  stop 
endif
! itracoption   Option for tracer logic. See inittrac()
! 1= distributions of tracers defined by dr and ndist
! 2= one markerchain with marker distance dr
! 9= synthetic tracer test without time integration
! ndist         Number of tracer distributions/markerchains
! ifollowchem   Option to follow chemistry in each tracer
! dr            Average distance between particles
! dr_min          remove markers that are closer to neighbor than dr_min
! chemwithrho   Take compressibility into account?
if (itracoption == 2) then
   nochain=1
   dm(1)=dr
   dm_min(1)=dr_min
endif
if (itracoption == 1 .and. ndist>NDISTMAX) then
   if (print_node) then
     write(irefwr,*) 'PERROR(compress): number of tracer distributions'
     write(irefwr,*) 'too large: ',idist,' > ndistmax = ',NDISTMAX
   endif
  stop
endif

if (print_node) write(irefwr,1112) itracoption,ndist,ifollowchem,dr,dr_min,chemwithrho

qwithrho = (iiqtype/100.eq.1)
iiiqtype = iiqtype - (iiqtype/100)*100
qnondim = (iiiqtype/10).eq.1
iqtype = iiqtype - (iiqtype/100)*100 - (iiiqtype/10)*10
if (iqtype == 2 .and. cyl) then
   if (print_node) then
     write(irefwr,*) 'Need to modify code for layer dependent heating'
     write(irefwr,*) ' in cyl '
   endif
   stop
endif
if (iqtype < 0 .or. iqtype > 5) then
   if (print_node) then
     write(irefwr,*) 'PERROR(compr_start): iqtype has incorrect value'
     write(irefwr,*) ' iqtype = ',iqtype
   endif
   stop
endif

if (iqtype==1) then
   if (qnondim) then
      q_layer=q_layer_d
   else
      q_layer=q_layer_d*QBSE_DIM*height_dim*height_dim/(cp_dim*rkappa_dim*deltaT_dim)
   endif
endif

if (print_node) write(irefwr,920) iqtype,qwithrho,qnondim,q_layer_d(1),q_layer(1)


if (nlay > 1 .and. zint_d(1) > 1.0_8) then
!  if zint_d(1)>1 it's surely dimensional
   ! Convert to non-dimensional depth
   if (axi)  then
      do i=1,nlay-1
         r_zint(i) = (R2_dim - zint_d(i))/height_dim
      enddo
   else if (cyl) then
      ! apply shrunken core scaling
      do i=1,nlay-1
         r_zint(i) = (R2_dim - zint_d(i))*(R2_dim - zint_d(i)) * r2/(R2_dim*R2_dim)
      enddo
   else ! Cartesan 
      do i=1,nlay-1
         r_zint(i) = (R2_dim - zint_d(i) - R1_dim)/height_dim
      enddo

   endif
else if (nlay > 1) then
!        !!! Old non-dimensional input - update r_zint
   do i=1,nlay-1
      r_zint(i) = r2-zint_d(i)
   enddo
endif

if (eos_type==0.and.nph>0.and.compress) then
   if (print_node) then
      write(irefwr,*) 'PERROR(read_namelist): nph>0 but eos_type=0 and compress'
      write(irefwr,*) 'comprP code not suited for phase changes without table look up'
   endif
   call instop
endif
   

if (print_node) write(irefwr,NML=comprP_nml)

if (print_node) write(irefwr,910) nlay,(zint_d(i),i=1,nlay-1)
if (print_node) write(irefwr,911) (r_zint(i),i=1,nlay-1)
if (print_node) write(irefwr,1115) Di,Grueneisen,deltaT_dim_pot,nph

if (myid==1) write(irefwr,*) 'myid1: Ra,Rb_local=',Ra,Rb_local

if (Rb_local<1e-5) CY85_Eq22=.false.


920   format('iqtype ................................... ',i10,/, 'Multiply heatproduction with rho(z)? ..... ',L10,/, &
     &       'Is q_layer nondimensional? ............... ',L10,/,'q_layer_d/q_layer .......................... ',2f8.3:)
910   format('nlay ..................................... ',i8,/,8x, 'zint_d ................................... ',10e10.3:)
911   format(8x,'r_zint ................................... ',10f10.3:)
1112  format(/,'Tracer option .................... ',i10,/, 'Number of tracer distributions.... ',i12,/, &
     &         'Choice for following chemicals.... ',i12,/, 'Distance between particles ....... ',f12.7,/, &
     &         'Minimum distance between particles ',f12.7,/, 'Take compressiblity into account?  ',L12)
1115  format(/'Extended Boussinesq parameters: ',/, '    Di    ..................... ',f8.3,/, &
     &        '    Grueneisen parameter ...... ',f8.3,/, &
     &        '    deltaT_dim ................ ',f8.3,/, &
     &        '    Number of phase changes ... ',i8)



end subroutine read_namelist_input
