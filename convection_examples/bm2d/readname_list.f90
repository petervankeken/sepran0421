! PvK July 2024
subroutine readname_list
use sepmodulecomio
use control
use coeff
use geometry
use convparam
implicit none
integer :: open_status,jtypv
namelist /bm2d_nml/ isolmethod9,jtypv,intrule900,intrule800,interpol900,interpol800,icoor800,icoor900, &
     & mcontv,penalty_parameter,cyl,Ra,wavel_perturb,ampini_perturb,itype_stokes,irestart,Tstartfile,UVstartfile, &
     & eps_convergence,b_eta,c_eta,relax,Di,nsteady_max,r1,r2,tosi15,tosi15_case,steady,sigma_y,eta_star,read_velocity, &
     & itop,ibottom,ipetsc8,ipetsc9,isolmethod8,isolmethod9,imatrix8,imatrix9

irestart=-1 ! start from conductive + perturbation
itype_stokes=900 ! set penalty function method as default 
isolmethod9=0 ! direct solution method
itypv=0 ! isoviscou
intrule900=0 ! let sepran decide on integration rule
intrule800=0
interpol900=0 ! no interpolation P2->P1
interpol800=0
icoor800=0 ! Cartesian coordinates
icoor900=0
mcontv=0 ! Boussinesq
penalty_parameter=1e-6_8 ! penalty function parameter
cyl=.false. ! Cartesian 
print_node=.true.  ! make sure to modify this in parallel
Ra=1e4  
wavel_perturb=1.0_8
ampini_perturb=0.05_8
Tstartfile='T_start.nf'
UVstartfile='UV_start.nf'
eps_convergence=1e-4_8
relax=0.0_8
Di=0.0_8
nsteady_max=20
r1=0.0_8
r2=1.0_8
steady=.true.
read_velocity=.false.
itop=3
ibottom=1

itypv=0
jtypv=0
ivl=.false.
ivt=.false.
ivn=.false.
tosi15=.false.
tosi15_case=0
tackley=.false.
b_eta=0.0_8
c_eta=0.0_8
sigma_y=0
eta_star=0
ipetsc8=0
ipetsc9=0
imatrix8=2
imatrix9=1
isolmethod9=0
isolmethod8=0


open(lu_nml,file='bm2d.nml',iostat=open_status)
if (open_status /= 0) then
   if (print_node) write(irefwr,*) 'PERROR(readname_list): error on opening bm2d.nml: ',open_status
   call instop
endif
read(lu_nml,NML=bm2d_nml)
close(lu_nml)

if (jtypv>0) then
   itypv=jtypv-(jtypv/100)*100
   tackley=(jtypv>=100)
else
   itypv=jtypv
   tackley=.false.
endif



!write(6,*) 'penalty_parameter: ',penalty_parameter

end subroutine readname_list
