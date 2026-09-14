! Solve Stokes equation with petsc
! This approach follows closely what we use in the production codes (e.g., comprS)

! PG=Programmer's Guide (ta.twi.tudelft.nl/NW/sepran/pg.pdf)
! SP=Standard Problems Guide (ta.twi.tudelft.nl/NW/sepran/sp.pdf)

! PvK August 2017
! PvK October 2017: expanded for parallel computing
subroutine comprP_start
#ifdef MPI
use mpi
#endif
use sepmoduleoldrouts
use sepmodulecomio  ! provides lu's for reading & writing
use sepmodulecpack  ! parallel computing interface
use sepmodulekmesh
use sepmoduleintmat
use sepran_arrays ! access to kmesh, kprob, isol, etc. and interfaces to subroutines
use sepran_interface 
use control
use geometry
use mtime
use coeff
use eos
use brandenburg
use tracers
use convparam
use mpetrac
use mparallel
use dtm_elem  ! controls building of pressure mass matrix
implicit none
integer :: modelv,mconv,nparm,nsup1,idum,iread,ihelp,icurvs(2),ierr
real(kind=4) :: timeA,timeB
!real(kind=8) :: u1lc(3),solve_rin(100)
!integer :: commat_in(14),build_in(10),iu1lc(3),solve_in(14)
integer :: idegfd,i,ichoice,isubr,jtypv
real(kind=8) :: umax,outputtime,tmin2,tmax2,p,q,rinputcr(20),vrms,vrms_d,rho_n,Rb_here
character(len=120) :: namedof(4),makedir
integer :: ichvc,iprob_here,ip,ipbuoy2,ipuser,inputcr(20),i1
logical :: namelist_input,postprocessing=.false.,noGMT_store,GMT
integer :: commat_in(20)=0,ipoint,alloc_status,numpar
real(kind=8) :: get_resmem


cpu_first=second()
irefwr=6 ! before start() is called
print_node=.true.

inquire(file='comprP.nml',exist=namelist_input)
if (.not.namelist_input) then
   if (print_node) write(irefwr,*) 'PERROR(comprP_start): comprP.nml is missing'
   stop
endif
   
if (namelist_input) then
   open(lu_nml,file='comprP.nml')
   call read_namelist_input()
   close(lu_nml)
endif

! Check whether MPI parallelism is used outside of Sepran
call check_parallel
write(lu_out,*) 'hi there! ',myid

if (solve_for_Tperturb.and.ibench_type/=3) then
   if (print_node) then
      write(irefwr,*) 'PERROR(comprP_start): solve_for_Tperturb is only allowed for King et al. 2010 benchmark'
      write(irefwr,*) '    but ibench_type = ',ibench_type
   endif
   call instop
endif

if (Di.le.1e-7.and.(icondtype==1.or.ialphatype==1)) then
   if (print_node) then
      write(irefwr,*) 'PERROR(comprP_start): Di=0: ',Di,'but icondtype or '
      write(irefwr,*) ' ialphatype are set to 1: ',icondtype,ialphatype
   endif
   call instop
endif

! define BA/EBA/ALA and various EOS parameters
call set_up_eos

if (itracoption==0) then
   tracerC=.false.
   fieldC=.false.
else
   tracerC=.not.fieldC
endif


if (axi.and.cyl) then
   if (print_node) then
      write(irefwr,*) '***************************************'
      write(irefwr,*) '*   AXISYMMETRIC SPHERICAL GEOMETRY   *'
      write(irefwr,*) '***************************************'
   endif
else if (.not.axi.and.cyl) then
   if (print_node) then
      write(irefwr,*) '****************************'
      write(irefwr,*) '*   CYLINDRICAL GEOMETRY   *'
      write(irefwr,*) '****************************'
   endif
else if (axi.and..not.cyl) then
   if (print_node) then
      write(irefwr,*) '****************************'
      write(irefwr,*) '*   AXISYMMETRIC GEOMETRY  *'
      write(irefwr,*) '****************************'
   endif
else
   if (print_node) then
      write(irefwr,*) '****************************'
      write(irefwr,*) '*   CARTESIAN GEOMETRY  *'
      write(irefwr,*) '****************************'
   endif
endif

read_velocity=.true.
if (krestart<0) then
   irestart=krestart
else
   if (krestart >= 1000) read_velocity=.false.
   irestart = krestart - 1000*(krestart/1000)
endif   
restart = (irestart==1).or.(irestart==2)

if (print_node) write(irefwr,*) 'restart  = ',restart

tmaxp = tmax_d * tscale_dim
dtoutp = dtout_d * tscale_dim
if (print_node) then
   write(irefwr,'(''Dimensional model time (Byr) : '',f15.7)') tmax_d
   write(irefwr,'(''Non-dimensional model time   : '',f15.7)') tmaxp
   write(irefwr,'(''Dimensional output time (Byr): '',f15.7)') dtout_d
   write(irefwr,'(''Non-dimensional output time  : '',f15.7)') dtoutp
endif

if (restart) then
   tstart_d = tstartp / tscale_dim
   if (print_node) then
      write(irefwr,'(''Dimensional start time (Byr) : '',g15.7)') tstart_d
      write(irefwr,'(''Non-dimensional start time   : '',g15.7)') tstartp
   endif
   time_now=tstartp
else
   time_now=0.0_8
endif

if (print_node) then
   call system("mkdir -p images GMT PLOTS VTK")
   if (netcdf) then
      call system("mkdir -p solutions")
   else
      if (mpi_partrac) then
         if (print_node) write(irefwr,*) 'adjust output to old tracer files for mpi_partrac'
         call instop
      else
         call system("mkdir -p Tracers_old solutions_old")
      endif
   endif
endif
if (mpi_partrac) then 
   write(makedir,'(''mkdir -p Tracers/'',i3.3)') myid
   call system(makedir)
else
   call system("mkdir -p Tracers")
endif

if (intrule900==1) then
   if (print_node) then
      write(irefwr,*) 'PWARN(comprP_start) :: intrule900=1'
      write(irefwr,*) 'which will lead to inaccurate results with sepran version 0619 or 1219'
   endif
   !stop  ! not instop here as Sepran has not been initialized
endif
     

! Make sure arrays are initialized
iuser_here=0
iinput=0
! List declared length in first entries of iuser and user
iuser_here(1)=NUSERMAX
iuser_here(2)=1 ! indicate that there is only one element group

! default settings for commat
!ipetsc=0  ! use petsc solve
!imatrix=6  ! use sepran build of matrix before generating petsc matrix
            ! currently only working option in parallel

#ifdef MPI
if (mpi_parallel) call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
if (print_node) then
   write(irefwr,*) 
   write(irefwr,*) 'PINFO(comprP_start): calling sepstr(); check sep_output.000 for any errors when running under mpi_parallel'
endif
! sepstr initiates sepran, reads mesh, defines problem and matrix structure
call sepstr(kmesh1,kprob1,intmat1)
! rerout stdout to lu 6 for process 0 (rest write to sep_output.XXX via lu_out
if (myid==0) irefwr=6
!if (myid==0) write(6,*) irefwr
!call instop

call get_checkmem_command
!write(6,*) get_resmem()
!call instop

! force renumbering
call topol(2)
matr1=0
! set print_node to avoid multiple outputs to stdout for sepran parallel
if (parallel .and. inodenr /= 1) print_node=.false.

call setsepranmodulevalues()  ! adjust a few default values 
! find dimension of space and number of nodal points
call pvk_allocate(user_here,100+7*npoint,'user_here')
call pvk_allocate(funcx,6+ndim*npoint,'funcx')
call pvk_allocate(funcy,6+ndim*npoint,'funcy')
funcx=0e0_8
funcy=0e0_8
user_here=0.0_8
user_here(1)=1.0_8*(100+7*npoint)
funcx(1)=1.0_8*(6+ndim*npoint)
funcy(1)=1.0_8*(6+ndim*npoint)
if (gable_plates) then
   allocate(GP(nplates,nplates),TP(nplates))
endif

if (itracoption>0) call sepfillkmesho

!if ((ipetsc9/=0.or.ipetsc8/=0).and.mpi_partrac) then
!   if (print_node) then 
!      write(irefwr,*) 'PERROR(comprP_start): use of Petsc with mpi_partrac creates segfaults'
!      write(irefwr,*) 
!      write(irefwr,*) 'for now at least...'
!   endif
!   call instop
!endif

!if (ipetsc9/=2) then
!   if (print_node) write(irefwr,'(''PERROR(comprP_start): do only use petsc/mumps for Stokes until'', &
!      & '' errors are fixed'')') 
!   call instop
!endif
   
  

! Redefine matrix structure  and create intmat2
commat_in=0
commat_in(1)=14  ! number of entries in commat_in
commat_in(2)=imatrix9 ! row compact storage for petsc (55=row compact; 56=row compact, no sepran mat)
commat_in(3)=0
commat_in(4)=0
commat_in(5)=1  ! iprob = problem number 
commat_in(13)=ipetsc9  ! use_petsc=petsc (2=mumps)
if (print_node) write(irefwr,'(''commat: '',14i5)') commat_in(1:14)
call matstruc(commat_in,kmesh1,kprob1,intmat1)
nested=.false.
commat_in(1)=14  ! number of entries in commat_in
commat_in(2)=imatrix8 ! row compact storage for petsc (55=row compact; 56=row compact, no sepran mat)
commat_in(3)=0
commat_in(4)=0
commat_in(5)=2  ! iprob = problem number 
commat_in(13)=ipetsc8  ! use_petsc=petsc (2=mumps)
if (print_node) write(irefwr,'(''commat: '',14i5)') commat_in(1:14)
call matstruc(commat_in,kmesh1,kprob1,intmat2)
! Emergency intervention for change in 0421. For some reason normal 
! profile matrix is seen as nested (see sepfillmatrstr1)
nested=.false.

! namelist dtm_elem_nml default
dtm_pressure_mass_matrix=.false.
dtm_only_pressure_mass_matrix=.false.
print_matrix = .false.
do_solve=.true.
! read namelist dtm_elem_nml
inquire(file='dtm_elem.nml',exist=do_dtm_elem)
if (do_dtm_elem) then
  open(11,file='dtm_elem.nml')
  read(11,NML=dtm_elem_nml)
  close(11)
  if (parallel .and. print_node) write(irefwr,NML=dtm_elem_nml)
endif


! at this stage we know whether we will use buoyancy of the tracers (tracerC)
! and have a mesh defined 
alloc_status=0
if (tracerC) then
  ! allocate necessary arrays for evaluation of f1 and f2 for
  ! Stokes problem
  numpar=1
  if (mpi_partrac.and. numprocs>1) numpar=2
  if (cyl) then
     allocate(f1_tracrhsd(7,nelem,numpar),stat=alloc_status)
     if (alloc_status /= 0) then
       write(irefwr,*) 'comprP_start: allocate on f1 failed: ',alloc_status,myid
       write(irefwr,*) 'nelem = ',nelem,myid
       write(irefwr,*) 'numpar = ',numpar,myid
       write(irefwr,*) 'cyl = ',cyl,myid
       call instop
     endif
     f1_tracrhsd=0.0_8
  endif
  allocate(f2_tracrhsd(7,nelem,numpar),stat=alloc_status)
  if (alloc_status /= 0) then
    write(irefwr,*) 'comprP_start: allocate on f2 failed: ',alloc_status,myid
    write(irefwr,*) 'nelem = ',nelem,myid
    write(irefwr,*) 'numpar = ',numpar,myid
    write(irefwr,*) 'cyl = ',cyl,myid
    call instop
  endif
  f2_tracrhsd=0.0_8
  if (print_node) write(irefwr,*) 'allocated f1/f2: ',nelem,numpar
endif

! The input above provides a relative density jump and dimensional values
! for Clapeyron slope, depth and thickness of the phase changes etc.
! Non-dimensionalize and set up Rb as a function of Ra and density change
if (nph.gt.0) then
   if (myid==0) then
      write(6,'(a4,7a16)') 'i','gamma','phz0','phdz','pht0','P','Rb','Rb1/Rb'
   endif
endif
do i=1,nph
!  Phase buoyancy number used in heat equation.
   drho_ph_rel(i)=drho_ph_d(i)/rho_dim
   Rb(i) = Ra*drho_ph_rel(i)/(alpha_dim*deltaT_dimK)
!  non-dimensional Clapeyron slope
!  NB: height_dim == R2_dim - R1_dim
   gamma(i)=gamma_d(i)*deltaT_dimK/(rho_dim*g_dim*height_dim)
!  non-dimensional thickness of phase transition
   phdz(i) = phdz_d(i)/height_dim
!  non-dimensional temperature is set above
!  non-dimensional depth of phase transition (make sure to
!  perform correct scaling for cylinder geometry; see
!  eq. 6 in Van Keken, PEPI, 2001.
   if (cyl.and..not.axi) then
      phz0(i)=r2*(1-(R2_dim-phz0_d(i))**2/(R2_dim**2))
   else  if (cyl.and.axi) then
      phz0(i)=r2-(R2_dim-phz0_d(i))/height_dim
   else
!     Cartesian
      phz0(i)=phz0_d(i)/height_dim
   endif
enddo
if (Kequivalent) then 
   pht0=pht0_d+273.0_8 ! not Ts_dimK as that will confuse the solution for CY85 with Ts>273
else if (delta1K) then
   pht0=pht0_d
else
   pht0=pht0_d/deltaT_dim
endif
! Rb1 is used in heat equation for compressible convection with phase changes
Rb1=Rb

!write(6,*) 'Rb    : ',Rb(1:2)
!write(6,*) 'Rb1   : ',Rb1(1:2)

if (compress) then
   if (print_node) write(irefwr,*) 'call update_rho_adia'
   if (ibench_type==8) then
!     Leng and Zhong 2010
      call update_rho_adia(-1)
   else
!     *** update reference density and adiabat if necessary
      if (eos_ph_latent_heat) then
         call update_rho_adia(2)
      else
         call update_rho_adia(1)
      endif
   endif
endif

!??!??pht0=0.0 ! follow definition of CY85 (4)

call print_final_eos()
! write(irefwr,*) 'stopping after print_final_eos'
! call instop


jtypv=abs(itypv)
ivl=jtypv.eq.1.or.jtypv.eq.3.or.jtypv.eq.5.or.jtypv.eq.7
if (ivl) jtypv = jtypv-1
ivt = jtypv.eq.2.or.jtypv.eq.6
if (ivt) jtypv = jtypv-2
ivn = jtypv.eq.4
ivb = jtypv.eq.10
if (print_node) then
   write(irefwr,*) 'Type of viscosity law: ',itypv
   write(irefwr,*) 'IVL, IVT, IVN: ',ivl,ivt,ivn,ivb
!  write(11,*) 'Type of viscosity law: ',itypv
!  write(11,*) 'IVL, IVT, IVN: ',ivl,ivt,ivn,ivb
endif

xi_eta_stored=.false.


! Create velocity and temperature vectors
idum=1
u1lc=0.0_8  ! all components are 0
iu1lc=0
ichvc=1  ! type of solution vector
iprob_here=1  ! first problem is velocity
call creavc(0,ichvc+(iprob_here-1)*1000,idum,isol1,kmesh1,kprob1,iu1lc,u1lc)
ichvc=1  ! type of solution vector
if (print_node) write(irefwr,*) 'ichoice_init_temp: ',ichoice_init_temp
iu1lc(1)=ichoice_init_temp  ! set dof 1 to func(1,x,y,z)
iprob_here=2  ! second problem is temperature
call creavc(0,ichvc+(iprob_here-1)*1000,idum,isol2,kmesh1,kprob1,iu1lc,u1lc)
!write(irefwr,*) 'temp done'
iu1lc=0

! Now add vectors for adiabat, density, etc.
! Do this before filling viscosity vector or before solving any equations as these arrays 
! may be required
call create_additional_vectors()
!write(6,*) 'stopping after create_additional_vectors'
!call instop

if (itracoption==1) then
   iu1lc=0
   u1lc=0.0_8
   ichvc=1
   iprob_here=2
   call creavc(0,ichvc+(iprob_here-1)*1000,idum,idens,kmesh1,kprob1,iu1lc,u1lc)
   if (fieldC) then
      call creavc(0,ichvc+(iprob_here-1)*1000,idum,iarea,kmesh1,kprob1,iu1lc,u1lc)
      call create_area_array()
   endif
   if (ratio_method) then
      if (print_node) write(irefwr,*) 'update creating idens for ratio_method'
      !call instop
      call creavc(0,ichvc+(iprob_here-1)*1000,idum,idens2,kmesh1,kprob1,iu1lc,u1lc)
   endif
endif


inquire(file='noGMT',exist=noGMT)
if (print_node) write(irefwr,*) 'noGMT: ',noGMT
inquire(file='GMT',exist=GMT)
if (print_node) write(irefwr,*) 'GMT: ',GMT
noGMT=.not.GMT

call determine_geometry(iuser_here,user_here)

! Set up the temperature boundary condition information
! Proper initialization of temperature happens below in initialize_temp_vel
! The total temperature contrast depends on potential temperature jump and
! adiabatic increase. The adiabat has been set in set_up_eos and may have been
! updated by update_rho_adia.
call setup_temperature_bc(iuser_here,user_here)

!write(6,*) 'stopping after setup_temperature_bc'
!call instop




if (print_node) then
     write(irefwr,'('' Total T contrast = '',f12.3)') Tjump_tot_nondim
     write(irefwr,'('' Tpot contrast    = '',f12.3)') Tjump_pot_nondim
     write(irefwr,'('' Ra_orig          = '',e12.3)') Ra_orig
     write(irefwr,'('' Ra reset to      : '',2e12.3)') Ra
endif

if (Ra==0) then
   DiRa=0
else
   DiRa=Di/Ra
endif

if (print_node) write(irefwr,*) 'call to set_up_tracers'
if (cyl) call detrminmax(kmesh1,kprob1,isol1)
call set_up_tracers()
inquire(file='postprocessing',exist=postprocessing)
if (postprocessing.and.print_node) then
   call comprP_post
   write(irefwr,*) 'stopping after option comprP_post'
   call instop
endif

ichoice=2
call initialize_temp_vel(ichoice)
!write(irefwr,*) 'after initialize_temp_vel'
!write(irefwr,*) '1: iuser_here = ',iuser_here(1)

!if (cartplume) then
   !if (print_node) write(irefwr,*) 'stopping in comprP_start after initialize_temp_vel'
   !call instop
!endif



! Create vector for viscosity after temperature has been created
inputcr(1)=7
inputcr(2)=1
inputcr(3)=2
inputcr(4)=1
inputcr(5)=1
inputcr(6)=0
inputcr(7)=-1
if (ivn.or.tackley) then
!  create vector containing the strain rate 
   if (print_node) write(irefwr,*) 'Create strain rate vector'
   call creatv(kmesh1,kprob1,isecinv,inputcr,rinputcr)
endif
!     *** create vector containing the viscosity
if (itypv/=0) then
  if (print_node) write(irefwr,*) 'Create viscosity vector(s)'
  call creatv(kmesh1,kprob1,ivisc,inputcr,rinputcr)
endif
if (tosi15 .or. tackley) then
   ! also create arrays for linear and plastic components to viscosity
   call creatv(kmesh1,kprob1,ivisclin,inputcr,rinputcr)
   call creatv(kmesh1,kprob1,iviscplast,inputcr,rinputcr)
endif



if (gable_plates.and.itracoption==1) then
   call set_up_chempost_logs
   if (initialize_crust) call makecrust(-1)
endif

if (itypv/=0) then
   if (print_node) write(irefwr,*) 'Fill by call to coefvis'
   call coefvis()
endif
i1=1
call algebr(6,1,isol2,i1,i1,kmesh1,kprob1,tmin2,tmax2,p,q,ipoint)
if (print_node) write(irefwr,*) 'temperature min/max: ',tmin2,tmax2
if (itypv>0) then
   call algebr(6,1,ivisc,i1,i1,kmesh1,kprob1,tmin2,tmax2,p,q,ipoint)
   if (print_node) write(irefwr,*) '  Viscosity min/max: ',tmin2,tmax2
endif

if (.not.stokes_is_updated.and.itracoption==0) then
   call stokesP(0)
   if (compress .and. .not.tala) then
      !call pevrms(vrms,isol1)
      !write(irefwr,*) 'vrms: ',vrms
      call stokesP(0) ! repeat to incorporate pressure effect in buoyancy
      !call pevrms(vrms,isol1)
      !write(irefwr,*) 'vrms: ',vrms
      !call stokesP(0) ! repeat to incorporate pressure effect in buoyancy
      !call pevrms(vrms,isol1)
      !write(irefwr,*) 'vrms: ',vrms
      !call stokesP(0) ! repeat to incorporate pressure effect in buoyancy
      !call pevrms(vrms,isol1)
      !write(irefwr,*) 'vrms: ',vrms
      !call instop
   endif
   call pevrms(vrms,isol1)
   vrms_d=vrms*rkappa_dim/height_dim*100*year_dim
   if (print_node) then
       write(irefwr,'(''Compute consistent Stokes solution IT1'')')
       write(irefwr,'('' vrms/vrms_d: '',2f12.3)') vrms,vrms_d
       write(irefwr,*) 'done with Stokes'
   endif
    if (print_node) write(irefwr,*) 'stop in comprP_start  after first stokes'
    call instop
endif

if (itracoption /= 0 .and. (.not.hannah)) then
   if (print_node) write(irefwr,*) 'vrms before stokes: ',vrms
   ibuoy_trac=1
   !call pevrms(vrms,isol1)
   !if (print_node) write(irefwr,'(''Before consistent Stokes solution with tracers CS1, vrms,cpu = '',e12.5,f12.2)') vrms,timeB-timeA
   call cpu_time(timeA)
   call stokesP(ibuoy_trac)
   call cpu_time(timeB)
   !  call pevrms(vrms,isol1)
   !  write(irefwr,*) 'vrms: ',vrms
   !  call stokesP(ibuoy_trac) ! repeat to incorporate pressure effect in buoyancy
   !  call pevrms(vrms,isol1)
   !  write(irefwr,*) 'vrms: ',vrms
   !  call stokesP(ibuoy_trac) ! repeat to incorporate pressure effect in buoyancy
   !  call pevrms(vrms,isol1)
   !  write(irefwr,*) 'vrms: ',vrms
   !  call stokesP(ibuoy_trac) ! repeat to incorporate pressure effect in buoyancy
   !  call pevrms(vrms,isol1)
   !  write(irefwr,*) 'vrms: ',vrms
   call pevrms(vrms,isol1)
   if (print_node) write(irefwr,'(''Compute consistent Stokes solution with tracers CS1, vrms,cpu = '',e12.5,f12.2)') &
      & vrms,timeB-timeA
   !if (print_node) write(irefwr,*) 'stop after CS1'
   !call instop
   if (use_tracers_from_nate) then
      if (print_node) write(irefwr,*) 'PWARN(comprP_start): stopping because use_tracers_from_nate=',use_tracers_from_nate
      call instop
   endif
   !call cpu_time(timeA)
   !call stokesP(ibuoy_trac)
   !call cpu_time(timeB)
   !call pevrms(vrms,isol1)
   !if (print_node) write(irefwr,'(''Before consistent Stokes solution with tracers CS1, vrms,cpu = '',e12.5,f12.2)') vrms,timeB-timeA
   !call bmout(iuser_here,user_here)
    !write(irefwr,*) 'stopping for test after CS1'
    !call instop
   ! determine velocity in tracers
   call tracvel(kmesh1,kprob1,isol1,user_here)
endif


if (stokes_only) then
   call bmout()
   !call check_stress_tensor(iuser_here,user_here)
   if (print_node) write(irefwr,*) 'stop in comprS because stokes_only=',stokes_only
   call instop
endif

if (.not.steady.and.start_from_steady) then
  ! find steady state solution at present Ra etc.
  ! can be used for benchmarks similar to JGR97
  Rb_here=Rb_local
  Rb_local=0
  if (print_node) write(irefwr,*) 'steady_iteration'
  call steady_iteration ! includes call to bmout
  Rb_local=Rb_here
endif

inout=0
noGMT_store=noGMT
noGMT=.true.
!noGMT=.false.
!petest=.true.
if (print_node) call peplotit
!noGMT=noGMT_store
if (petest) then
    if (print_node) write(irefwr,*) 'stop after peplotit',steady,start_from_steady
    call instop
endif
noGMT=noGMT_store

call cpu_time(cpu_after_start)
!write(irefwr,*) 'Gamma: ',Grueneisen
!write(irefwr,*) 'density: ',rho_n(0.0_8),rho_n(1.0_8)
!call instop


end subroutine comprP_start

subroutine setsepranmodulevalues()  ! adjust a few default values 
use sepmoduleplot
implicit none

jkader=-1
jmark=5
irotat=1 ! no rotation is allowed

end subroutine setsepranmodulevalues

subroutine get_checkmem_command
use control
#ifdef IFORT
use IFPORT
#endif
implicit none
integer :: iret
#ifdef MPI
return
#endif
!write(checkmem_command,'(''echo 0 $(awk '',a1,''/Rss/ {print "+", $2}'',a1,'' /proc/'',I0,''/smaps) | bc '')') "'","'",getpid()
!iret=system(checkmem_command)
write(checkmem_command,'(''echo 0 $(awk '',a1,''/Rss/ {print "+", $2}'',a1,'' /proc/'',I0,''/smaps) | bc > resmem'')') &
  &  "'","'",getpid()
write(6,'(a120)') checkmem_command
!iret=system(checkmem_command)
!iret=system("ls -l")

end subroutine get_checkmem_command

real(kind=8) function get_resmem()
#ifdef IFORT
use IFPORT
#endif
use control
implicit none
integer :: iret,isize
#ifdef MPI
return
#endif
iret=system(checkmem_command)
open(9,file='resmem')
read(9,*) isize
close(9)
get_resmem=isize*1.0_8/1024/1024
end function get_resmem
