! CREATE_ADDITIONAL_VECTORS
! 
! Set up vectors for density, expansivity etc.
! PvK 070511
subroutine create_additional_vectors()
use sepmoduleoldrouts
use coeff
use convparam
use sepran_arrays
use sepmodulecomio
use geometry
use control
implicit none
integer :: iu1(10),ipoint,i1,lu,i,ival
real(kind=8) ::  u1(10),tmin2,tmax2,p,q,avrho(6),avalpha(6)
real(kind=8) ::  avcond(6),avcp(6),rcond1,rcond2,rho1,rho2,funccf
real(kind=8) ::  Ks1,Ks2,Ta1,Ta2,cp1,cp2,alpha1,alpha2
real(kind=8) ::  pefvis,eta_top,eta_bot,rhos1,rhos2

! Create density vector (even if rho=1) and fill with analytical values
! With two problems on a single mesh choose these to become a solution 
iu1=0
! vector for temperature (otherwise errors appear in deriv and volint)
if (compress) then
!  use func(3)=funccf(3)
   iu1(1) = 3
else
!  density is constant
   iu1(1) = 0 
    u1(1) = 1d0
endif
call creavc(0,1001,1,irho,kmesh1,kprob1,iu1,u1,iu1,u1)


iu1(1) = 7
 u1(1) = 0d0
!write(6,*) 'PINFO(create_additional_vectors: ',iadiabat,Ts_eos0
call creavc(0,1001,1,iadia,kmesh1,kprob1,iu1,u1,iu1,u1)
iu1(1) = 0
call creavc(0,1001,1,ipress,kmesh1,kprob1,iu1,u1,iu1,u1)
iu1(1)=0 
 u1(1)=0d0
! ichcrv = ichcr + (iprob-1)*1000
call creavc(0,1001,1,icompwork,kmesh1,kprob1,iu1,u1,iu1,u1)

! Find reference values for density and overwrite rho with EOS if necessary
! For EBA and ALA we use (with eos_type=0) a functional description from Jarvis&McKenzie, 1980
! which allows for analytical estimates. 
if (Di > 0d0) then
    if (compress) then
      call algebr(6,1,irho,i,i,kmesh1,kprob1,rhos1,rhos2,p,q,i)
      rho1 = funccf(3,0d0,r1,0d0)
      rho2 = funccf(3,0d0,r2,0d0)
      rho_z0 = rho2
      rhomax = max(rho1,rho2)
      if (print_node) then
         write(irefwr,'(''            density at top and bottom: '',4e15.7)') rho1,rho2,rhos1,rhos2
         write(irefwr,'(''Dimensional density at top and bottom: '',2f12.1)') rho_dim*rho1,rho_dim*rho2
         write(irefwr,'(''Corresponding to r1,r2               : '',2f12.4)') r1,r2
      endif
    endif !compress
endif ! Di>0

!write(irefwr,*) 'stopping after <rho>'
!call instop

! iphi is used to find viscous dissipation in getvisdip()
call creavc(0,1,1,iphi,kmesh1,kprob1,iu1,u1,iu1,u1)
call creavc(0,1001,1,iphi2,kmesh1,kprob1,iu1,u1,iu1,u1)

call creavc(0,1001,1,iheat,kmesh1,kprob1,iu1,u1,iu1,u1)


! Find reference values for density and overwrite rho with EOS if necessary
! For EBA and ALA we use (with eos_type=0) a functional description from Jarvis&McKenzie, 1980
! which allows for analytical estimates. 

! Create expansivity vector (even if alpha=1) and fill with analytical values
if (ialphatype > 0) then
!  use func(105)=funccf(5)
   iu1(1) = 105
else
!  expansivity is constant
   iu1(1) = 0 
    u1(1) = 1d0
endif
call creavc(0,1001,1,ialpha,kmesh1,kprob1,iu1,u1,iu1,u1)

!  Create conductivity vector (even if cond=1) and fill with analytical values
!  Note that we need to have rho_av at this stage for the functional description in funccf
if (icondtype > 0.and.print_node) then
   if (rho_av < 1.0e-3) then
      write(irefwr,*) 'PERROR(compr_start): icondtype > 0 but rho_av=0'
      write(irefwr,*) 'icondtype, rho_av = ',icondtype,rho_av
   endif
!  use func(104)=funccf(4)
   iu1(1) = 104
else
!  conductivity is constant
   iu1(1) = 0 
    u1(1) = 1d0
endif
call creavc(0,1001,1,icond,kmesh1,kprob1,iu1,u1,iu1,u1)

iu1(1) = 0
 u1(1) = 1d0
call creavc(0,1001,1,icp,kmesh1,kprob1,iu1,u1,iu1,u1)
call creavc(0,1001,1,ibulkmod,kmesh1,kprob1,iu1,u1,iu1,u1)
!write(6,*) 'icp = ',icp

if (compress .and. eos_type/=0) then
!  Overwrite existing vectors that contain density, expansivity, conductivity
!  Generate vectors with specific heat and bulk modulus.
!  Note that when reading from tables (e.g., sfo05) the parameters
!  have been non-dimensionalized using values in convparam_mod
   !write(irefwr,*) 'irho: '
   call init_eos(1,irho,kmesh1,kprob1)
!  overwrite analytical estimate of average density with real value
   rho_av = eos_rhoav
   call init_eos(2,ialpha,kmesh1,kprob1)
   call init_eos(3,icp,kmesh1,kprob1)
   call init_eos(4,ibulkmod,kmesh1,kprob1)
   ! PvK 110420 updating icond via look up tables requires update of look up tables...
   ! call init_eos(5,icond,kmesh1,kprob1)
   if (print_node) write(irefwr,*) 'PINFO: eos_T(1) = ',eos_T(1)
   call init_eos(6,iadia,kmesh1,kprob1)
endif
!write(6,*) 'icp=',icp
!call instop

if (Di > 0 .and. print_node) then
!  Find min/max/average of temperature
   call algebr(6,1,iadia,i,i,kmesh1,kprob1,Ta1,Ta2,p,q,i)
    do lu=irefwr,irefwr
     write(lu,'(''Adiabatic temperature is defined: '')') 
     write(lu,'(''Non-dimensional range of values: '',2f12.3)') Ta1,Ta2
     write(lu,'(''Dimensional range of values    : '',2e12.5)') Ts_dimK+Ta1*deltaT_dim,Ts_dimK+Ta2*deltaT_dim
    enddo
endif
if (icondtype /= 0 .and. print_node) then
!  Conductivity is stored in icond
   call algebr(6,1,icond,i,i,kmesh1,kprob1,rcond1,rcond2,p,q,i)
    do lu=irefwr,irefwr
     write(lu,*)
     write(lu,'(''Conductivity is variable. '')') 
     write(lu,'(''Non-dimensional range of values: '',2f12.3)') rcond1,rcond2
     write(lu,'(''Dimensional range of values    : '',2e12.5)') rcond1*rk_dim,rcond2*rk_dim
    enddo
endif ! icondtype /= 0
if (ialphatype /= 0 .and. print_node) then
   call algebr(6,1,ialpha,i,i,kmesh1,kprob1,alpha1,alpha2,p,q,i)
    do lu=irefwr,irefwr
     write(lu,*)
     write(lu,'(''Expansivity is variable. '')')
     write(lu,'(''Non-dimensional range of values: '',2f15.3)') alpha1,alpha2
     write(lu,'(''Dimensional range of values    : '',2e15.5)') alpha1*alpha_dim,alpha2*alpha_dim
    enddo
endif ! ialphatype /= 0

if (compress .and. print_node) then
!     *** some output specifically for compressible convection
     call algebr(6,1,icp,i,i,kmesh1,kprob1,cp1,cp2,p,q,i)
     call algebr(6,1,ibulkmod,i,i,kmesh1,kprob1,Ks1,Ks2,p,q,i)
     call algebr(6,1,irho,i,i,kmesh1,kprob1,rho1,rho2,p,q,i)
     rhomax=max(rho1,rho2)
      do lu=irefwr,irefwr
       write(lu,*)
       write(lu,'(''Specific heat is variable '')')
       write(lu,'(''Non-dimensional range of values: '',2f15.3)') cp1,cp2
       write(lu,'(''Dimensional range of values    : '',2e15.5)') cp1*cp_dim,cp2*cp_dim
       write(lu,*)
       write(lu,'(''Bulk modulus is variable '')')
       write(lu,'(''Non-dimensional range of values: '',2f15.3)') Ks1,Ks2
       write(lu,'(''Dimensional range of values    : '',2e15.5)') Ks1*Ks_dim,Ks2*Ks_dim

       write(lu,*)
       write(lu,'(''Density is specified through irho'')')
       write(lu,'(''Non-dimensional range of values: '',2f15.3)') rho1,rho2
       write(lu,'(''Dimensional range of values    : '',2e15.5)') rho1*rho_dim,rho2*rho_dim
      enddo
endif


!if (itypv.ne.0) then
!    do lu=irefwr,irefwr
!      write(lu,*)
!      write(lu,*) 'Viscosity is variable, itypv = ',itypv
!      if (itypv == 2 .or. itypv == 3) then
!          eta_bot=pefvis(0d0,r1,T_bot,T_bot_adia,ival,0d0)
!          eta_top=pefvis(0d0,r2,T_top,T_top_adia,ival,0d0)
!        write(lu,'(''For itypv = '',i2,'': '')') itypv
!        write(lu,'(''Non-dimensional range of values: '',2f12.3)') eta_bot,eta_top
!      endif
!    enddo
!endif

end subroutine create_additional_vectors
