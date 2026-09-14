! SET_UP_EOS
! 
! Error correction February 2021:
! Converged solution should be invariant to choice of Tbar.
! Hence a chance in dynamics by switching from Tbars_dim=273 K to Tbars_dim=1600 K is an artefact.
! When solving for Tpot the difference Tbars_nondimK-Ts_nondimK enters the work term. 
! When solving for full T it is better to set Tbar=0.
!
! First pass in setting up the EoS parameters.
! For compressible convection (ALA) this includes tables with
! dimensional quantities (eos_alpha_d, etc.) and their
! nondimensional equivalents (eos_alpha, etc.). 
! 
! Figure out the EoS. 
! eos_type = 0  is WA with constant (King2010) 
! We'll only allow solve_for_Tperturb (i.e., for solving T' for the King2010 benchmark).
! Ohterwise we solve for full temperature with Tbar==0.

! PvK October 2011
! Updated with Kequivalent PvK February 2021
subroutine set_up_eos
use sepmodulecomio
use convparam
use coeff
use mparallel
use control
use geometry
use eos
implicit none
real*8 rho_max,rho1,rho2,z,rl,height,dz
real(kind=8) :: adiabat
real(kind=8) :: funccf,y,alphal,cpl,rhol
real(kind=8) :: dzero=0.0_8
integer :: lu,find_eos_iz,iz,i,iph,jbench_type
character(len=80) :: eosname

rho_av = 1d0
rho_z0 = 1d0
rhomax = 1d0

if (eos_type /= 0) then
   if (print_node) then
     write(irefwr,*) 'PERROR(set_up_eos): use eos_type=0 for now'
     write(irefwr,*) 'code needs to be updated for other choices'
   endif
   !call instop
endif



compress=.false.
BA=.false.
EBA=.false.
write(irefwr,*) 'set_up_eos: ',pemcont,Di
! Check for EBA, TALA, or ALA
mcontv=10*pemcont
if (pemcont==1) then
   if (print_node) then
      write(irefwr,*) 'Use compressible formulation. TALA = ',TALA
   endif
   EBA=.false.
   compress=.true.
   if (print_node) write(irefwr,*) 'Di = ',Di
else if (Di>0) then
   if (print_node) write(irefwr,*) 'Use EBA'
   EBA=.true.
   compress=.false.
endif
if (abs(Grueneisen)>1e-7) then
  DiG = Di/Grueneisen
else
  DiG = Di
endif

if (Di<1e-7.and.pemcont==0)  BA=.true.



write(irefwr,*) 'EBA: ',BA,EBA,TALA,compress


if (compress) then
     absolute1=.true.
endif ! compress

! Set temperature_in_C and Ts_nondimK, Tbars_nondimK
if (delta1K.and..not.Kequivalent) then
   temperature_in_C=.true.
else 
   temperature_in_C=.false.
endif
T_bc_top=0.0_8

! Note that Tbar=0 for full T moving forward from Feb 2021
Ts_nondimK=0.0_8
if (compress) then
   if (solve_for_Tperturb) then
      ! Force adiabatic surface T to 273
      if (delta1K) then
         Ts_nondimK=Ts_dimK
         Tbars_nondimK=Tbars_dimK
      else
         Ts_nondimK=Ts_dimK/deltaT_dim
         Tbars_nondimK=Tbars_dimK/deltaT_dim
      endif
   else  
      ! full temperature
      if (.not.Kequivalent) then
         if (delta1K) then
            Ts_nondimK=Ts_dimK
         else
            Ts_nondimK=Ts_dimK/deltaT_dim
         endif
      endif
   endif
endif
if (EBA.or.BA) then
   if (.not.Kequivalent) then
      if (delta1K) then
         Ts_nondimK=Ts_dimK
      else
         Ts_nondimK=Ts_dimK/deltaT_dim
      endif
   else
      Ts_nondimK=Ts_dimK-273.0 ! pretty much a patch for CY85 'weird' Ts
   endif
endif
dCY85_nondim=dCY85/deltaT_dimK
!write(6,*) 'dCY85: ',dCY85,dCY85_nondim+Ts_nondimK


! Figure out which approximation we're using
jbench_type=ibench_type-(ibench_type/10)*10
if (TALA .and. jbench_type/=3 .and. ibench_type/=17) then
!  Shouldn't use TALA since it is inaccurate but make an exception for the King et al. benchmark
   if (myid==0) then
     write(irefwr,*) 'PERROR(set_up_eos): TALA should not '
     write(irefwr,*) '   be used since it is a bit inaccurate'
   endif
!  call instop()
endif

! For EBA and ALA we use (with eos_type=0) a functional description from Jarvis&McKenzie, 1980
! which allows for analytical estimates. 
if (Di > 0d0) then
  if (eos_type == 0) then
!   Analytical function allows for analytical estimates of <rho>
    if (print_node) then
      do lu=irefwr,irefwr
      write(lu,*)
      write(lu,'(''Analytical density values for eos_type == 0'')')
      enddo
    endif
    if (cyl) then
      if (.not.axi) then
!        <rho> in cylindrical
         rho_av = (r2/DiG + 1d0/(DiG*DiG))*exp(-DiG*r2)
         rho_av = rho_av -(r1/DiG +1d0/(DiG*DiG))*exp(-DiG*r1)
         rho_av = -2*rho_av *exp(DiG*r2)/(r2*r2-r1*r1)
      else ! axi
!        <rho> in axisymmetric spherical
         rho_av = (1d0/DiG*r1*r1 + 2d0/(DiG*DiG)*r1 +2d0/(DiG*DiG*DiG))*exp(-DiG*r1) - &
                & (1d0/DiG*r2*r2 + 2d0/(DiG*DiG)*r2 +2d0/(DiG*DiG*DiG))*exp(-DiG*r2)
         rho_av = 3*exp(DiG*r2)/(r2*r2*r2-r1*r1*r1)*rho_av
      endif
    else ! .not.cyl
!           *** Cartesian geometry
      rho_av  = exp(DiG)/DiG - 1d0/DiG
    endif ! cyl

    if (print_node) then
      do lu=irefwr,irefwr
        write(lu,'(''Average normalized density (rho0=1) is '',f15.3)') rho_av
      enddo
    endif
  endif
endif

! Figure out the EoS. 
! eos_type = 0  is WA with constant (King2010) or variable alpha (BvK13)
! eos_type = 1  is sfo05
! eos_type = -1 is same as 0 but specified in table with Ts_eos0=273
! eos_type = -2 is for WA 1600*exp(Di*z/Gamma) (default)

if (compress) then
!    need to read info from file
!    read just coordinates when using eos_type=0
     if (eos_type == 1) then
        eos_data_base = 'sfo05'
     else if (eos_type == 2) then
!       version in which we use sfo05 for alpha,cp but 
!       independently set up rho  and adiabat
        eos_data_base = 'sfo05'
     else if (eos_type == -1 .or. eos_type==0) then
!       Jarvis 1980 with Ts=273 K
        eos_data_base = 'Jarvi'
     else if (eos_type == -2) then
!       Jarvis 1980 with Ts=1600 K
        eos_data_base = 'Jar16'
     else 
        write(irefwr,*) 'PERROR(compr_start): eos_type = ',eos_type
        write(irefwr,*) 'is undefined. 0=exp(Di*z/gamma), 1=sfo05'
        write(irefwr,*) '-1 = Jarvis80, -2 = Jarvis80 with Ts=1600 K'
        call instop()
     endif
!    first check if global data file is available
     write(eosname,'(''/work/pvankeken/legato/work2/comprS/'',a5,''.dat'')') eos_data_base
     inquire(file=eosname,exist=eos_exists)
     if (eos_exists) then
        open(9,file=eosname) 
     else
!       check existence of local file
        write(eosname,'(a5,''.dat'')') eos_data_base
        inquire(file=eosname,exist=eos_exists)
        if (eos_exists) then
           open(9,file=eosname) 
        else
           if (myid==0) then
           write(irefwr,*) 'PERROR(compr_start): eos_type = ',eos_type
           write(irefwr,*) 'but I cannot find correct eos_data file'
           write(irefwr,*) 'seeking for ',eos_data_base
           write(irefwr,*) 'locally or in /work/pvankeken/legato/work2/comprS'
           endif
           call instop
        endif! eos_exists - 1
     endif  ! eos_exists - 2
     if (myid==0) write(irefwr,*) 'read EoS information from ',eosname
     read(9,*) eos_np
     if (eos_np > EOS_NPMAX) then
        if (print_node) then
          write(irefwr,*) 'PERROR(compr_start): eos_np > EOS_NPMAX' 
          write(irefwr,*) 'eos_np, EOS_NPMAX: ',eos_np,EOS_NPMAX
          write(irefwr,*) 'from file : ',eosname
        endif
        call instop
     endif
     do i=1,eos_np
        read(9,*) eos_z_d(i),eos_p_d(i),eos_T_d(i),eos_alpha_d(i),eos_cp_d(i),eos_K_d(i),eos_rho_d(i)
     enddo
     close(9)

     ! convert eos_z_d to SI
     eos_z_d=eos_z_d*1e3_8
     ! non-dimensionalize
     eos_z=eos_z_d/eos_z_d(eos_np)

     if (ialphatype/=0.and.ialphatype/=1.and.ialphatype/=3) then
        if (print_node) then
          write(irefwr,*) 'PERROR(set_up_eos): ialphatype should be '
          write(irefwr,*) ' 0,1,3  but is : ',ialphatype
        endif
        call instop
     endif


     if (eos_type == 0) then
       ! overwrite table with entries using analytical expressions
       eos_rho_d(1)=rho_dim
       eos_T_d(1)=Tbars_nondimK
       eos_alpha_d(1)=alpha_dim
       eos_cp_d=1250.0_8 ! assumed constant throughout
       do i=2,eos_np
          eos_T_d(i) = get_Tbar_dim(eos_z_d(i))
          eos_rho_d(i) =get_rhobar_dim(eos_z_d(i))
          eos_alpha_d(i)= get_alpha_dim(eos_z_d(i),eos_rho_d(i))
       enddo
   
    endif
    eos_T = eos_T_d-Ts_dimK
    eos_z_d(1)=0d0
    eos_z_d(eos_np)=height_dim
    write(irefwr,'(''z_d    : '',6e15.7)')  eos_z_d(1:eos_np:eos_np/5)
    write(irefwr,'(''alpha_d: '',6e15.7)')  eos_alpha_d(1:eos_np:eos_np/5)
    write(irefwr,'(''rho_d  : '',6e15.7)')  eos_rho_d(1:eos_np:eos_np/5)
    write(irefwr,'(''T_d    : '',6e15.7)')  eos_T_d(1:eos_np:eos_np/5)
    write(irefwr,'(''T      : '',6e15.7)')  eos_T(1:eos_np:eos_np/5)
    write(irefwr,*) 'iadiabat: ',iadiabat


     if (ibench_type == 8) then
!      special case for Leng&Zhong 2010 redo; eos_T is set up in update_rho_adia
       ialphatype=3
       do i=1,eos_np
           z = eos_z_d(i)/height_dim
           alphal = get_alpha_dim(eos_z_d(i),eos_rho_d(i))/alpha_dim
           y = R2-z
           eos_alpha_d(i)=alphal*alpha_dim*4e0/3  !! LZ10 have 4e-5 as alpha0
           eos_rho_d(i)=exp(0.5*z)*rho_dim*3.4/3.3 !! LZ use rho0=3400; LZ10 eq 8.
       enddo
       eos_cp_d = cp_dim*1000.0_8/1250.0_8   !! LZ10 have Cp=1000
       write(irefwr,*) 'LENG_ZHONG: eos_alpha/rho: '
       write(irefwr,*) R2-eos_z_d(1)/height_dim,eos_alpha_d(1),eos_rho_d(1)
       write(irefwr,*) R2-eos_z_d(eos_np)/height_dim,eos_alpha_d(eos_np),eos_rho_d(eos_np)
     endif
   open(9,file='eos_local.dat')
   write(9,'(6a15)') 'y','zdim','rhodim','Tdim','alphadim','cond_dim'
   do i=1,eos_np
      y=R2-eos_z_d(i)/height_dim
      write(9,'(6e15.7)') y,eos_z_d(i),eos_rho_d(i),eos_T_d(i),eos_alpha_d(i),funccf(4,dzero,y,dzero)
   enddo
   close(9)
endif  ! compress 


if (compress) then
!    Convert selected quantities to non-dimensional quantities, scaled by surface
!    values specified in convparam_mod. Note that tables have
!    bulk modulus specified in bar
     do i=1,eos_np
        eos_alpha(i) = eos_alpha_d(i)/alpha_dim
        eos_cp(i) = eos_cp_d(i)/cp_dim
        if (eos_type > 0) then
           eos_K(i) = eos_K_d(i)*1.01e5/Ks_dim
        else 
           eos_K(i) = eos_K_d(i)/Ks_dim
        endif
        eos_rho(i) = eos_rho_d(i)/rho_dim
        if (ibench_type==4) then
           z=eos_z_d(i)/height_dim
           eos_rho(i) = funccf(3,dzero,1.0_8-z,z)
        endif
     enddo
     eos_T=0.0_8 ! force adiabatic / reference Tbar to 0
endif ! compress



if (print_node) then
   write(irefwr,'(''   After setting up EoS: '')')
   write(irefwr,'(''   Ts_dimK      =    '',2f12.3)') Ts_dimK,Ts_nondimK
   write(irefwr,'(''   deltaT_dim  =    '',f12.3)') deltaT_dim
   write(irefwr,'(''   iadiabat/eostype:'',2i12)') iadiabat,eos_type
endif

!write(6,*) 'stopping in set_up_eos'
!call instop

end subroutine set_up_eos

