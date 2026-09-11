subroutine setup_temperature_bc(iuserh,userh)
use sepmoduleoldrouts
use sepmodulecomio
use convparam
use control
use coeff
use geometry
use sepran_arrays
implicit none
integer :: iuserh(*)
real(kind=8) :: userh(*)
integer :: jbench_type

if (delta1K) then
   if (compress) then
      if (Kequivalent) then
         T_top=Ts_nondimK
         T_bot=Ts_nondimK+T_bot*deltaT_dim
      else
         T_top=0.0_8
         T_bot=T_bot*deltaT_dim
      endif
   else ! if (EBA.or.BA) then
      T_top=0.0_8
      T_bot=T_bot*deltaT_dim
      if (Kequivalent) then
         ! pretty much a patch for CY85 with weird Ts
         T_top=273.0_8 
         T_bot=273.0_8+T_bot
      endif
   endif
else
   T_top=0.0_8
endif


! Feb 2021: this has become obsolete with Tbar==0
!jbench_type=ibench_type-(ibench_type/10)*10
!if (compress .and. jbench_type/=3 .and. jbench_type/=7 .and. (jbench_type/=5.and.Di>1e-7).and..not.cartplume) then
!!  Update T_bot to reflect adiabatic temperature increase
!!  Ignore this for the King et al. GJI benchmark (ibench_type==3)
!!  and the cylindrical benchmark (ibench_type==7)
!!  as well as the JGR97 style benchmark with compressibility
!   T_bot_orig = T_bot
!   T_top_orig = T_top
!   T_bot = T_bot + Ts1_nondim - Tso_nondim
!   if (print_node) then
!      write(irefwr,'('' T_bot is updated from '',f12.3, '' to '',5f12.3)') T_bot_orig,T_bot, &
!         & Ts1_nondim,Tso_nondim,deltaT_dim_TBL,deltaT_dim_LBL
!      write(irefwr,'('' T_top is updated from '',f12.3, '' to '',f12.3)') T_top_orig,T_top
!   endif
!endif

! ?? See above
!if (Kequivalent) then
!  T_top=T_top+273.0_8
!   T_bot=T_bot+273.0_8
!endif

if (cartplume) then
   ! overwrite for case without top boundary layer
   if (Kequivalent) then
      T_top=1600
      T_bot=1600*(-palpha*Di+1)**(-1.0_8/palpha)
   else if (delta1K) then
      T_top=1600.0_8-273.0_8
      T_bot=1600*(-palpha*Di+1)**(-1.0_8/palpha)+deltaT_dim_LBL-273.0_8
   else 
      if (print_node) then
         write(irefwr,*) 'PERROR(setup_temperature_bc): for cartplume use either delta1K or Kequivalent' 
         call instop
      endif
   endif
   if (print_node) then
      write(irefwr,'('' T_bot is updated to '',f12.3)') T_bot ! ,1600*(-palpha*Di+1),1600*(-palpha*Di+1)**(-1.0_8/palpha)
      write(irefwr,'('' T_top is updated to '',f12.3)') T_top
   endif
endif


if (print_node) then
  write(irefwr,*) 'Initial condition for temperature, including b.c.'
  write(irefwr,'(''PINFO(setup_temperature): T_top, T_bot= '',2f12.3)') T_top,T_bot
  write(irefwr,'(''                          deltaT_dim  = '',f12.3)') deltaT_dim
  write(irefwr,*)
endif

if (compress.and.solve_for_Tperturb) then
   T_bot=T_bot-Tbars_nondimK*(exp(Di)-1)
endif


if (.not.insulbot) then
   call bvalue(0,1002,kmesh1,kprob1,isol2,T_bot,ibottom,ibottom,1,0)
   if (print_node) then
     write(irefwr,*) 'set temperature at curve ',ibottom,' to ',T_bot
   endif
endif
call bvalue(0,1002,kmesh1,kprob1,isol2,T_top,itop,itop,1,0)
if (print_node) then
  write(irefwr,*) 'set temperature at curve ',itop,' (top) to ',T_top
endif

end subroutine setup_temperature_bc
