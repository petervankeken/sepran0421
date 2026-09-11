! ichoice: 1 = initialize T
!          2 = initialize T+V
subroutine initialize_temp_vel(ichoice)
use sepmoduleoldrouts
use sepmodulecomio
use sepmodulecpack
use sepran_arrays
use sepran_interface
use convparam
use control
use coeff
use brandenburg
use geometry
use tracers
implicit none
integer,intent(in) :: ichoice
real(kind=8) :: rinvec(2)
integer :: ichvc,idum,iprob_here,ipert,iinvec(11),inpcre(20)
real(kind=8) :: tmin_pert,tmax_pert,p,q,contln(11),vrms,u1(5),rincre(20)
real(kind=8) :: q_a(10),veloc(10),veloc_d(10),gnus1,gnus2,temp1,temp00
real(kind=8) :: tmin2,tmax2,vrms_d,tmin3,tmax3,format,T_bot_here,T_top_here,Rb_here
integer :: i1,ipoint,ncntln=0,np,nt,i,icurvs(5),ic,iu1(5)
logical :: solve_for_conductive
real :: timeB,timeA





if (gable_plates) then
   if (nplates<=0) then
      if (print_node) write(irefwr,*) 'PERROR: gable_plates is on but nplates <=0: ',nplates
      call instop
   endif
endif

if (pedebug.and.print_node) write(irefwr,*) 'irestart: ',irestart
if (irestart < 0) then
   ! initialize with conductive solution (irestart==-1) 
   ! For compress with absolute T
   ! use potential temperature contrast to set up conductive 
   ! solution for T'. Tbar is added below
   if (pedebug.and.print_node) write(irefwr,*) 'before bvalue isol1,isol2: ',isol1,isol2
! Feb2021 This has become obsolete: Tbot for solve_for_Tperturb has already been set to the reduced value
!   if (compress.and.solve_for_Tperturb) then
!      T_bot_here = T_bot - Ts1_nondim + Tso_nondim
!      T_top_here = T_top ! - Tso_nondim
!   else
!      T_top_here = T_top
!      T_bot_here = T_bot
!   endif
   if (.not.insulbot) then
      !call bvalue(0,1002,kmesh1,kprob1,isol2,T_bot_here, ibottom,ibottom,1,0)
      call bvalue(0,1002,kmesh1,kprob1,isol2,T_bot, ibottom,ibottom,1,0)
   endif
   if (print_node) write(irefwr,*) 'sset temperature at curve ',ibottom,' (bot) to ',T_bot
   !Feb2021 call bvalue(0,1002,kmesh1,kprob1,isol2, T_top_here,itop,itop,1,0)
   call bvalue(0,1002,kmesh1,kprob1,isol2, T_top,itop,itop,1,0)
   if (print_node) write(irefwr,*) 'set temperature at curve ',itop, ' (top) to ',T_top
   if (irestart == -1) then
     ! conductive solution + perturbation
     solve_for_conductive=.true.
     if (debugTbars) solve_for_conductive=.false.
     if (solve_for_conductive) then
        ! find conductive solution by solving heat equation without advective terms
        conductive=.true.
        if (print_node) write(irefwr,*) 'solve for conductive solution'
        call steady_heatP()
        conductive=.false.
     else 
       ! create initial condition using func() with conductive profile
       iprob_here=2  ! second problem is temperature
       ichvc=1  ! type of solution vector
       iu1lc(1)=12  ! set dof 1 to func(12,x,y,z)
       call creavc(0,ichvc+(iprob_here-1)*1000,idum,isol2,kmesh1,kprob1,iu1lc,u1lc)
     endif
     stokes_is_updated=.false.
     if (print_node) then
       ncntln=0
       write(ourplotname,'(''PLOTS/COND.0000'')') 
       call plotc1(1,kmesh1,kprob1,isol2,contln,ncntln,15.0_8,1.0_8,1)
        call algebr(6,1,isol2,i1,i1,kmesh1,kprob1,tmin_pert,tmax_pert,p,q,ipoint)
        write(irefwr,*) 'icond: ',tmin_pert,tmax_pert
     endif

     ! create perturbation
     iprob_here=2  ! second problem is temperature
     ichvc=1  ! type of solution vector
     iu1lc(1)=1  ! set dof 1 to func(2,x,y,z)
     call creavc(0,ichvc+(iprob_here-1)*1000,idum,ipert,kmesh1,kprob1,iu1lc,u1lc)
     if (print_node) then
        call algebr(6,1,ipert,i1,i1,kmesh1,kprob1,tmin_pert,tmax_pert,p,q,ipoint)
        write(irefwr,*) 'ipert: ',tmin_pert,tmax_pert
        ncntln=0
        write(ourplotname,'(''PLOTS/PERT.0000'')') 
        call plotc1(1,kmesh1,kprob1,ipert,contln,ncntln,15.0_8,1.0_8,1)
     endif
     ! add perturbation to conductive solution
     iinvec=0
     iinvec(1)=11
     iinvec(2)=27
     iinvec(11)=2 ! iprob
     rinvec(1)=1.0_8
     rinvec(2)=1.0_8
     call manvec(iinvec,rinvec,isol2,ipert,isol2,kmesh1,kprob1)
     if (print_node) then
        call algebr(6,1,isol2,i1,i1,kmesh1,kprob1,tmin_pert,tmax_pert,p,q,ipoint)
        write(irefwr,*) 'isol2+ipert: ',tmin_pert,tmax_pert
        write(ourplotname,'(''PLOTS/COND_PERT.0000'')') 
        call plotc1(1,kmesh1,kprob1,isol2,contln,ncntln,15.0_8,1.0_8,1)
     endif
     stokes_is_updated=.false.
  else if (irestart<=-2) then
!    Compute initial condition through func(-irestart).
     if (print_node) then
        write(irefwr,'(''irestart<=-2 compute initial condition '',$)')
        write(irefwr,'(''from func with ichoice='',i3)') -irestart
     endif
!    length_inpre, no_vec, type_of_vec, iprob, ivec, icompl, ifill, 
!    Total number of entries
     inpcre(1)=10
!    Number of vectors to be created
     inpcre(2)=1
!    Type of vector  (0 or 1=solution vector)
     inpcre(3)=0
!    problem number
     inpcre(4)=2
!    sequence number of array of special structure
     inpcre(5)=0
!    0=real vector
     inpcre(6)=0
!    IFILL: filling is defined by IFILL commands
     inpcre(7)=1
!    fill first degree of freedom
     inpcre(8)=0
!    fill using FUNC with ichoice=-irestart
     inpcre(9)=-irestart
!    fill all nodes
     inpcre(10)=0
     call creatv(kmesh1,kprob1,isol2,inpcre,rincre) 
     call algebr(6,1,isol2,i1,i1,kmesh1,kprob1,tmin2,tmax2,p,q,ipoint)
     if (print_node) then
        write(irefwr,'(''minmax(isol2) after creatv: '',2f12.3)') tmin2,tmax2
     endif
     !call instop
     stokes_is_updated=.false.
  endif ! irestart<=-2

  if (compress.and. (.not.solve_for_Tperturb)) then
        ! add adiabat
        if (print_node) write(irefwr,*) 'add adiabat and solve_for_Tperturb=',solve_for_Tperturb
        iinvec(1)=2
        iinvec(2)=27
        rinvec(1)=1d0
        rinvec(2)=1d0
        if (print_node) then
           call algebr(6,1,isol2,i1,i1,kmesh1,kprob1,tmin2,tmax2,p,q,ipoint)
           write(irefwr,'(''minmax(isol2) before addition: '',4f12.3)') &
                & tmin2,tmax2,tmin2*deltaT_dim,tmax2*deltaT_dim
           call manvec(iinvec,rinvec,isol2,iadia,isol2,kmesh1,kprob1)
           write(irefwr,*) 'bottom: ',T_bot,ibottom
           write(irefwr,*) 'top   : ',T_top,itop
           if (.not.insulbot) then
             call bvalue(0,1002,kmesh1,kprob1,isol2,T_bot,ibottom,ibottom,1,0)
           endif
           call bvalue(0,1002,kmesh1,kprob1,isol2,T_top,itop,itop,1,0)
           call algebr(6,1,isol2,i1,i1,kmesh1,kprob1,tmin2,tmax2,p,q,ipoint)
           call algebr(6,1,iadia,i1,i1,kmesh1,kprob1,tmin3,tmax3,p,q,ipoint)
           if (print_node) then
             write(irefwr,'(''minmax(iadia)               : '',4f12.3)') tmin3,tmax3,tmin3*deltaT_dim,tmax3*deltaT_dim
             write(irefwr,'(''minmax(isol2) after addition: '',4f12.3)') tmin2,tmax2,tmin2*deltaT_dim,tmax2*deltaT_dim
           endif
           write(ourplotname,'(''PLOTS/COND_PERT_ADIA.0000'')') 
           call plotc1(1,kmesh1,kprob1,isol2,contln,ncntln,15.0_8,1.0_8,1)
        endif
     endif
     call nusseltP(q_a)
      if (print_node) then
         write(irefwr,'(''q_a              : '',3f10.3)') q_a(1),q_a(2),q_a(3)
         write(ourplotname,'(''PLOTS/T_INIT.001'')')
         format=10.0_8
         call plotc1(1,kmesh1,kprob1,isol2,contln,ncntln,format,1.0_8,1)
         call averages(user_here)
      endif
      !write(irefwr,*) 'stopping after adiabat addition' 
      !call instop

else ! irestart >= 0 
   if (print_node) write(irefwr,*) 'Initial condition for T  is read from: ', Tstartfile(1:len_trim(Tstartfile))
   !call sepactsolbf1(isol2)
   call readbs_netcdf(Tstartfile,isol2)
   if (read_velocity) then
      if (print_node) write(irefwr,*) 'Initial condition for UV is read from: ', UVstartfile(1:len_trim(UVstartfile))
      !call sepactsolbf1(isol1)
      call readbs_netcdf(UVstartfile,isol1)
      stokes_is_updated = .true.
   else
      if (print_node) then
         write(irefwr,*) 'Velocity is not read from file'
         write(irefwr,*) 'compute consistent velocity'
      endif
      stokes_is_updated=.false.
   endif
endif

! Let's make sure the boundary conditions are updated for all cases
if (.not.insulbot) then
   call bvalue(0,1002,kmesh1,kprob1,isol2,T_bot,ibottom,ibottom,1,0)
endif
call bvalue(0,1002,kmesh1,kprob1,isol2,T_top,itop,itop,1,0)

call copyvc(isol1,isolold1)
call copyvc(isol1,isolold1(2))
call copyvc(isol1,isolold1(3))
call copyvc(isol2,isolold2)


! Initial estimate (without chemical buoyancy)
if (.not.stokes_is_updated.and.itracoption==0) then
!     ! compute consistent Stokes solution
!   !  write(irefwr,*) '3: iuser_here = ',iuser_here

      call cpu_time(timeA)
      if (print_node) write(irefwr,*) 'call stokesP: ',do_not_do_Stokes
      if (.not.do_not_do_Stokes) then
          call stokesP(1)
      endif
      call cpu_time(timeB)
      if (print_node) write(irefwr,*) 'first stokes solve cpu=',timeB-timeA
      if (debugTbars) then
         call bmout()
         call instop
      endif  
      if (compress.and. .not.TALA) then
         ! repeat to incorporate effect of pressure
         if (.not.do_not_do_Stokes) call stokesP(1)
      endif
!      call cpu_time(timeB)
!      if (pedebug) write(irefwr,*) 'consistent Stokes solve took: ',timeB-timeA
!   !  write(irefwr,*) '4: iuser_here = ',iuser_here
      call pevrms(vrms,isol1)
!   !  write(irefwr,*) '5: iuser_here = ',iuser_here
      stokes_is_updated=.true.
      if (print_node) write(irefwr,*) 'vrms = ',vrms
      !if (print_node) then
      !   write(irefwr,*) 'stop after stokes in initialize_temp_vel'
      !endif
      !call instop
  
endif
!if (stokes_only) then
!   if (print_node) write(irefwr,*) 'stop in initialize_tempvel because stokes_only=',stokes_only
!   call bmout()
!   call check_stress_tensor(iuser_here,user_here)
!   write(irefwr,*) 'do_not_do_Stokes: ',do_not_do_Stokes
!   call instop
!endif

!write(irefwr,*) 'pvk stop in initialize_temp_vel 230'
!call instop


call nusseltP(q_a)
temp1=q_a(3)
temp00=q_a(4)
gnus1=q_a(1)/(temp1-temp00)  ! update for delta1K
gnus2=q_a(2)/(temp1-temp00)
call pevrms(vrms,isol1)
if (print_node) then
     write(irefwr,*)
     write(irefwr,'(''q_a              : '',4f10.3)') gnus1,gnus2,q_a(3),q_a(4)
     write(irefwr,'(''vrms             : '',f8.3)') vrms
     call surfacevel(kmesh1,kprob1,isol1,veloc,veloc_d,-1)
     write(irefwr,'(''surface velocity : '',2f15.7)') veloc(1),veloc(2)
endif




call copyvc(isol1,isolold1)
call copyvc(isol1,isolold1(2))
call copyvc(isol1,isolold1(3))
call copyvc(isol2,isolold2)

if (gable_plates) then
   iu1=0
    u1=0.0_8
   do i=1,nplates
      call copyvc(isol1,plate_sol(i))
      if (print_node) write(irefwr,'(''Create Gable solution '',2i5)') i,plate_sol(1)
!     write(irefwr,*) 'plate_sol: ',i,plate_sol(i)
   enddo
   call creavc(0,1,1,irotation,kmesh1,kprob1,iu1,u1,iu1lc,u1lc)
   if (print_node) write(irefwr,*) 'irotation: ',irotation
   call force_balance_setup()
   if (print_node) then
     do i=1,nplates+1
        write(irefwr,'(''plate boundary '',i5,'' is at '',f15.7)') i,plate_boundaries(i)
     enddo
   endif
endif

if (print_node) then
   ncntln=0
   write(ourplotname,'(''PLOTS/T_INIT.0000'')') 
   call plotc1(1,kmesh1,kprob1,isol2,contln,ncntln,15.0_8,1.0_8,1)
endif


end subroutine initialize_temp_vel
