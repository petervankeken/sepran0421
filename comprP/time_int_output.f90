subroutine time_int_output()
use sepmodulesol
use sepmodulecomio  ! provides lu's for reading & writing
use sepmodulecpack  ! parallel computing interface
use sepran_arrays ! access to kmesh, kprob, isol, etc. and interfaces to subroutines
use sepran_interface
use brandenburg
use mpetrac
use tracers
use coeff
use geometry
use convparam
use control
use mtime
use dtm_elem  ! controls building of pressure mass matrix
implicit none
integer :: ndim,nphys,npoint,nunkp
integer :: modelv,mconv,nparm,nsup1,idum,iread,ihelp,icurvs(2)
integer :: idegfd,i,ichoice,isubr,ncntln
real(kind=8) :: rayleigh,umax,outputtime,contln(10),gnus1,gnus2
character(len=120) :: namedof(4)
integer :: ichvc,iprob_here,ip,ipbuoy2,ipuser
real(kind=8) :: q_a(20),anorm,vrms,t_d,veloc(10),veloc_d(10),avphi(10),avwork(10),get_resmem
logical :: first=.true.
integer :: ih,im,is,imelt

if (.not.print_node) return


if (first .and. print_node) then
   if (restart) then
     open(LU_VRMS,file='vrms.dat',status='old',access='append')
     open(LU_NU,file='nusselt.dat',status='old',access='append')
     open(LU_NU_BOT,file='nusselt_bot.dat',status='old',access='append')
     open(LU_HEATFLOW,file='heatflow.dat',status='old',access='append')
     if (cyl) open(LU_ROTATION,file='rotation.dat',status='old',access='append')
     open(LU_SURFACEVEL,file='surfacevel.dat',status='old',access='append')
     open(LU_WORK,file='work.dat',status='old',access='append')
     open(LU_PHI,file='phi.dat',status='old',access='append')
     open(LU_MEM,file='mem.dat',status='old',access='append')
   else
     open(LU_VRMS,file='vrms.dat')
     open(LU_NU,file='nusselt.dat')
     open(LU_NU_BOT,file='nusselt_bot.dat')
     open(LU_HEATFLOW,file='heatflow.dat')
     if (cyl) open(LU_ROTATION,file='rotation.dat')
     open(LU_SURFACEVEL,file='surfacevel.dat')
     open(LU_WORK,file='work.dat')
     open(LU_PHI,file='phi.dat')
     open(LU_MEM,file='mem.dat')
   endif
endif
first=.false.

call nusseltP(q_a)
call pevrms(vrms,isol1)
call surfacevel(kmesh1,kprob1,isol1,veloc,veloc_d,-1)
gnus1=q_a(1)/q_a(3)
gnus2=q_a(2)/q_a(3)
call getvisdip
call averageT(2,iuser_here,user_here,avphi)
call compressionwork(kmesh1,kprob1,isol1,isol2,icompwork)
call averageT(3,iuser_here,user_here,avwork)


if (print_node) then
   cpu_then = cpu_now
   cpu_now=second()
   dcpu=cpu_now-cpu_then
   cpu_total=cpu_now-cpu_first
   ih=int(cpu_total/3600)
   im=int((cpu_total-ih*3600)/60)
   is=int(cpu_total-ih*3600-im*60)
   if (use_varRa) then
     if (print_node) write(irefwr,'(''time='',e15.7,2e12.5,4e12.5,f10.2,i5,'':'',i2.2,'':'',i2.2)')  &
      & time_now,Ra,gnus1,vrms,gamma(1),phz0(1),pht0(1),dcpu,ih,im,is
   else 
     if (print_node) write(irefwr,'(''time='',e15.7,e12.5,4e12.5,f10.2,i5,'':'',i2.2,'':'',i2.2)')  &
      & time_now,gnus1,vrms,gamma(1),phz0(1),pht0(1),dcpu,ih,im,is
   endif
   dcpu=cpu_total-cpu_after_start
   ! write(irefwr,'(10x,''cpu = '',4f8.1,5x,4f5.1)') dcpu,cpu_stokes,cpu_heat,cpu_tracers,cpu_stokes/dcpu*100,cpu_heat/dcpu*100,cpu_tracers/dcpu*100, & 
   !   & (cpu_stokes+cpu_heat+cpu_tracers)/dcpu*100
   if (itracoption>0) then
     if (its_CH94) then
        if (print_node) write(irefwr,'(''          T: '',2f12.7,11i10:)') tracer(1)%x,tracer(1)%y, &
              & (imeltnow(1,imelt),imelt=1,nmeltzone), &
              & (imeltnow(2,imelt),imelt=1,nmeltzone)
     else if (gable_plates) then
        if (print_node) write(irefwr,'(''          T: '',2f12.7,5i10)') tracer(1)%x,tracer(1)%y, &
              & sum(imeltnow(1,1:nplates)),&
              & sum(imeltnow(2,1:nplates)),needed_neighbors,needed_neighbors_neighbors,missed
     else if (itracoption==2) then
        if (print_node) write(irefwr,'(''          T: '',2f12.7,5i10,:)') tracer(1)%x,tracer(1)%y,imark(1)
     else
        if (print_node) write(irefwr,'(''          T: '',2f12.7,3i10)') tracer(1)%x,tracer(1)%y, &
              & needed_neighbors,needed_neighbors_neighbors,missed
     endif
   endif
   !if (print_node) write(LU_MEM,*) time_now,nextfreesol-1
   if (gable_plates) then
      t_d=time_now/tscale_dim
      if (print_node) write(LU_TSTEP,'(2e15.7)') time_now,t_d
      if (print_node) call flush(LU_TSTEP)
   endif
   if (print_node) then
      write(LU_NU,*) time_now,gnus1
      write(LU_NU_BOT,*) time_now,gnus2
      write(LU_VRMS,*) time_now,vrms
      write(LU_SURFACEVEL,*) time_now,veloc(1),veloc(2)
      write(LU_PHI,*) time_now,avphi(1)*DiRa
      write(LU_WORK,*) time_now,avwork(1)*DiRa
      write(LU_MEM,'(e15.7,f15.7,i5)') time_now,get_resmem(),nextfreesol-1
   endif
   if (print_node) then
      call flush(irefwr)
      call flush(LU_VRMS)
      call flush(LU_NU)
      call flush(LU_NU_BOT)
      call flush(LU_SURFACEVEL)
      call flush(LU_PHI)
      call flush(LU_WORK)
      call flush(LU_MEM)
   endif
endif

end subroutine time_int_output
