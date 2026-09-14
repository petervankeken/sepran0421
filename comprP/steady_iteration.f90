subroutine steady_iteration
use sepmoduleoldrouts
use sepmodulecomio  ! provides lu's for reading & writing
use sepmodulecpack  ! parallel computing interface
use sepran_arrays ! access to kmesh, kprob, isol, etc. and interfaces to subroutines
use sepran_interface
use brandenburg
use control
use dtm_elem  ! controls building of pressure mass matrix
implicit none
integer :: ndim,nphys,npoint,nunkp
integer :: modelv,mconv,nparm,nsup1,idum,iread,ihelp,icurvs(2)
integer :: idegfd,i,ichoice,isubr,ncntln
real(kind=8) :: rayleigh,umax,outputtime,contln(10),tempmax
character(len=120) :: namedof(4)
integer :: ichvc,iprob_here,ip,ipbuoy2,ipuser,ipee,ih,im,is,iinvec(11)
real(kind=8) :: q_a(10),gnus1,gnus2,vrms,rinvec(2)
real(kind=4) :: t0,timeA,timeB

ichoice=0
gable_stokes_choice=0
niter=0
!call cpu_time(timeA)
!call stokesP(1)
!call cpu_time(timeB)
!write(irefwr*) 'initial stokes solution took: ',timeB-timeA
!if (gable_plates) then
!   call cpu_time(timeA)
!   call do_gable_solve(iuser_here,user_here,niter,1)
!   call cpu_time(timeB)
!   write(irefwr*) 'gable stokes solution took: ',timeB-timeA
!endif

call cpu_time(t00)
do
   call cpu_time(t0)
   niter=niter+1
   call copyvc(isol1,isolold1(1))
   call copyvc(isol2,isolold2)
   call steady_heatP()
!  call bmout()
!  if (print_node) write(irefwr,*) 'stop after bmheat/bmout'
!  call instop

   if (relax > 1e-6_8) then
      iinvec=0
      iinvec(1)=11
      iinvec(2)=27 ! linear combination
      iinvec(11)=2 ! iprob
      rinvec(1)=1.0_8-relax
      rinvec(2)=relax
      call manvec(iinvec,rinvec,isolold2,isol2,isol2,kmesh1,kprob1)
      !call algebr(3,0,isolold2,isol2,isol2,kmesh1,kprob1,1d0-relax,relax,pee,quu,ipee)
   endif
   gable_stokes_choice=0
   call cpu_time(timeA)
   if (old_gable_plates) then
      call stokesP(1)
      call pevrms(vrms,isol1)
      !write(irefwr*) 'initial vrms: ',vrms
   else
      call stokesP_new(1)
      !call bmout()
      !write(irefwr*) 'stop after stokesP_new/bmout'
      !call instop
   endif
   !call getvisdip
   !call prinrv(iphi,kmesh1,kprob1,1,1,'work')
   !call instop
   call cpu_time(timeB)
   if (pedebug.and.print_node) then
       call pevrms(vrms,isol1)
       write(irefwr,*) 'initial stokes solution took: ',timeB-timeA,vrms
   endif
   ! ourplotname='plate.000'
   !call plotvc(1,2,isol1,isol1,kmesh1,kprob1,15.0_8,1.0_8,0.0_8)


   if (gable_plates.and.old_gable_plates) then
      !call cpu_time(timeA)
      call do_gable_solve(iuser_here,user_here,niter,1)
      call cpu_time(timeB)
      !do i=1,nplates
      !   write(ourplotname,'(''plate.'',i3.3)') i
      !   call plotvc(1,2,plate_sol(i),plate_sol(i),kmesh1,kprob1,15.0_8,1.0_8,0.0_8)
      !enddo
   endif
      !debug=.true.
      if (debug.and.print_node) write(irefwr,*) 'total gable stokes solution took: ',timeB-timeA
      debug=.false.
   !call instop
   call nusseltP(q_a)
   gnus1=q_a(1)/(q_a(3)-q_a(4))
   dif2=anorm(0,3,0,kmesh1,kprob1,isol2,isolold2,ihelp)
   tempmax=anorm(1,3,0,kmesh1,kprob1,isol2,isol2,ihelp)
   dif2=dif2/tempmax
   call pevrms(vrms,isol1)
   call cpu_time(t1)
   ih=(t1-t00)/3600
   im=((t1-t00)-ih*3600)/60
   is=(t1-t00)-ih*3600-im*60
   cpu_total=t1-t00
   if (print_node) write(irefwr,'('' iter, dif, Nu, vrms: '',i5,3f12.5,f12.2,i3,'':'',i2.2,'':'',i2.2)') &
         niter,dif2,gnus1,vrms,t1-t0,ih,im,is
   !write(irefwr*) 'stop after first stokes' 
   !call instop
   if (niter<min_iter_steady) then
      if(print_node) write(irefwr,*) 'steady_iteration: going again: ',niter,min_iter_steady
      cycle
   endif
   if (dif2<eps .or. niter>=nsteady_max) then
       if (dif2<eps) then
          did_not_converge=.false.
       else
          did_not_converge=.true.
       endif
       exit
   endif
enddo
!t2=second()
!write(irefwr*) 'time: ',t2-t1
inout=1
call bmout()
end subroutine steady_iteration
