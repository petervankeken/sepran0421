subroutine steady_iteration()
use sepmodulecomio
use sepmoduleoldrouts
use control
use sepran_arrays
use coeff
implicit none
integer :: iinvec(20)=0.0_8,nitertosi=0,isoltosi=0,ihelp
real(kind=8) :: rinvec(10)=0.0_8,diftosi,veltosimax

niter=0
do
  niter=niter+1
  call bmstokes()
  nitertosi=0
  if (tosi15.and.(tosi15_case==2.or.tosi15_case==4)) then
     do
        nitertosi=nitertosi+1
        call copyvc(isol(1),isoltosi)
        call bmstokes()
        diftosi=anorm(0,3,0,kmesh,kprob,isol(1),isoltosi,ihelp)
        veltosimax=anorm(1,3,0,kmesh,kprob,isol(1),isoltosi,ihelp)
        diftosi=diftosi/veltosimax
        write(6,'(''      tosi: '',i5,f12.7)') nitertosi,diftosi
!       write(6,*) 'stopping after first Tosi iteration'
!       call instop
        if (diftosi<1e-4 .or. nitertosi>30) exit
        if (relaxtosi>1e-3) then
             iinvec=0
             iinvec(1)=11
             iinvec(2)=27 ! linear combination
             rinvec(1)=1.0_8-relaxtosi
             rinvec(2)=relaxtosi
             call manvec(iinvec,rinvec,isoltosi,isol(1),isol(1),kmesh,kprob)
        endif
     enddo
  endif


  ! solve steady state heat equation
  call steady_heatP()

  if (relax>1e-7_8) then
      iinvec=0
      iinvec(1)=11
      iinvec(2)=27 ! linear combination
      iinvec(11)=2 ! iprob
      rinvec(1)=1.0_8-relax
      rinvec(2)=relax
      call manvec(iinvec,rinvec,islold(2),isol(2),isol(2),kmesh,kprob)
  endif
  call intermediate_output()

  call copyvc(isol(1),islold(1))
  call copyvc(isol(2),islold(2))
  
  if (niter>=nsteady_max.or.dif<eps_convergence) exit
end do
call bmout()
end subroutine steady_iteration
