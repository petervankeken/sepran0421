module sepran_interface
implicit none
interface
   subroutine sol2vtu(ichoice,isol,basename,idegfd,namedof)
      integer,intent(in) :: ichoice ! 1=single scalar; 2=(u,v,p) and output (u,v) as vector, p as scalar
      integer,intent(in) :: isol  ! standard sepran array
      character(len=*),intent(in) :: basename,namedof(:)
   end subroutine sol2vtu
 
   subroutine pefilbuoy(ic,fbuoy,isol,adia,ipress,idens,ialpha)
     integer,intent(in) :: ic,isol,adia,ipress,idens,ialpha
     real(kind=8),intent(inout) :: fbuoy(:)
   end subroutine pefilbuoy

   subroutine sepgetprobinfo(nphys_here,nunkp_here,nusol_here,isol_here)
     integer :: nphys_here, nunkp_here, isol_here, nusol_here
   end subroutine sepgetprobinfo

   subroutine sepgetmeshinfo(ndim_here,npoint_here)
     integer :: ndim_here,npoint_here
   end subroutine sepgetmeshinfo

   subroutine sepgetrhsdinfo(irhsd)
     integer :: irhsd
   end subroutine sepgetrhsdinfo

   subroutine pvk_allocate(array,length_array,name_array)
     real(kind=8),allocatable :: array(:)
     integer,intent(in) :: length_array
     character(len=*),intent(in) :: name_array
   end subroutine pvk_allocate

  subroutine pecopy(ichoice,user,isol)
     integer,intent(in) :: ichoice,isol
     real(kind=8),intent(inout) :: user(*)
  end subroutine pecopy

  subroutine surfacevel(kmesh,kprob,isol,veloc,veloc_d,noutje)
     integer :: kmesh,kprob,isol,noutje
     real(kind=8) :: veloc(:),veloc_d(:)
  end subroutine surfacevel

  subroutine averageT(ichoice,iuser,user,avT)
     integer :: ichoice,iuser(:)
     real(kind=8) :: user(:),avT(:)
  end subroutine averageT

  subroutine phifromgrad(kmesh1,kprob1,isol1,ivisdip)
     integer :: kmesh1,kprob1,isol1,ivisdip
  end subroutine phifromgrad

  logical function checkinelem(ichoice,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,ielh)
       integer :: ichoice,nodno(6),nodlin(3),isub,ielh
       real(kind=8) :: xn(6),yn(6),xm,ym,rl(3),xi,eta,phiq(6)
  end function checkinelem

  subroutine velintmark_cart(xm,ym,u,v,ielh,xi,eta,user,calling_routine)
       integer :: iehl
       real(kind=8) :: xm,ym,u,v,xi,eta
       character(len=*) :: calling_routine
       real(kind=8) :: user(:)
  end subroutine velintmark_cart
! subroutine move_tracers4(ichoice,kmesh,kprob,isolold,isolnew,user,tstep,tfach)
!      integer :: ichoice,isolold,isolnew,kmesh,kprob
!      real(kind=8) :: user(:),tstep,tfach
! end subroutine move_tracers4

  integer function icheckinelem(xn,yn,xm,ym)
        real(kind=8) :: xn(:),yn(:),xm,ym
  end function icheckinelem

end interface
end module sepran_interface
