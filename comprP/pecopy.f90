subroutine pecopy(ichoice,user,isol)
use sepmodulecomio
use sepmodulemain
use sepmodulekmesh
use sepmodulekprob
use sepmodulesol
use control
implicit none
interface
  subroutine pecopy01(ichoice,ndim,npoint,user,nusol,usol,indprpi,kprobpi)
     integer,intent(in) :: ichoice,ndim,npoint,nusol,indprpi,kprobpi(:,:)
     real(kind=8),intent(in) :: usol(nusol)
     real(kind=8),intent(inout) :: user(*)
  end subroutine pecopy01
end interface     
integer,intent(in) :: ichoice,isol
real(kind=8),intent(inout) :: user(*)

! Make sure everything points to isol
call sepactsolbf1(isol)

if (ichoice < 0 .or. ichoice==1 .or. ichoice > 3) then
   if (print_node) write(irefwr,*) 'PERROR(pecopy) not suited for ichoice other than 0,2,3 but: ',ichoice
   call instop
endif
 
if (npelm/=6) then
   if (print_node) write(irefwr,*) 'PERROR(pecopy) not suited for npelm other than 6'
   call instop
endif
if (ndim/=2) then
   if (print_node) write(irefwr,*) 'PERROR(pecopy) not suited for ndim other than 2'
   call instop
endif
if (debug) write(irefwr,*) 'pecopy: indprpi=',indprpi
if (indprpi==0 .and. ichoice>0) then
   if (print_node) write(irefwr,*) 'PERROR(pecopy) expect kprob P here but indprpi = ',indprpi
   if (print_node) write(irefwr,*) 'ichoice = ',ichoice
   call instop
endif
call pecopy01(ichoice,ndim,npoint,user,nusol,ks(isol)%sol,indprpi,kprobpi)
end subroutine pecopy

subroutine pecopy01(ichoice,ndim,npoint,user,nusol,usol,indprp,kprobp)
use control
use sepmodulecomio
implicit none
integer,intent(in) :: ichoice,ndim,npoint,nusol,indprp,kprobp(:,:)
real(kind=8),intent(in) :: usol(nusol)
real(kind=8),intent(inout) :: user(*)
integer :: i,ip

if (ichoice==0) then
   !do i=1,npoint
   !   user(i)=usol(i)
   !enddo
   user(1:npoint) = usol(1:npoint)
   !if (pedebug) write(irefwr,*) 'min/max user: ',minval(user(1:npoint)),maxval(user(1:npoint))
else
  do i=1,npoint
     ip=kprobp(i,ichoice-1)
     user(i)=usol(ip)
  enddo
endif

end subroutine pecopy01
