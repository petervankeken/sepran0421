!   TRACVEL
!
!   Determine velocity in markers. Adapted from t4detvel
!
!   PvK 073004 Modified to work with markerchain method (itracoption=2)
!
!   PvK 950508/990407
!   Adapted for use with cylindrical coordinates
!
!   PvK/AKM 051501
!   Corrected for curved elements
!   All elements are treated as straight element; if an element
!   is curved a local transformation to a straight element is
!   performed.
subroutine tracvel(kmesh1,kprob1,isol1,user_here)
use sepmodulekmesh
use control
use coeff
use msper01
use mpetrac
use tracers
use mparallel
use sepran_interface
implicit none
integer :: kmesh1,kprob1,isol1
real(kind=8) :: user_here(*)

integer :: ntot,itrac,j,k,iel
real(kind=8) :: xm,ym,xn(6),yn(6),shapef(6),un(6),vn(6),u,v
integer :: nodno(6),nmark
logical :: edge
integer :: nodlin(3),icorrect,isub,imissed
real(kind=8) :: r,th,rl(3)
logical :: midflag,dotflag
real(kind=8) :: xi,eta
integer :: ndist2,ntrac2
real(kind=4) :: t0


call cpu_time(t0)
! make velocity available in user()
call pecopy(2,user_here(10:10+npoint-1),isol1)
call pecopy(3,user_here(10+npoint:10+2*npoint-1),isol1)
call cpu_time(t1)
if (petest.and.print_node) write(irefwr,*) 'pecopy: ',t1-t0

! find total number of markers
if (itracoption == 2) then
   nmark=0
   do ichain=1,nochain
      nmark = nmark + imark(ichain)
   enddo
   ntot=nmark
else !! includes itracoption=3,9
   ntot=0
   do idist=1,ndist
      ntot = ntot + ntrac(idist)
   enddo
endif
!write(lu_out,*) 'ntot = ',myid,ntot,ntrac(1:ndist)

if (.not.xi_eta_stored) then
!  tracer element info is out of date
   call detelemtrac(1,'tracvelc')
endif
call cpu_time(t1)
 
do itrac=1,ntot
   iel = tracer(itrac)%ielem
   xi  = tracer(itrac)%xi !xi_eta(1,itrac)
   eta = tracer(itrac)%eta !xi_eta(2,itrac)
   !write(lu_out,*) itrac,iel,xi,eta
   call sper01(kmeshc,coor,nodno,xn,yn,iel)
   call sper06(user_here(10:10+npoint-1),un,nodno)
   call sper06(user_here(10+npoint:10+2*npoint-1),vn,nodno)
   call getshape6_xi_eta(xi,eta,shapef)
   u=0
   v=0
   do k=1,6
      u = u + shapef(k)*un(k)
      v = v + shapef(k)*vn(k)
   enddo
   ! write(irefwr,*) xi,eta,iel,u,v
   tracer(itrac)%u=u
   tracer(itrac)%v=v
enddo

end subroutine tracvel
