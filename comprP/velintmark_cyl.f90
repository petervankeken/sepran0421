! ioffset: extra offset in array user for velocity information 
subroutine velintmark_cyl(ioffset,xm,ym,u,v,ielh,xi,eta,callingroutine,itrac)
implicit none
integer :: ioffset
real(kind=8) :: xm,ym,u,v
integer :: ielh,k,itrac
character(len=*) callingroutine
real(kind=8) :: un(6),vn(6),xn(6),yn(6),shapef(6)
integer :: nodno(6),icorrect,imissed,ifound
integer :: nodlin(3),isub
real(kind=8) :: rl(3)
logical :: midflag,dotflag
logical :: guess_first
real(kind=8) :: xi,eta

!     *** find velocity and coordinates for the element in which
!     *** this tracer is located
call pedeteltrac(1,ioffset,xm,ym,nodno,nodlin,rl,xn,yn,un,vn,ielh,imissed,ifound, &
     &     shapef,xi,eta,guess_first,callingroutine,itrac)

u=0
v=0
do k=1,6
   u=u+shapef(k)*un(k)
   v=v+shapef(k)*vn(k)
enddo

end subroutine velintmark_cyl
