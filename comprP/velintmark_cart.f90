subroutine velintmark_cart(xm,ym,u,v,ielh,xi,eta,user,callingroutine)
use sepmodulekmesh
use sepmodulecomio
use control
use msper01
implicit none
real(kind=8) :: xm,ym,u,v,user(:),xi,eta
integer :: ielh,k,ludcmp_error
character(len=*) :: callingroutine
real(kind=8) :: un(6),vn(6),xn(6),yn(6),shapef(6)
integer :: nodno(6),icheckinelem
logical :: fail

!write(6,*) 'pedetel'
call pedetel(1,xm,ym,ielh)
!write(6,*) 'sper01'
call sper01(kmeshc,coor,nodno,xn,yn,ielh)
if (icheckinelem(xn,yn,xm,ym).ne.1.and.print_node) then
   write(irefwr,*) 'PERROR(velintmark): sper01 missed this one'
   write(irefwr,*) 'xm, ym, ielh: ',xm,ym,ielh
   write(irefwr,*) 'calling routine: ',callingroutine
   write(irefwr,*) 'npoint: ',npoint
   write(irefwr,*) 'nodno: ',nodno
   write(irefwr,*) 'xn ',xn
   write(irefwr,*) 'yn ',yn
   call instop
endif

!write(6,*) 'sper06'
call sper06(user(10:10+npoint-1),un,nodno)
call sper06(user(10+npoint:10+2*npoint-1),vn,nodno)
ludcmp_error=0
call find_xi_eta(xn,yn,xm,ym,xi,eta,shapef,fail,ludcmp_error)
if (ludcmp_error/=0.and.print_node) then
   write(irefwr,*) 'PERROR(velintmark_cart): find_xi_eta failed on error: ',ludcmp_error
   call instop
endif
u=0
v=0
do k=1,6
   u=u+shapef(k)*un(k)
   v=v+shapef(k)*vn(k)
enddo
!write(6,*) 'u,v: ',u,v

end subroutine velintmark_cart
