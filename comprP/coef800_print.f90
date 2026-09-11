
subroutine coef800_print(npoint,coor,u,v,beta,right,rhocp,temp,cond)
use sepmodulecomio
implicit none
integer :: npoint,i
real(kind=8),dimension(*) :: u,v,beta,right,rhocp,temp,cond
real(kind=8),dimension(2,*) :: coor

write(irefwr,'(''coef800: '',9a12 )') 'x','y','u','v','beta','right','rhocp','temp','cond'
do i=1,npoint,npoint/10
   write(irefwr,'(''         '',9f12.3)') coor(1,i),coor(2,i),u(i),v(i),beta(i),right(i),rhocp(i),temp(i),cond(i)
enddo

end subroutine coef800_print
