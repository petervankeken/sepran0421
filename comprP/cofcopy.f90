!   cofcopy
!
!   Copy first dof of isol2 into first dof of isol1
!
!   Note: as of 1995 or so, solution vectors are not renumbered
!   anymore...
!
!   PvK 012210: expanded for quadratic Taylor-Hood triangle
subroutine cofcopy(ichoice,isol1,isol2,kmesh1,kprob1,factor)
use sepmodulekmesh
use sepmodulekprob
use sepmodulesol
implicit none
integer :: ichoice,isol1,isol2,kmesh1(*),kprob1(*)
real(kind=8) :: factor
integer :: nunkp1,nunkp2,nusol1,nusol2,nphys1,nphys2
 
!call sepgetmeshinfo(ndim,npoint)
call sepgetprobinfo(nphys1,nunkp1,nusol1,isol1)
call sepgetprobinfo(nphys2,nunkp2,nusol2,isol2)
call sepactsolbf1(isol1)
 
call cofcop01(ks(isol1)%sol,ks(isol2)%sol,indprfi,indprpi,kprobfi,kprobpi,npoint,nphys1,nunkp1,nunkp2,factor) 

return
end subroutine cofcopy

!   COFCOP01
!   Copy temperature from sol2 to sol1. Sol1 is used in 
!   the specification of the coefficients.
!   PvK 9906
subroutine cofcop01(sol1,sol2,indprf1,indprp1,kprobf1,kprobp1,npoint,nphys,nunkp1,nunkp2,factor)
use control
use sepmodulecomio
implicit none
real(kind=8) :: sol1(*),sol2(*)
integer :: npoint,nunkp1,nunkp2,indprf1,kprobf1(*)
integer :: kprobp1(npoint,nphys),nphys,indprp1
real(kind=8) :: factor
integer :: i,j,j1,j2

if (indprf1 == 0 .and. indprp1== 0) then
  do i=1,npoint
     j1 = (i-1)*nunkp1 + 1
     sol1(j1) = factor*sol2(i)
  enddo
else if (indprf1 /= 0 .and. indprp1 == 0) then
! variable number of dofs
  do i=1,npoint
     j1 = kprobf1(i) + 1
     if (print_node) write(irefwr,'('' cofcop01: '',3i7)') i,j1,1+(i-1)*3
     sol1(j1) = factor*sol2(i)
  enddo
else
! write(irefwr,*) 'use kprob p'
  do i=1,npoint
!    write(irefwr,'(''KPROB P: '',3i10)') (kprobp1(i,j),j=1,nphys)
     j1 = kprobp1(i,1)
     sol1(j1) = factor*sol2(i)
  enddo
     
endif
!do i=1,npoint
!   write(irefwr,*) 'sol: ',sol1(i),sol2(i)
!enddo

return
end subroutine cofcop01

subroutine f12copy(islol1,isol2,isol1,ipress,kmesh1,kprob1)
use sepmodulecomio
use sepmodulekmesh
use sepmodulekprob
use sepmodulesol
implicit none
integer :: islol1,isol2,isol1,ipress,kmesh1(*),kprob1(*)
integer :: nunkp1,nphys1,nusol1

call sepgetprobinfo(nphys1,nunkp1,nusol1,islol1)
call sepactsolbf1(islol1)

call f12cop01(ks(islol1)%sol,ks(isol2)%sol,ks(isol1)%sol,ks(ipress)%sol,indprfi,kprobfi, & 
     & indprpi,kprobpi,nphys1,coor,npoint,nunkp1)

end subroutine f12copy

!   F12COP01
!   Copy f1+f2 into 1st and 2nd dof of sol1. Sol1 is used in 
!   the specification of the coefficients.
!   sol2 contains the temperature.
!   PvK 050200
subroutine f12cop01(slol1,sol2,sol1,press,indprf1,kprobf1,indprp1,kprobp1,nphys,coor,npoint,nunkp1)
use sepmodulecomio
use control
use coeff
implicit none
real(kind=8) :: slol1(*),sol1(*),sol2(*),coor(2,*),press(*)
integer :: npoint,indprf1,kprobf1(*),nunkp1,nphys
integer :: kprobp1(npoint,nphys),indprp1
real(kind=8) :: x,y,pefbuoy1,pressure,f1,f2,fmag
integer :: i,j,j1,j2,j3
real(kind=8) :: xmin,xmax,ymin,ymax,f1min,f1max,f2min,f2max
real(kind=8) :: fmmax,fmmin,loadmax1,loadmax2
data fmmin,xmin,ymin,f1min,f2min/5*1e9_8/
data fmmax,xmax,ymax,f1max,f2max/5*0e0_8/
data loadmax1,loadmax2/2*0d0/


if (indprf1 == 0 .and. indprp1 == 0) then
  do i=1,npoint
     x=coor(1,i)
     y=coor(2,i)
     j1 = (i-1)*nunkp1 + 1
     pressure=0
     slol1(j1) = pefbuoy1(1,x,y,sol2(i),pressure)
     j1 = j1+1
     slol1(j1) = pefbuoy1(2,x,y,sol2(i),pressure)
     f1=slol1(j1-1)
     f2=slol1(j1)
     if (i/100*100==i.and.print_node) write(irefwr,'(''F1/2a: '',5e15.7)') x,y,slol1(j1-1),slol1(j1),pressure
  enddo
else if (indprp1 /= 0 ) then
! *** Use kprob part P for renumbered vectors
  do i=1,npoint
     x=coor(1,i)
     y=coor(2,i)
     j1 = kprobp1(i,1)
     j2 = kprobp1(i,2)
     j3 = kprobp1(i,3)
     pressure = 0.0_8
     slol1(j1) = pefbuoy1(1,x,y,sol2(i),pressure)
     slol1(j2) = pefbuoy1(2,x,y,sol2(i),pressure)
!    f1=slol1(j1)
!    f2=slol1(j2)
     if (i/100*100==i.and.print_node) write(irefwr,'(''F1/2b: '',i5,6e15.7)') i,x,y,slol1(j1),slol1(j2),pressure,sol2(i)
  enddo
else 
! indrprf /= 0
  do i=1,npoint
     x=coor(1,i)
     y=coor(2,i)
     j1 = kprobf1(i) + 1 
     pressure = 0.0_8
     slol1(j1) = pefbuoy1(1,x,y,sol2(i),pressure)
     j1 = j1+1
     slol1(j1) = pefbuoy1(2,x,y,sol2(i),pressure)
     if (i/100*100==i.and.print_node) write(irefwr,'(''F1/2c: '',5e15.7)') x,y,slol1(j1-1),slol1(j1),pressure
!    f1=slol1(j1-1)
!    f2=slol1(j1)
  enddo
endif

!do i=1,npoint
!   write(irefwr,'(''f12: '',4e15.7)') coor(1,i),coor(2,i),sol2(i),slol1(i)
!enddo

end subroutine f12cop01
