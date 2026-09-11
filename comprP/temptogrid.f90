!   TEMPTOGRID
!
!   Determine temperature in the regular grid nodal points
!   from kmesh1
!
!   See velintmarkc() for a description of the interpolation logic
!   PvK 990908
!
!   Updated to allow for grids with variable angular extent
!   (added check_cyl_angular_extent).
!   PvK 091605
subroutine temptogrid_with_trac(kmesh,kprob,isolh,user_here,temp)
use sepmodulecomio
use coeff
use geometry
use control
implicit none
integer :: nth,nr
integer :: isolh,kmesh,kprob
real(kind=8) :: user_here(*),temp(*)
real(kind=8) :: thmin,thmax
integer :: ith,ir,ip
integer :: nodno(6),nodlin(6),imissed,isub,nelem,nelgrp,iel,ichoice,ioffset,ifound,k,ielh,itrac
real(kind=8) :: xi,eta,r,th,temperature,dth,dr,un(6),vn(6),shapef(6),xn(6),yn(6),rl(6),xm,ym
logical :: first,guess_first
save first
data first / .true. /

call pecopy(0,user_here,isolh)

if (first) then
   call check_cyl_angular_extent(kmesh,kprob,isolh,thmin,thmax)
   first=.false.
endif
nth=nth_output
nr=nr_output
dth=dth_output
dr=dr_output

dth=(thmax-thmin)/nth
dr=(r2-r1)/nr
guess_first=.false.

do ir=1,nr
   r=r1+(ir-0.5)*dr
   do ith=1,nth
      th=thmin+(ith-0.5)*dth
      xm=r*sin(th)
      ym=r*cos(th)
      ichoice=2
      ioffset=0
      call pedeteltrac(ichoice,ioffset,xm,ym,nodno,nodlin,rl,xn,yn,un,vn,ielh,imissed, &
     &          ifound,shapef,xi,eta,guess_first,'detelemtrac',itrac)
      if (ierror/=0) then
        if (print_node) write(irefwr,*) 'ierror=',ierror,' in temptogrid'
        call instop
      endif
      temperature=0
      do k=1,6
         temperature=temperature+shapef(k)*un(k)
      enddo
      ip=(ith-1)*nth+ir
      temp(ip)=temperature
   enddo
enddo


end subroutine temptogrid_with_trac


subroutine temptogrid(kmesh,kprob,isolh,user_here,temp)
use sepmodulecomio
use coeff
use geometry
use control
implicit none
integer :: nth,nr
integer :: isolh,kmesh,kprob
real(kind=8) :: user_here(*),temp(*)
integer,parameter :: NMAX=401*401
real(kind=8) :: thmin,thmax,coor(2,NMAX)
integer :: iinmap(2),map(5),NDIM,NTOT,NUNKP,ir,ith,ip
real(kind=8) :: dth,dr,xm,ym,th,r
logical :: first=.true.
save iinmap,map,first,coor

call sepactsolbf1(isolh)

iinmap(1)=2
iinmap(2)=2
if (first) iinmap(2)=1
nth=nth_output
nr=nr_output
if (print_node) write(irefwr,*) 'nth,nr: ',nth,nr
dth=(thmax-thmin)/nth
dr=(r2-r1)/nr
if (nth*nr>NMAX) then
   if (print_node) write(irefwr,*) 'PERROR(temptogrid): nth*nr>NMAX: ',nth,nr,nth*nr,NMAX
   call instop
endif

if (first) then
  ip=0
  do ir=1,nr
     r=r1+(ir-0.5)*dr
     do ith=1,nth
        ip=ip+1
        th=thmin+(ith-0.5)*dth
        xm=r*sin(th)
        ym=r*cos(th)
        coor(1,ip)=xm
        coor(2,ip)=ym
     enddo
  enddo
  first=.false.
endif
NTOT=nth*nr
NDIM=2
NUNKP=1
if (print_node) write(irefwr,*) 'ntot=',NTOT,nr,nth
call intcoor(kmesh,kprob,isolh,temp,coor,NUNKP,NTOT,NDIM,iinmap,map)

end subroutine temptogrid

subroutine check_cyl_angular_extent(kmesh,kprob,isol,thmin,thmax)
use sepmodulecomio
use geometry
use control
implicit none

integer  :: kmesh,kprob,isol
real(kind=8) :: thmin,thmax
integer, parameter :: NMAX=9000
real(kind=8) :: funcx(NMAX),funcy(NMAX)
real(kind=8) :: x0,y0,x1,y1,th0,th1,r01,r11
integer :: icurvs(2),npoint_on_curve

funcx(1)=NMAX
funcy(1)=NMAX

icurvs(1)=0
icurvs(2)=itop

icurvs(1)=0
icurvs(2)=itop
call compcr(-1,kmesh,kprob,isol,0,icurvs,funcx,funcy)
npoint_on_curve = funcx(5)/2
x0 = funcx(6)
y0 = funcx(7)
x1 = funcx(4+2*npoint_on_curve)
y1 = funcx(5+2*npoint_on_curve)
r01 = sqrt(x0*x0+y0*y0)
r11 = sqrt(x1*x1+y1*y1)
th0 = acos(y0/r01)
th1 = acos(y1/r11)
thmin = min(th0,th1)
thmax = max(th0,th1)

icurvs(1)=0
icurvs(2)=ibottom
call compcr(-1,kmesh,kprob,isol,0,icurvs,funcx,funcy)
npoint_on_curve = funcx(5)/2
x0 = funcx(6)
y0 = funcx(7)
x1 = funcx(4+2*npoint_on_curve)
y1 = funcx(5+2*npoint_on_curve)
r01 = sqrt(x0*x0+y0*y0)
r11 = sqrt(x1*x1+y1*y1)
th0 = acos(y0/r01)
th1 = acos(y1/r11)
thmin = min(th0,th1,thmin)
thmax = max(th0,th1,thmax)

end subroutine check_cyl_angular_extent


