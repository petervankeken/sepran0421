! compute azimuthally averaged value of ivec (supposed to be vector of type temperature)
subroutine averages(user)
use coeff
use geometry
use sepmodulekmesh
use sepran_arrays
implicit none
integer :: ipadia
real(kind=8) :: user(*)

if (cyl) then
   !write(irefwr,*) 'PERROR(averages) :: not yet adapted for cyl'
!  write(irefwr,*) 'temptogrid : ',isol2
!  call temptogrid(kmesh1,kprob1,isol2,user,user(10+npoint))
!  ipadia=10+npoint+nth_output*nr_output
!  write(irefwr,*) 'temptogrid : ',iadia
!  call temptogrid(kmesh1,kprob1,iadia,user,user(10+ipadia))
!  call averages_cyl(kmesh1,kprob1,isol2,user(10+npoint),user(10+ipadia))
else
   call pefilxy(1)
   call averages_cart(user(10+npoint))
endif

end subroutine averages

subroutine averages_cart(temp)
use sepmodulecomio
use geometry
use coeff
use sepran_arrays
use control
implicit none
integer :: npoint
real(kind=8) :: temp(*)
integer :: i,ix,iy,ip,it1,it2
real(kind=8) :: x1,x2
character(len=80) :: fname
!real(kind=8), allocatable, dimension(:) :: avetemp
integer, parameter :: NRMAXH=3000
real(kind=8),dimension(NRMAXH) :: avetemp

if (ny>NRMAXH) then
   if (print_node) write(irefwr,*) 'PERROR(averages_cart): enlarge NRMAXH. ny=',ny
   call instop
endif
call pecopy(0,temp,isol2)

avetemp=0
do iy=1,ny
   do ix=1,nx-1
      x1=xc(ix)
      x2=xc(ix+1)
      it1=(iy-1)*nx+ix    
      it2=(iy-1)*nx+ix+1
      avetemp(iy)=avetemp(iy)+0.5*(x2-x1)*(temp(it1)+temp(it2))
   enddo
enddo
avetemp=avetemp/(xcmax-xcmin)

if (nout>0.and.print_node) then
   write(fname,'(''GMT/avetemp.'',i4.4)') nout
   open(78,file=fname)
   do i=1,ny
      write(78,*) avetemp(i),yc(i)
   enddo
   close(78)
endif

end subroutine averages_cart

! AVERAGES_CYL
! Compute and output radial averages of T, eta.
! Then computer average T (and average T in upper/lower layer
! if rmid is defined).
! Should be used after temptogrid has been called
! 
! PvK 990908
subroutine averages_cyl(kmesh2,kprob2,isol2,temp,adia,noutje)
use sepmodulecomio
use geometry
use convparam
use coeff
use geometry
use control
implicit none
integer :: kmesh2(*),kprob2(*),isol2(*),noutje
real(kind=8) :: temp(*),adia(*)
real(kind=8) :: r,th,avetemp(500),avevis(500),dth,dr,x,y,pefvis
real(kind=8) :: secsqrdummy,viscosity,r1l,r0l,davTh,avTh,temperature,avThl,avThu
real(kind=8) :: avevistop,avevisbot,avetempbot,avetemptop,temp_mid
real(kind=8) :: temp_a
integer :: ir,ith,ival,ivalfind1,nrmid,ip,nth,nr
character(len=80) gname

nr  = nr_output
nth = nth_output
dth = dth_output
dr  = dr_output
if (nr.gt.500) then
   write(irefwr,*) 'PERROR(averages) nr > 500: ',nr
   call instop
endif

if (noutje.ge.0) then
   write(gname,'(''GMT/avetemp.'',i4.4)') noutje
   open(78,file=gname)
   write(78,*) 2,nr+2
   if (itypv.gt.0) then
     write(gname,'(''GMT/avevis.'',i4.4)') noutje
     open(77,file=gname)
     write(77,*) 2,nr+2
   endif
endif
do ir=1,nr
   avetemp(ir) = 0
   avevis(ir) = 0
   avevisbot = 0
   avevistop = 0
   r = r1 + (ir-0.5)*dr
   do ith=1,nth
     th = (ith-0.5)*dth
      x = r*sin(th)
      y = r*cos(th)
     ival=1
     ip = (ith-1)*nth+ir
     temperature = temp(ip)
     temp_a = adia(ip)
     if (ivl) ival = ivalfind1(x,y)
     viscosity = pefvis(x,y,temperature,temp_a,ival,secsqrdummy)
     avetemp(ir) = avetemp(ir) + temperature*dth
     avevis(ir) = avevis(ir) + viscosity*dth
   enddo
   avetemp(ir) = avetemp(ir)/(nth*dth)
   avetemptop = 0d0
   avetempbot = t_bot
   avevis(ir) = avevis(ir)/(nth*dth)
   avevistop = avevistop/(nth*dth)
   avevisbot = avevisbot/(nth*dth)
enddo

!     *** bottom layer
r = r1
th = 0
x = r*sin(th)
y = r*cos(th)
ival=1
if (ivl) ival = ivalfind1(x,y)
viscosity = pefvis(x,y,t_bot,t_bot_adia,ival,secsqrdummy)
avevisbot = viscosity
!     *** top layer
r = r2
th = 0
x = r*sin(th)
y = r*cos(th)
ival=1
if (ivl) ival = ivalfind1(x,y)
viscosity = pefvis(x,y,0d0,0d0,ival,secsqrdummy)
avevistop = viscosity

if (noutje.ge.0) then
   if (itypv.gt.0) then
      write(77,*) avevisbot,r1
      do ir=1,nr
        write(77,*) avevis(ir),r1+(ir-0.5)*dr
      enddo
      write(77,*) avevistop,r2
      close(77)
   endif
   write(78,*) avetempbot,r1
   do ir=1,nr
     write(78,*) avetemp(ir),r1+(ir-0.5)*dr
   enddo
   write(78,*) avetemptop,r2
   close(78)
endif


if (axi) return

! FIND AVERAGE T
! Use cylindrical version of trapezoid rule.
! T is defined in the middle of interval.
! 
! <<T>> =      int pi <T(r)> r dr
!              ------------------
!                 int pi r dr

! *** first half interval
r0l = r1
r1l = r1+0.5*dr
davTh = (1d0/6d0*r1l**3-1d0/2*r0l*r0l*r1l+1d0/3d0*r0l*r0l*r0l)*t_bot/(0.5*dr)
davTh = davTh + (1d0/6d0*r0l**3-1d0/2*r0l*r1l*r1l+1d0/3d0*r1l*r1l*r1l)*avetemp(1)/(0.5*dr)
avTh = pi*davTh
do ir=1,nr-1
   r0l = r1 + (ir-0.5)*dr
   r1l = r1 + (ir+0.5)*dr
   davTh = (1d0/6d0*r1l**3-1d0/2*r0l*r0l*r1l+1d0/3d0*r0l*r0l*r0l)*avetemp(ir)/dr
   davTh = davTh + (1d0/6d0*r0l**3-1d0/2*r0l*r1l*r1l+1d0/3d0*r1l*r1l*r1l)*avetemp(ir+1)/dr
   avTh = avTh + pi*davTh
enddo
! last half interval
r0l = r2-0.5*dr
r1l = r2
davTh = (1d0/6d0*r1l**3-1d0/2*r0l*r0l*r1l+1d0/3d0*r0l*r0l*r0l)*avetemp(nr)/(0.5*dr)
davTh = davTh + (1d0/6d0*r0l**3-1d0/2*r0l*r1l*r1l+1d0/3d0*r1l*r1l*r1l)*avetemptop/(0.5*dr)
avTh = avTh + pi*davTh
avTh = avTh / (pi*frac*(r2**2-r1**2))

! *************************************************************
! *   Compute average T in upper/lower layer if necessary
! *
! *   R1               Rmid              R2
! *   |                  |                |
! *     |---|---|---|---|-++|+++|+++|+++|
! *     1   2   3     nrmid 
! *************************************************************
if (rmid.gt.0) then
   nrmid = (rmid-(r1+0.5*dr))/dr+1
!  temperature at Rmid
   temp_mid = avetemp(nrmid)+(avetemp(nrmid+1)-avetemp(nrmid))/dr*(rmid-(nrmid+0.5)*dr)
!  Lower layer
!  First half interval
   r0l = r1
   r1l = r1+0.5*dr
   davTh = (1d0/6d0*r1l**3-1d0/2*r0l*r0l*r1l+1d0/3d0*r0l*r0l*r0l)* t_bot/(0.5*dr)
   davTh = davTh + (1d0/6d0*r0l**3-1d0/2*r0l*r1l*r1l+1d0/3d0*r1l*r1l*r1l)*avetemp(1)/(0.5*dr)
   avThl = davTh*pi
   do ir=1,nrmid-1
     r0l = r1 + (ir-0.5)*dr
     r1l = r1 + (ir+0.5)*dr
     davTh =(1d0/6d0*r1l**3-1d0/2*r0l*r0l*r1l+1d0/3d0*r0l*r0l*r0l)*avetemp(ir)/dr
     davTh = davTh + (1d0/6d0*r0l**3-1d0/2*r0l*r1l*r1l+1d0/3d0*r1l*r1l*r1l)*avetemp(ir+1)/dr
     avThl = avThl + pi*davTh
   enddo
   r0l = r1 + (nrmid-0.5)*dr
   r1l = rmid
   if (r1l-r0l.gt.0) then
    davTh =(1d0/6d0*r1l**3-1d0/2*r0l*r0l*r1l+1d0/3d0*r0l*r0l*r0l)* avetemp(nrmid)/(r1l-r0l)
    davTh = davTh + (1d0/6d0*r0l**3-1d0/2*r0l*r1l*r1l+1d0/3d0*r1l*r1l*r1l)*temp_mid/(r1l-r0l)
    avThl = avThl + pi*davTh
   endif

!  upper layer
   r0l = rmid
   r1l = r1 + (nrmid+0.5)*dr
   if (r1l-r0l.gt.0) then
    davTh = (1d0/6d0*r1l**3-1d0/2*r0l*r0l*r1l+1d0/3d0*r0l*r0l*r0l)*temp_mid/(r1l-r0l)
    davTh = davTh + (1d0/6d0*r0l**3-1d0/2*r0l*r1l*r1l+1d0/3d0*r1l*r1l*r1l)*avetemp(nrmid)/(r1l-r0l)
    avThu = pi*davTh
   endif
   do ir=nrmid+1,nr-1
     r0l = r1 + (ir-0.5)*dr
     r1l = r1 + (ir+0.5)*dr
     davTh =(1d0/6d0*r1l**3-1d0/2*r0l*r0l*r1l+1d0/3d0*r0l*r0l*r0l)* avetemp(ir)/dr
     davTh = davTh + (1d0/6d0*r0l**3-1d0/2*r0l*r1l*r1l+1d0/3d0*r1l*r1l*r1l)*avetemp(ir+1)/dr
     avThu = avThu + pi*davTh
   enddo
   r0l = r2-0.5*dr
   r1l = r2
   davTh = (1d0/6d0*r1l**3-1d0/2*r0l*r0l*r1l+1d0/3d0*r0l*r0l*r0l)*avetemp(nr)/(0.5*dr)
   davTh = davTh + (1d0/6d0*r0l**3-1d0/2*r0l*r1l*r1l+1d0/3d0*r1l*r1l*r1l)*avetemptop/(0.5*dr)
   avThu = avThu + pi*davTh
   avThl = avThl / (pi*frac*(rmid**2-r1**2))
   avThu = avThu / (pi*frac*(r2**2-rmid**2))
endif

end subroutine averages_cyl
