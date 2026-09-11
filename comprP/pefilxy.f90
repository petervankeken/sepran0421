! PEFILXY
! 
! Fill the coordinates of the primary nodes of a rectangular
! mesh in common /pexcyc/
!    
! 
! PvK 2/1/99 Updated to allow for more than one curve along
!            either the bottom or side boundary. Information
!            on the curve definition should be stored in 
!            cfilxy.inc
subroutine pefilxy(pvk_ishape)
use sepmodulecomio
use geometry
implicit none
integer :: pvk_ishape

real(kind=8) :: funcy,dx,dy,dxmin,dymin
integer :: number,i,ic,nc
integer :: icurvs(12)

call oldpefilxy(pvk_ishape)

end subroutine pefilxy

! OLDPEFILXY
! 
! Fill the coordinates of the primary nodes of a rectangular
! mesh in common /pexcyc/
! 
! 
! PvK 10-8-89
subroutine oldpefilxy(pvk_ishape)
use sepmodulecomio
use sepmodulecpack
use control
use geometry
use sepran_arrays
implicit none
integer :: pvk_ishape
integer,parameter :: NMAX=10000
real(kind=8) ::funcx1(NMAX),funcy1(NMAX)

real(kind=8) :: dx,dy,dxmin,dymin
integer :: icurvs(12),number=0,i,ic

!pedebug=.true.
pedebug=.false.

funcx1(1) = NMAX
funcy1(1) = NMAX
icurvs(1)=0
icurvs(2)=1
call compcr(-1,kmesh1,kprob1,isol1,number,icurvs,funcx1,funcy1)
if (pedebug.and.print_node) write(irefwr,*) 'curve1: ',funcx1(5),inodenr
if (nint(funcx1(5)/2).gt.NXCYCMAX) then
   if (print_node) then
      write(irefwr,*) 'PERROR(pefilxy) Mesh is too large in x-direction.'
      write(irefwr,*) 'Limit number of x-nodal points to NXCYCMAX'
      write(irefwr,*) 'nx       = ',funcx1(5)/2
      write(irefwr,*) 'NXCYCMAX = ',NXCYCMAX
   endif
endif
if (pvk_ishape==1) then
!  linear elements
   nx = nint(funcx1(5)/2)
   do i=1,nx
      xc(i) = funcx1(4+2*i)
   enddo
else if (pvk_ishape==2) then
!  quadratic elements
   if (pedebug.and.print_node) write(irefwr,*) 'funcx1/4 = ',funcx1(5)/4
   nx = int(funcx1(5)/4+1)
   do i=1,nx
      xc(i) = funcx1(2+4*i)
   enddo
endif

icurvs(1)=0
icurvs(2)=2
call compcr(-1,kmesh1,kprob1,isol1,number,icurvs,funcx1,funcy1)
if (funcx1(5)/2.gt.NXCYCMAX) then
   if (print_node) then
      write(irefwr,*) 'PERROR(pefilxy) Mesh is too large in y-direction.'
      write(irefwr,*) 'Limit number of y-nodal points to NXCYCMAX'
      write(irefwr,*) 'ny       = ',funcx1(5)/2
      write(irefwr,*) 'NXCYCMAX = ',NXCYCMAX
   endif
   call instop
endif
if (pedebug.and.print_node) write(irefwr,*) 'pvk_ishape= ',pvk_ishape,funcx1(5),inodenr
if (pvk_ishape==1) then
!  linear elements
   ny = nint(funcx1(5)/2)
   do i=1,ny
      yc(i) = funcx1(5+2*i)
   enddo
else if (pvk_ishape==2) then
!  quadratic elements
   ny = int(funcx1(5)/4+1)
   do  i=1,ny
      yc(i) = funcx1(3+4*i)
   enddo
endif
if (pedebug) write(irefwr,*) 'nx = ',nx
if (pedebug) write(irefwr,*) 'xc: ',xc(1:nx)
if (pedebug) write(irefwr,*) 'ny = ',ny
if (pedebug) write(irefwr,*) 'yc: ',yc(1:ny)
      
xcmin = xc(1)
ycmin = yc(1)
xcmax = xc(nx)
ycmax = yc(ny)
if (pedebug) write(irefwr,*) 'xcmax, ycmax=',xcmax,ycmax

dxmin = xc(2)-xc(1)
do ic=2,nx-1
   dx = xc(ic+1)-xc(ic)
   dxmin = min(dxmin,dx)
enddo
 
dymin = yc(2)-yc(1)
do ic=2,ny-1
   dy = yc(ic+1)-yc(ic)
   dymin = min(dymin,dy)
enddo
      
if (pvk_ishape==2) then 
   dxmin = dxmin/2 
   dymin = dymin/2
endif
! assume no grid refinement in horizontal direction
dx_two_elements=4*(xcmax-xcmin)/nx
!write(irefwr,*) 'dx_two_elements = ',dx_two_elements,xcmax,xcmin,nx

if (pedebug) then
   pedebug=.false.
!  call instop
endif

end subroutine oldpefilxy

