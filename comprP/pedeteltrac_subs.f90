!   ELEMENTCURVED
!   Check whether an element is curved or straight
logical function elementcurved(x,y,eps)
implicit none
real(kind=8) :: x(*),y(*),eps
real(kind=8) :: fx,fy,dx,dy

elementcurved = .false.
fx = abs( x(1)-2*x(2)+x(3) )
fy = abs( y(1)-2*y(2)+y(3) )
dx = abs(x(1)-x(3))
dy = abs(y(1)-y(3))
if (max(fx,fy).gt.max(dx,dy)*eps) then
   elementcurved = .true.
   return
endif
fx = abs( x(3)-2*x(4)+x(5) )
fy = abs( y(3)-2*y(4)+y(5) )
dx = abs(x(3)-x(5))
dy = abs(y(3)-y(5))
if (max(fx,fy).gt.max(dx,dy)*eps) then
   elementcurved = .true.
   return
endif
fx = abs( x(5)-2*x(6)+x(1) )
fy = abs( y(5)-2*y(6)+y(1) )
dx = abs(x(5)-x(1))
dy = abs(y(5)-y(1))
if (max(fx,fy).gt.max(dx,dy)*eps) then
   elementcurved = .true.
   return
endif
 
end

!   PRINTKMESHO
!   Utility to print neighboring elements
subroutine printkmesho(kmesho1,kmesho2,nelem,nelgrp)
use sepmodulecomio
use control
implicit none
integer :: kmesho1,kmesho2(*),nelem,nelgrp
integer :: ip,i

if (print_node) write(irefwr,*) 'kmesho1 (3): ',kmesho1,nelgrp
ip = 0
do i=1,nelem
   if (print_node) write(irefwr,*) kmesho2(ip+i),kmesho2(ip+i+1),kmesho2(ip+i+2)
   ip=ip+3
enddo
return
end


!   CHECK_NEIGHBORS
!
!   Utility that is used if the normal element number lookup
!   mechanism (through ielfromgrid) fails. This can happen
!   when the interpolation grid (defined in /pielgrid/) cannot 
!   accurately represent the fine details of the finite element grid.
!   ielfromgrid() tests 4 different candidate elements (stored in 
!   iel_now(*); it is likely that the correct element is a neighbor
!   of one of these elements. This subroutine provides the 
!   mechanisms to test the boundary elements which are found 
!   by using kmesh part O.
!
!   xm,ym         i     tracer coordinates
!   kmeshc        i     nodal point numbers associated with elements
!   coor          i     coordinates of the nodal points
!   nneighbors    i     Standard number of neighbors for this element
!                       (3 for a triangle)
!   kmesho2       i     kmesh part O2. This contains nneighbors*nelem
!                       entries, where each entry provides the 
!                       element number of the neighbor and in some
!                       cases the boundary number as well (see below)
!   nelem         i     number of elements in the mesh
!   nodno         o     nodal point numbers of correct element (if found)
!   xn,yn         o     coordinates of correct elements (if found)
!   iel_now       i     4 candidate elements provided by ielfromgrid
!   ielh          o     number of correct element (if found)
!   rl            o     linear shape functions in xm,ym
!
!   PvK 012804
subroutine check_neighbors(xm,ym,kmeshc,coor,kmesho1,kmesho2,nelem, nodno,nodlin,xn,yn,iel_now,nnow,ielh, &
    & current_neighbors,ncurrent,rl,phiq,xi,eta)
use sepmodulecomio
use geometry
use mtime
use mpetrac
use sepran_interface
use control
implicit none
real(kind=8) :: xm,ym,coor(2,*),xn(*),yn(*)
real(kind=8) :: rl(*),phiq(*),xi,eta
integer :: kmeshc(*),kmesho1(*),kmesho2(*),nelem,nodno(*),nodlin(*)
integer :: iel_now(*),ielh,ip1,ip,current_neighbors(*),ncurrent,nnow
integer :: inpelm,i,j,k,iel_neighbor,isub,ko,nneighbors
logical :: out

nneighbors=kmesho1(1)
ip=0
ip1=0
if (time_now>0.1682150E-04.and.tracer_number==2645188) then
   verbose=.true.
else
   verbose=.false.
endif
!verbose=.true.
inpelm = 6
if (verbose.and.print_node) write(irefwr,'(''in check_neighbors: '',5i10)') (iel_now(i),i=1,4),nneighbors
i=1
ncurrent=0
100   continue
   ! the nnow possible elements indicated in iel_now
   if (verbose.and.print_node) write(irefwr,'(''j/i/nnow = '',3i15,'' iel_now: '',i15)') tracer_number,i,nnow,iel_now(i)
   if (iel_now(i).gt.0) then
!     neighboring elements
      j=1
200     continue
        ip1 = (iel_now(i)-1)*nneighbors+j
        if (verbose.and.print_node) write(irefwr,*) 'ip1=',ip1
        if (verbose.and.print_node) write(irefwr,*) 'kmesho2: ',kmesho2(1)
        ko = kmesho2(ip1)
!       number of jth neighbor of ielnow(i)
!       either ielem or 
!       (nelem+1)*kboundary
!       of kmesh in Sepran Programmer's Guide
250       continue
          if (ko.gt.nelem) then
             if (verbose.and.print_node) write(irefwr,*) 'ko > nelem: ',ko,nelem
             ko = ko - nelem - 1
             goto 250
          endif
          iel_neighbor=ko
          ncurrent=ncurrent+1
          current_neighbors(ncurrent)=iel_neighbor
          if (verbose.and.print_node) write(irefwr,*) 'iel_neighbor = ',i,j,iel_neighbor,iel_now(i) 
          if (iel_neighbor.gt.0) then
             ip = (iel_neighbor-1)*inpelm
             do k=1,inpelm
                nodno(k) = kmeshc(ip+k)
                xn(k) = coor(1,nodno(k))
                yn(k) = coor(2,nodno(k))
             enddo
     if (verbose.and.print_node) then
       write(109,'(''>'')') 
       do k=1,6
          write(109,*) xn(k),yn(k)
       enddo
       write(109,*) xn(1),yn(1)
     endif
             if (curved_elem(iel_neighbor)) then 
                out = .not.checkinelem(3,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,iel_neighbor)
             else 
                out = .not.checkinelem(1,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,iel_neighbor)
             endif
          else
             out = .true.
          endif
          if (out) then
             if (j.lt.3) then
!               continue on inner loop
                j=j+1
                if (verbose.and.print_node) write(irefwr,*) 'continuing on inner loop: ',i,j,iel_now(i),ip,ip1
                goto 200
             else
!               search for this iel_now element has failed
                ielh = -1
             endif
          else
!            Return correct element number
             ielh = iel_neighbor
             return
          endif
   endif
   if (i.lt.nnow) then
!     outer loop
      i=i+1
      goto 100
   endif

end subroutine check_neighbors

!   FIND_XI_ETA
!   Find the barycentric coordinates (xi,eta) of the point
!   (xr,yr) in the triangle defined by xnodr,ynodr; provide 
!   quadratic shapefunctions through phiq
!
!   Parameters
!      xnodr,ynodr      i    nodal points coordinates
!      xr,yr            i    coordinates  of the interpolation point
!      xi,eta           o    barycentric coordinates 
!      phiq             o    quadratic shapefunctions
!
!   PvK 013004
subroutine find_xi_eta(xnodr,ynodr,xr,yr,xi,eta,phiq,fail,ludcmp_error)
implicit none
real(kind=8) :: xnodr(*),ynodr(*),xr,yr,xi,eta,phiq(*)
logical :: fail
real(kind=8) :: tolx
integer :: ntrial,ludcmp_error

ntrial = 30
tolx = 1e-10
! At input: xi,eta should contain initial estimate
! use mid point of the triangle for initial estimate
xi=1d0/3
eta=1d0/3
call xi_eta_newton(ntrial,xi,eta,tolx,xnodr,ynodr,xr,yr,phiq,fail,ludcmp_error)

end subroutine find_xi_eta

!   XI_ETA_NEWTON
!
!   Find barycentric coordinates (xi,eta) and quadratic 
!   shapefunctions (phiq) of the point (xr,yr) in the 
!   curved quadratic triangle defined by (xnodr,ynodr)
!   using a Newton iteration
!
!   PvK 013004
subroutine xi_eta_newton(ntrial,xi,eta,tolx,xnodr,ynodr,xr,yr,phiq,fail,ludcmp_error)
use sepmodulecomio
use control
implicit none
integer :: ntrial 
real(kind=8) :: xi,eta,tolx,xnodr(*),ynodr(*),xr,yr,phiq(*)
logical :: fail
real(kind=8) :: p(2)
integer :: ND,ludcmp_error
parameter(ND=2)
real(kind=8) :: fvec(ND),fjac(ND,ND),d,errf
integer :: indx(ND)=0,k


k=1
100   continue
  call get_f_and_jac(xi,eta,fvec,fjac,ND,xnodr,ynodr,xr,yr,phiq)
  errf = abs(fvec(1))+abs(fvec(2))
  if (verbose.and.print_node) write(irefwr,'(''Newton: '',6e15.7)') errf,tolx,fvec(1),fvec(2),xi,eta
  if (errf.lt.tolx) then
     fail = .false.
     return
  endif

  p(1) = -fvec(1)
  p(2) = -fvec(2)
  !if (verbose.and.print_node) write(irefwr,*) 'before ludcmp: ',p(1),p(2)
  call ludcmp(fjac,2,ND,indx,d,ludcmp_error)
  if (ludcmp_error/=0) then
     if (print_node) write(irefwr,*) xr,yr,xnodr(1:5:2),ynodr(1:5:2)
     return
  endif
  call lubksb(fjac,2,ND,indx,p)
  xi  = xi  + p(1)
  eta = eta + p(2)
  if (k.lt.ntrial) then
    k = k+1 
    goto 100
  else
    fail = .true.
!         write(irefwr,*) 'PERROR(xi_eta_newton): Newton iteration' 
!         write(irefwr,*) 'failed to converge after ',ntrial,' iterations' 
!         write(irefwr,*) 'xi, eta: ',xi,eta
  endif
 
return
end subroutine xi_eta_newton

!   GET_F_AND_JAC
!
!   Compute the mismatch vector and Jacobian matrix 
!   PvK 013004
subroutine get_f_and_jac(xi,eta,fvec,fjac,n,xnodr,ynodr,xr,yr,phiq) 
implicit none
integer :: n
real(kind=8) :: xi,eta,fvec(n),fjac(n,n),xnodr(*),ynodr(*),xr,yr,phiq(*)
real(kind=8) :: phil(3),dphidksi(6),dphideta(6)
integer :: i,j

call getshape_quad_xi_eta(xi,eta,phil,phiq,dphidksi,dphideta)
fvec(1)   = 0
fvec(2)   = 0
fjac(1,1) = 0
fjac(1,2) = 0
fjac(2,1) = 0
fjac(2,2) = 0
do j=1,6
   fvec(1)    =  fvec(1) + phiq(j)*xnodr(j)
   fvec(2)    =  fvec(2) + phiq(j)*ynodr(j)
   fjac(1,1)  =  fjac(1,1) + dphidksi(j)*xnodr(j)
   fjac(1,2)  =  fjac(1,2) + dphideta(j)*xnodr(j)
   fjac(2,1)  =  fjac(2,1) + dphidksi(j)*ynodr(j)
   fjac(2,2)  =  fjac(2,2) + dphideta(j)*ynodr(j)
enddo   
fvec(1) = fvec(1) - xr
fvec(2) = fvec(2) - yr

end subroutine get_f_and_jac

!   getshape_quad_xi_eta
!   Get shape functions and derivatives of quadratic reference 
!   triangle in point (xi,eta)
!
!   Parameters:
!       xi,eta      i    trial point
!       phil        o    linear shape functions 
!       phiq        o    quadratic shape functions
!       dphidksi    o    dphiq/dx in 
!       dphideta    o    dphiq/dy in
!       nl          i    number of linear shapefunctions (3)
!       nq          i    number of quadratic shapefunctions (6)
!
!   PvK 013004
subroutine getshape_quad_xi_eta(xc,yc,phil,phiq,dphidksi,dphideta)
implicit none
real(kind=8) :: xc,yc
real(kind=8) :: phil(3),phiq(6),dphidksi(6),dphideta(6)
real(kind=8) :: dldks1,dldks2,dldks3,dldet1,dldet2,dldet3
integer :: i
real(kind=8) :: a(3),b(3),c(3)
! coefficients for linear shape functions of reference element
data (a(i),i=1,3)/ 1d0, 0d0, 0d0/
data (b(i),i=1,3)/-1d0, 1d0, 0d0/
data (c(i),i=1,3)/-1d0, 0d0, 1d0/

do i=1,3
   phil(i) = a(i)+b(i)*xc+c(i)*yc
enddo
! quadratic shapefunctions
phiq(1) = phil(1)*(2*phil(1)-1)
phiq(3) = phil(2)*(2*phil(2)-1)
phiq(5) = phil(3)*(2*phil(3)-1)
phiq(2) = 4*phil(1)*phil(2)
phiq(4) = 4*phil(2)*phil(3)
phiq(6) = 4*phil(3)*phil(1)
! slope of linear shape functions on reference triangle
dldks1 = b(1)
dldks2 = b(2)
dldks3 = b(3)
dldet1 = c(1)
dldet2 = c(2)
dldet3 = c(3)
! derivatives of shape functions with respect to xi
dphidksi(1) = (4*phil(1)-1d0)*dldks1
dphidksi(3) = (4*phil(2)-1d0)*dldks2
dphidksi(5) = (4*phil(3)-1d0)*dldks3
dphidksi(2) = 4*(dldks1*phil(2) + dldks2*phil(1))
dphidksi(4) = 4*(dldks2*phil(3) + dldks3*phil(2))
dphidksi(6) = 4*(dldks3*phil(1) + dldks1*phil(3))
! and with respect to eta
dphideta(1) = (4*phil(1)-1d0)*dldet1
dphideta(3) = (4*phil(2)-1d0)*dldet2
dphideta(5) = (4*phil(3)-1d0)*dldet3
dphideta(2) = 4*(dldet1*phil(2) + dldet2*phil(1))
dphideta(4) = 4*(dldet2*phil(3) + dldet3*phil(2))
dphideta(6) = 4*(dldet3*phil(1) + dldet1*phil(3))

end subroutine getshape_quad_xi_eta

!   getshape_quad_xy
!   Get shape functions and derivatives of quadratic triangle
!   in point (xc,yc)
!
!   Parameters:
!       xc,yc       i    trial point
!       x,y         i    coordinates of nodal points
!       phil        o    linear shape functions  in (xc,yc)
!       phiq        o    quadratic shape functions in (xc,yc)
!       dphidksi    o    dphiq/dx in (xc,yc)
!       dphideta    o    dphiq/dy in (xc,yc)
!       nl          i    number of linear shapefunctions (3)
!       nq          i    number of quadratic shapefunctions (6)
!
!   PvK 013004
subroutine getshape_quad_xy(xc,yc,x,y,phil,phiq,dphidksi,dphideta,nl,nq)
implicit none
integer :: nl,nq
real(kind=8) :: xc,yc,suml
real(kind=8) :: x(nq),y(nq),phil(nl),phiq(nq),dphidksi(nq),dphideta(nq)
real(kind=8) :: a(3),b(3),c(3),xl(3),yl(3),area
real(kind=8) :: dldks1,dldks2,dldks3,dldet1,dldet2,dldet3
integer :: i

! get coefficients for linear shapefunctions
xl(1) = x(1)
yl(1) = y(1)
xl(2) = x(3)
yl(2) = y(3)
xl(3) = x(5)
yl(3) = y(5)
call trilin(xl,yl,a,b,c,area)
! linear shapefunctions 
suml = 0
do i=1,3
   phil(i) = a(i)+b(i)*xc+c(i)*yc
   suml = suml + phil(i)
enddo
! quadratic shapefunctions
phiq(1) = phil(1)*(2*phil(1)-1)
phiq(3) = phil(2)*(2*phil(2)-1)
phiq(5) = phil(3)*(2*phil(3)-1)
phiq(2) = 4*phil(1)*phil(2)
phiq(4) = 4*phil(2)*phil(3)
phiq(6) = 4*phil(3)*phil(1)
! slope of linear shape functions on reference triangle
dldks1 = b(1)
dldks2 = b(2)
dldks3 = b(3)
dldet1 = c(1)
dldet2 = c(2)
dldet3 = c(3)

! derivatives of shape functions with respect to xi
dphidksi(1) = (4*phil(1)-1d0)*dldks1
dphidksi(3) = (4*phil(2)-1d0)*dldks2
dphidksi(5) = (4*phil(3)-1d0)*dldks3
dphidksi(2) = 4*(dldks1*phil(2) + dldks2*phil(1))
dphidksi(4) = 4*(dldks2*phil(3) + dldks3*phil(2))
dphidksi(6) = 4*(dldks3*phil(1) + dldks1*phil(3))
! and with respect to eta
dphideta(1) = (4*phil(1)-1d0)*dldet1
dphideta(3) = (4*phil(2)-1d0)*dldet2
dphideta(5) = (4*phil(3)-1d0)*dldet3
dphideta(2) = 4*(dldet1*phil(2) + dldet2*phil(1))
dphideta(4) = 4*(dldet2*phil(3) + dldet3*phil(2))
dphideta(6) = 4*(dldet3*phil(1) + dldet1*phil(3))

end subroutine getshape_quad_xy

subroutine getshape7_xi_eta(xi,eta,phiq)
implicit none
real(kind=8) :: xi,eta,phiq(7)
real(kind=8) :: a(3),b(3),c(3),phil(3),phil123
integer :: i

!     *** coefficients for linear shape functions of reference element
data (a(i),i=1,3)/ 1d0, 0d0, 0d0/
data (b(i),i=1,3)/-1d0, 1d0, 0d0/
data (c(i),i=1,3)/-1d0, 0d0, 1d0/

do i=1,3
   phil(i) = a(i)+b(i)*xi+c(i)*eta
enddo
phil123 = phil(1)*phil(2)*phil(3)

!     *** shapefunctions of extended quadratic element
phiq(1) = phil(1)*(2*phil(1)-1) + 3*phil123
phiq(3) = phil(2)*(2*phil(2)-1) + 3*phil123
phiq(5) = phil(3)*(2*phil(3)-1) + 3*phil123
phiq(2) = 4*phil(1)*phil(2) - 12*phil123
phiq(4) = 4*phil(2)*phil(3) - 12*phil123
phiq(6) = 4*phil(3)*phil(1) - 12*phil123
phiq(7) = 27*phil123

return
end subroutine getshape7_xi_eta

subroutine getshape6_xi_eta(xi,eta,phiq)
implicit none
real(kind=8) :: xi,eta,phiq(6)
real(kind=8) :: a(3),b(3),c(3),phil(3)
integer :: i

! coefficients for linear shape functions of reference element
data (a(i),i=1,3)/ 1d0, 0d0, 0d0/
data (b(i),i=1,3)/-1d0, 1d0, 0d0/
data (c(i),i=1,3)/-1d0, 0d0, 1d0/

phil(1) = a(1)+b(1)*xi+c(1)*eta
phil(2) = a(2)+b(2)*xi+c(2)*eta
phil(3) = a(3)+b(3)*xi+c(3)*eta

! shapefunctions of quadratic element
phiq(1) = phil(1)*(2*phil(1)-1) 
phiq(3) = phil(2)*(2*phil(2)-1) 
phiq(5) = phil(3)*(2*phil(3)-1) 
phiq(2) = 4*phil(1)*phil(2) 
phiq(4) = 4*phil(2)*phil(3)
phiq(6) = 4*phil(3)*phil(1) 

end subroutine getshape6_xi_eta

    
! From Numerical Recipes
subroutine ludcmp(a,n,np,indx,d,ludcmp_error)
use sepmodulecomio
use control
implicit none
integer :: n,np,indx(n),NMAX
real(kind=8) :: d,a(np,np),TINY
parameter(NMAX=500,TINY=1.0e-20)
integer :: i,imax=0,j,k,ludcmp_error
real(kind=8) :: aamax,dum,sum,vv(NMAX)

!if (verbose.and.print_node) then
!   write(irefwr,*) 'a: ',a(1:np,1:np)
!   write(irefwr,*) 'indx: ',indx(1:n)
!   write(irefwr,*) 'd: ',d
!endif

ludcmp_error=0
d=1.
do 12 i=1,n
  aamax=0.
  do 11 j=1,n
    if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
  if (aamax.eq.0.) then
     if (print_node) write(irefwr,*) 'singular matrix in ludcmp'
     ludcmp_error=1
     return
  endif 
  vv(i)=1./aamax
12    continue
do 19 j=1,n
  do 14 i=1,j-1
    sum=a(i,j)
    do 13 k=1,i-1
      sum=sum-a(i,k)*a(k,j)
13        continue
    a(i,j)=sum
14      continue
  aamax=0.
  do 16 i=j,n
    sum=a(i,j)
    do 15 k=1,j-1
      sum=sum-a(i,k)*a(k,j)
15        continue
    a(i,j)=sum
    dum=vv(i)*abs(sum)
    if (dum.ge.aamax) then
      imax=i
      aamax=dum
    endif
16      continue
  if (j.ne.imax)then
    do 17 k=1,n
      dum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum
17        continue
    d=-d
    vv(imax)=vv(j)
  endif
  indx(j)=imax
  if(a(j,j).eq.0.)a(j,j)=TINY
  if(j.ne.n)then
    dum=1./a(j,j)
    do 18 i=j+1,n
      a(i,j)=a(i,j)*dum
18        continue
  endif
19    continue
return
end subroutine ludcmp
!  (C) Copr. 1986-92 Numerical Recipes Software '%1&9p#!.

subroutine lubksb(a,n,np,indx,b)
integer :: n,np,indx(n)
real(kind=8) :: a(np,np),b(n)
integer :: i,ii,j,ll
real(kind=8) :: sum
ii=0
do 12 i=1,n
  ll=indx(i)
  sum=b(ll)
  b(ll)=b(i)
  if (ii.ne.0)then
    do 11 j=ii,i-1
      sum=sum-a(i,j)*b(j)
11        continue
  else if (sum.ne.0.) then
    ii=i
  endif
  b(i)=sum
12    continue
do 14 i=n,1,-1
  sum=b(i)
  do 13 j=i+1,n
    sum=sum-a(i,j)*b(j)
13      continue
  b(i)=sum/a(i,i)
14    continue
return
end subroutine lubksb
!  (C) Copr. 1986-92 Numerical Recipes Software '%1&9p#!.
