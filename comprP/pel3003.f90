! PEL3003
!  
! Modified from el3003.  Intended for buoyancy given by tracers
! in primitive variable formulation and Taylorhood
!  
! called from elm900rhsd_pvk.f that replaces elm900rhsd from sepran library
! ichoice=1 : f1
! ichoice=2 : f2
!  
! PvK 990408
! Updated for version 0220 09222020
subroutine pel3003(g,f,w,phi,mdecl,ndecl,work,x,ichoice,xgauss)
use sepmodulecomio
use tracers
use control
use coeff
use geometry
implicit none
integer mdecl,ndecl,ichoice,i
real(kind=8) :: g(mdecl),f(ndecl),w(mdecl),phi(ndecl,mdecl)
real(kind=8) :: work(mdecl),x(ndecl,*),xgauss(*)

if (itracoption == 1 .and. tracerC) then
   ! Tracer method using buoyancy in tracers
   if (cyl) then
      call pel3003_cylt(g,f,w,phi,mdecl,ndecl,work,ichoice)
   else
      if (ichoice==2) then
        ! fill only for upward buoyancy force
        call pel3003_cartt(g,f,w,phi,mdecl,ndecl,work,ichoice)
      endif
   endif

else if (itracoption == 2) then
   ! Markerchain method

   if (cyl) then
      ! call pel3003_cylm(g,f,w,phi,mdecl,ndecl,work,x,ichoice)
      if (print_node) then 
        write(irefwr,*) 'PERROR(pel3003): needs update for markerchain method & cyl'
      endif
      call instop
   else if (ichoice.eq.2) then
        ! fill only for upward buoyancy force
       call pel3003_cartm(g,f,w,phi,mdecl,ndecl,work,x)
   endif

endif

end subroutine pel3003


! PEL3003_CART_T
! 
! Add tracer buoyancy contribution (prepared by tracrhsd
! and stored in /petrac/ f_tracrhsd
!  
! ichoice=1 : f1
! ichoice=2 : f2
! PvK 111505
subroutine pel3003_cartt( g,f,w,phi,mdecl,ndecl,work,ichoice)
use mpetrac
use tracers
use coeff
use sepmodulecactl
implicit none
integer mdecl, ndecl, ichoice
real(kind=8) :: g(mdecl), f(ndecl), w(mdecl), phi(ndecl,mdecl),work(mdecl)
integer :: i

! first buoyancy component
if (jadd.eq.0) then
   do i=1,n
      f(i) = 0d0
   enddo
endif

! second buoyancy component
do i=1,n
   f(i) = f(i) + f2_tracrhsd(i,ielem,1)  ! always use full (first part of) this buffer
enddo
!if ((ielem/100)*100==ielem) write(irefwr,'(''f2: '',2i5,7e15.3:)') ielem,n,f(1:n)
!write(irefwr,*) 'in pel3003_cartt'
!call instop

end subroutine pel3003_cartt

! PEL3003_CYLT
!  
! Add tracer buoyancy contribution (prepared by tracrhsd
! and stored in /petrac/ f_tracrhsd
!  
!     ichoice=1 : f1
!     ichoice=2 : f2
! PvK 111505
subroutine pel3003_cylt( g,f,w,phi,mdecl,ndecl,work,ichoice)
use mpetrac
use tracers
use coeff
use sepmodulecactl
implicit none
integer mdecl, ndecl, ichoice
real(kind=8) :: g(mdecl), f(ndecl), w(mdecl), phi(ndecl,mdecl), work(mdecl)
integer :: i

if (jadd.eq.0) then
   do i=1,n
      f(i) = 0d0
   enddo
endif

if (ichoice==1) then
   ! horizontal buoyancy
   do i=1,n
      f(i) = f(i) + f1_tracrhsd(i,ielem,1)  ! always use full (first part of) this buffer
   enddo
else
   do i=1,n
      f(i) = f(i) + f2_tracrhsd(i,ielem,1)  ! always use full (first part of) this buffer
   enddo
endif
end subroutine pel3003_cylt

! PEL3003_CART
! 
! Modified from el3003.  Intended for buoyancy given by tracers
! in primitive variable formulation.
! 
! Note: neutral layer should be indicated in pixel grid by 
! the value '0'; the compositionally distinct layer by 1.
! Note that values different than 1 are ignored in this version.
! 
! Modified for the markerchain method
!  
! PvK 990408
! Update to sepran 1020 PvK 20201013
subroutine pel3003_cartm( g,f,w,phi,mdecl,ndecl,work,x)
use sepmodulecomio
use tracers
use coeff
use msper01
use sepmodulecactl
use control
implicit none
integer mdecl, ndecl
real(kind=8) :: g(mdecl),f(ndecl),w(mdecl),phi(ndecl,mdecl),work(mdecl),x(ndecl,*)

real(kind=8) :: h, sum, xn(7),yn(7),shapef(7),xm,ym,rl(3),xp,yp
integer :: i, k,itrac,pefirst=0
integer :: ip,ipx,ipy,ipxmin,ipxmax,ipymin,ipymax
integer :: nodlin(3),nodno(6),isub
real(kind=8) :: xcmin,xcmax,ycmin,ycmax,funccf
logical :: out
real(kind=8),parameter ::  epsh=1.e-8_8
save pefirst


if (jadd.eq.0) then
   do i=1,n
      f(i) = 0d0
   enddo
endif

if (jdiag.eq.1) then
   write(irefwr,*) 'PERROR(pel3003): '
   write(irefwr,*) '      not suited for Newton-Cotes integration'
   call instop
endif
if (jconsf.eq.1) then
   write(irefwr,*) 'PERROR(pel3003): not suited for jconsf=1'
   call instop
endif

do i=1,ndecl
   xn(i) = x(i,1)
   yn(i) = x(i,2)
enddo

! Assume we have a triangle with straight edges here
! (ok for cartesian box). Determine smallest rectangular
! area around the element
xcmin = min(xn(1),xn(3),xn(5))
ycmin = min(yn(1),yn(3),yn(5))
xcmax = max(xn(1),xn(3),xn(5))
ycmax = max(yn(1),yn(3),yn(5))
ipxmin = xcmin/dxpix + 1.0
ipymin = ycmin/dypix + 1.0
ipxmax = xcmax/dxpix + 1.5
ipymax = ycmax/dypix + 1.5
ipxmin = max(1,ipxmin)
ipymin = max(1,ipymin)
ipxmax = min(nxpix,ipxmax)
ipymax = min(nypix,ipymax)

! loop over pixels in this area; 
!   if 1) pixel is in heavy layer and
!      2) pixel is in element then
!         compute shapefunction and add to vector
do ipy=ipymin,ipymax
   yp = (ipy-0.5)*dypix
   do ipx=ipxmin,ipxmax
      xp = (ipx-0.5)*dxpix
      ip = (ipy-1)*nxpix + ipx
      if (pefirst.eq.0) write(irefwr,*) 'xp,yp: ',xp,yp,ip,ipix(ip)
      if (ipix(ip).eq.1) then
!         we are in the compositionally distinct layer
!         add buoyancy contribution
          call findrl(xn,yn,xp,yp,rl,isub,nodno,nodlin)
          if (pefirst.eq.0) write(irefwr,*) '   : ',xp,yp,rl(1),rl(2),rl(3)
          out = (rl(1).ge.-eps).and.(rl(2).ge.-eps).and.(rl(3).ge.-eps)
          if (out) then
             if (pefirst.eq.0) write(irefwr,*) 'in triangle: ',ipx,ipy
!            calculate the shape function in (xp,yp)
             call detshape7(xn,yn,xp,yp,shapef)
             if (compress) then 
                do i=1,7
                   f(i) = f(i) - funccf(3,xp,yp,yp)*Rb_local*shapef(i)*dxpix*dypix
                enddo
             else
                do i=1,7
                   f(i) = f(i) - Rb_local*shapef(i)*dxpix*dypix
                enddo
             endif
          endif
      endif
   enddo
enddo
if (pefirst==0) write(irefwr,'(''f: '',7e15.7)') f(1:7)
pefirst=1
end subroutine pel3003_cartm
