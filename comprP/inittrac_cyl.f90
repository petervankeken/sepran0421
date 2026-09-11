!   INITTRAC_CYL
subroutine inittrac_cyl()
use sepmodulecomio
use convparam
use blobs
use mpetrac
use tracers
use brandenburg
use coeff
use geometry
use control
implicit none
integer :: nxtrac,nytrac,i,j,ip,ntot,ny_now
real(kind=8) ::  rlamtr,dx,dy,dxn,dyn,x,y,r,arclength,xm,ym,dyn_here
real(kind=8) ::  He3losszone(3),He4losszone(3),Ar40losszone(3)
real(kind=8) ::  xp,yp,yend,xa,xb,xbl,yb,theta,xtmin,xtmax,ytmin,ytmax
integer :: ip1,ix,iy,itrac,nchemmax
character(len=80) :: cally
save He3losszone,He4losszone,Ar40losszone
    

if (itracoption == 9) then
   ntot = 1
   ndist = 1
   ntrac(1) = 1
   tracer(1)%x = 0.5+r1
   tracer(2)%y = 0.0159221 - 0.5
   if (quart) then 
!     rotate 45 degrees counterclockwise
      xp = (tracer(1)%x-tracer(1)%y)*sqrt(0.5d0)
      yp = (tracer(1)%x+tracer(1)%y)*sqrt(0.5d0)
      tracer(1)%x = xp
      tracer(1)%y = yp
   endif
endif


if (itracoption == 1 .and. ibench_type/=101) then
   !  Uniform tracer distribution in [r0,r1]x[theta0,theta1]
   ! first figure out how many tracers we will have; allocate arrays; then
   ! set up tracer
   !open(99,file='tracers.000.asc')
   ntot = 0
   do idist=1,ndist
      dx = x1tr(idist)-x0tr(idist)
      dy = (y1tr(idist)-y0tr(idist))* 0.5*(x1tr(idist)+x0tr(idist))
      if (dx.eq.0.or.dy.eq.0) then
         if (print_node) then
            write(irefwr,*) 'PERROR(inittrac_cyl): dx or dy = 0'
            write(irefwr,*) '  dx = ',dx,' dy = ',dy
         endif
         call instop
      endif
      nxtrac = int(dx/dr)
      dxn = dx/nxtrac 
      volume_dist(idist)= (x1tr(idist)*x1tr(idist)-x0tr(idist)*x0tr(idist)) * (y1tr(idist)-y0tr(idist)) * 0.5
      if (print_node) then
         write(irefwr,'(''INITTRAC: Distribute particles uniformly'')')
         write(irefwr,'(''          between r     = '',f9.6, '' and '',f9.6)') x0tr(idist),x1tr(idist)
         write(irefwr,'(''          and     theta = '',f9.6, '' and '',f9.6)') y0tr(idist),y1tr(idist)
         write(irefwr,'(''          with spacing  dr = '',f9.5, '' ('',f9.1,'' km) '')') dxn,dxn*height_dim*1e-3
         write(irefwr,'(''INITTRAC: volume of this distribution ='', f9.5)') volume_dist(idist)
      endif

      ip  = ntot
      ip1 = ip
!     The calculation of the coordinates of the tracers.
!     (x,y) take the role of (r,theta). tracer%x,y is filled
!     with (r,theta)
      do i=1,nxtrac 
!        radius of current arc
         r = x0tr(idist) + (i-0.5)*dxn
         dyn = dxn/r
!        length of arc
         arclength = r*(y1tr(idist)-y0tr(idist))
         ny_now = int(arclength/dyn/r)
         ip = ip+ny_now
      enddo 
      ntrac(idist) = ip - ntot
      ntot = ntot + ntrac(idist)
      eps_tracer_bound(idist) = dxn/2
   enddo  ! idist=1,ndist

   if (duplicate_tracers) then
      if (ndist>1) then
         if (print_node) write(irefwr,*) 'PERROR(inittrac_cyl): duplicate_tracers but ndist>1: ',ndist
         call instop
      endif
      ntot=2*ntrac(1)
   endif
   if (print_node) write(irefwr,*) 'ntot = ',ntot,tracers_allocated

   if (tracers_allocated) then
      if (ntot>ntrac_allocated) then
         if (print_node) then
            write(irefwr,*) 'PERROR(inittrac_cyl): number of allocated tracers = ',ntrac_allocated
            write(irefwr,*) '              but required number of tracers ntot = ',ntot
         endif
         call instop
      endif
   else
      nchemmax=0
      cally='inittrac_cyl'
      call allocate_tracers(ntot,nchemmax,cally)
   endif


   ! now initialize tracer density
   if (print_node) write(irefwr,*) 'assign coordinates ',ndist
   if (ndist>1) then
      if (print_node) then
        write(irefwr,*) 'PERROR(inittrac_cyl): ndist should be 1 on input'
        write(irefwr,*) 'with neutral tracers creatd through duplicate_tracers'
      endif
      call instop
   endif
   ntot = 0
   do idist=1,ndist
      dx = x1tr(idist)-x0tr(idist)
      dy = (y1tr(idist)-y0tr(idist))* 0.5*(x1tr(idist)+x0tr(idist))
      if (dx.le.1e-12.or.dy.eq.1e-12) then
         if (print_node) then
            write(irefwr,*) 'PERROR(inittrac_cyl): dx or dy = 0'
            write(irefwr,*) '  dx = ',dx,' dy = ',dy
         endif
         call instop
      endif
      nxtrac = int(dx/dr)
      dxn = dx/nxtrac 

      ip = ntot
      ip1 = ip

!     The calculation of the coordinates of the tracers.
!     (x,y) take the role of (r,theta). tracer is filled
!     with (r,theta)
      xtmin=100
      ytmin=100
      xtmax=-100
      ytmax=-100
!     write(irefwr,*) 'x0tr,dxn: ',x0tr,dxn,r1,r2
      do i=1,nxtrac 
!        radius of current arc
         r = x0tr(idist) + (i-0.5)*dxn
         dyn = dxn/r
!        length of arc
         arclength = r*(y1tr(idist)-y0tr(idist))
         ny_now = int(arclength/dyn/r)
!        write(irefwr,'(''arclength: '',3f12.6,i10)') arclength,dyn, r,ny_now
         x = r
         ! dyn_here=dyn ! temporary restored for JP reproduction
         dyn_here = arclength/(r*ny_now) ! <<< this really is better
         do j=1,ny_now
            if (ip+j>ntrac_allocated) then
               if (print_node) write(irefwr,*) 'ip+j>ntrac_allocated: ',ip+j,ntrac_allocated
               call instop 
            endif
            y = y0tr(idist) + (j-0.5)*dyn_here
            xm = x*sin(y)
            ym = x*cos(y)
            call checkbounds_cyl(1,xm,ym)
            tracer(ip+j)%x=xm
            tracer(ip+j)%y=ym
            !write(99,*) xm,ym
         enddo
         ip = ip+ny_now
      enddo 
      ntrac(idist) = ip - ntot
      ntot = ntot + ntrac(idist)
      if (ntot.gt.NTRACMAX) then
         if (print_node) write(irefwr,*) 'PERROR(inittrac_cyl): ntot > NTRACMAX'
         call instop
      endif

!     scale density in tracers with volume of tracers
      volume_per_tracer(idist) = volume_dist(idist)/ntrac(idist)
!     modify eps_tracer_bound so that tracer centers can't get to the edge of the domain
      do i=1,ntrac(idist)
         tracer(ip1+i)%density = -Rb_local*volume_per_tracer(idist)
      enddo
      !close(99)

   enddo  ! idist=1,ndist

   if (duplicate_tracers) then
      if (print_node) write(irefwr,*) 'duplicate tracers'
      ! this is by definition -now- the neutral set of tracers
      if (ndist>1) then
         if (print_node) write(irefwr,*) 'PERROR(inittrac_cyl): duplicate_tracers but ndist>1: ',ndist
         call instop
      endif
      ndist=2
      ntrac(2)=ntrac(1)
      ntot=ntrac(1)+ntrac(2)
      if (print_node) write(irefwr,*) 'duplicate tracers',ntot,ntrac_allocated,ntrac(1),ntrac(2)
      do i=1,ntrac(1)
         tracer(ntrac(1)+i)%x=tracer(i)%x
         tracer(ntrac(1)+i)%y=tracer(i)%y
      enddo
      volume_per_tracer(2)=volume_per_tracer(1)
      eps_tracer_bound(2)=eps_tracer_bound(1)
      ! do not set density for idist=2
   endif
   if (print_node) then
      write(irefwr,'(''INITTRAC: ntrac: '',i7,2x,10i7,:)') ntot,(ntrac(idist),idist=1,ndist)
      write(irefwr,'(''INITTRAC: for first set of tracers: '')') 
      write(irefwr,'(''INITTRAC: x and y ranges: '',5f12.5)') minval(tracer(1:ntot)%x), maxval(tracer(1:ntot)%x), &
          & minval(tracer(1:ntot)%x), maxval(tracer(1:ntot)%y),eps_tracer_bound(1)
   endif
   if (its_CH94 .or. gable_plates) then
      tracer(1:ntot)%meltable=.true.
      tracer(1:ntot)%ingasable=.true.
      tracer(1:ntot)%in_cmb=.false.
   endif

  else if (itracoption==1.and.ibench_type==101) then
 
    !! "Stokes sphere" benchmark

    !! first compute number of tracers
    itrac=0
    do iblob=1,nblob
       nypix=int((2*r_blob(iblob)-dr)/dr+1)
       do iy=nypix/2,-nypix/2,-1 ! start at the top
          yb=y_blob(iblob)+iy*dr ! global coordinates
          theta=acos(iy*dr/r_blob(iblob))
          xbl=r_blob(iblob)*sin(theta)
          nxpix=nint(2*xbl/dr) + 1 ! one in center
!         write(irefwr,*) 'blobs: ',yb,xbl,nxpix,theta,ip+itrac
          itrac=itrac+2*(nxpix/2)+1
       enddo
    enddo

    if (print_node) write(irefwr,*) 'call allocate tracers: ',itrac
    ntrac_allocated=itrac
    call allocate_tracers(ntrac_allocated,1,'inittrac')

    icenterblob=0
    ip=0
    do iblob=1,nblob
       itrac=0
       ! one tracer is at the center of the blob and has local
       ! coordinates (centerpoint) (0,0).
       nypix=int((2*r_blob(iblob)-dr)/dr+1)
       do iy=nypix/2,-nypix/2,-1 ! start at the top
          yb=y_blob(iblob)+iy*dr ! global coordinates
          theta=acos(iy*dr/r_blob(iblob))
          xbl=r_blob(iblob)*sin(theta)
          nxpix=nint(2*xbl/dr) + 1 ! one in center
!         write(irefwr,*) 'blobs: ',yb,xbl,nxpix,theta,ip+itrac
          do ix=-nxpix/2,nxpix/2
             xb=x_blob(iblob)+ix*dr
             itrac=itrac+1
             if (ip+itrac>ntrac_allocated) then
                if (print_node) then
                   write(irefwr,*) 'error: ip+itrac>ntrac_allocated'
                   write(irefwr,*) iblob,ip+itrac,ntrac_allocated
                endif
                call instop
             endif
             tracer(ip+itrac)%x=xb
             tracer(ip+itrac)%y=yb
             if (ix==0.and.iy==0) icenterblob(iblob)=ip+itrac
          enddo
!         if (nxpix==1) write(irefwr,*) 'xb,yb: ',xb,yb
       enddo
       ntrac(iblob)=itrac
       volume_per_tracer(iblob)=pi*r_blob(iblob)*r_blob(iblob)/ ntrac(iblob) 
       if (print_node) write(irefwr,*) 'fill density',iblob,ip+1,ip+ntrac(iblob)
       do i=1,ntrac(iblob)
          tracer(ip+i)%density=-Rb_local*delta_rhop(iblob)* volume_per_tracer(iblob)
       enddo
       ip=ip+ntrac(iblob)
    enddo
    if (print_node) then
       open(9,file='blobs.dat')
       do i=1,sum(ntrac(1:nblob))
         write(9,'(3e15.7)') tracer(i)%x,tracer(i)%y,tracer(i)%density
       enddo
       close(9)
       open(9,file='centerblobs.dat') 
       do iblob=1,nblob
         write(9,'(3i10,3e15.7)') iblob,icenterblob(iblob),ntrac(1), tracer(icenterblob(iblob))%x-x_blob(iblob), &
          &      tracer(icenterblob(iblob))%y-y_blob(iblob), volume_per_tracer(iblob)
       enddo
       close(9)
    endif
    ndist=1
    ntrac(1)=sum(ntrac(1:nblob)) ! collect all blobs as active tracers into first distribution
    ntot=ntrac(1)
endif ! itractoption == 1

if (itracoption == 2) then
!  markerchain at r=r1
   if (ndist.gt.1) then
      if (print_node) then
         write(irefwr,*) 'PERROR(cylaxm/inittrac_cyl): itracoption=2' 
         write(irefwr,*) '  but ndist > 1. The markerchain method works' 
         write(irefwr,*) '  only for one markerchain right now.'
      endif
      call instop
   endif
!  Note: x is radius r and y is angle theta here
   x  = y0m(1)
!  dy is length of markerchain 
   dy = 2*pi*frac*x
   if (print_node) write(irefwr,*) 'PINFO(inittrac_cyl): radius/length of chain: ',x,dy 
   if (dy<1e-12) then
       if (print_node) then
          write(irefwr,*) 'PERROR(inittrac_cyl): dy = 0'
          write(irefwr,*) 'itracoption = 2'
          write(irefwr,*) 'y0m = ',y0m
          write(irefwr,*) 'frac = ',frac
       endif
       call instop
   endif
   if (dm(1)==0) then
      if (print_node) then
         write(irefwr,*) 'PERROR(cylaxm/inittrac_cyl): itracoption=2'
         write(irefwr,*) 'but dm=0; set dr to desired dm'
      endif
      call instop
   endif
!  ntrac = Number of tracers in this chain
   ntrac(1) = int(dy/dm(1))
   imark(1)=ntrac(1)
   dyn = dy/(x*(ntrac(1)-1))
   if (print_node) then 
      write(irefwr,*) 'PINFO(inittrac_cyl): ntrac/imark: ', ntrac(1),imark(1)
      write(irefwr,*) 'PINFO(inittrac_cyl): ainit(1): ',ainit(1),dyn
   endif
!  dyn = angular separation between tracers
   if (ntrac(1).gt.NTRACMAX) then
      if (print_node) then
        write(irefwr,*) 'PERROR(inittrac_cyl): ntrac > NTRACMAX'
        write(irefwr,*) '  ntrac    = ',ntrac
        write(irefwr,*) '  NTRACMAX = ',NTRACMAX
      endif
      call instop
   endif
!  start chain at theta=0
   yend = (ntrac(1)-1)*dyn
   do j=1,ntrac(1)
      y = (j-1)*dyn
      xa = x+ainit(1)*cos(pi*y/yend)
      xm = xa*sin(y)
      ym = xa*cos(y)
      call checkbounds(1,xm,ym)
      tracer(j)%x=xm
      tracer(j)%y=ym
   enddo
   ntot = 0
   imark(1)=ntrac(1)
   do ichain=1,nochain
      ntot = ntot + imark(ichain)
   enddo
   if (print_node) then
      write(irefwr,*) 'PINFO(inittrac_cyl): nmark = ',ntot
      write(irefwr,*) 'PINFO(inittrac_cyl): first marker: ',tracer(1)%x,tracer(1)%y
      write(irefwr,*) 'PINFO(inittrac_cyl): last marker : ', tracer(ntot)%x,tracer(ntot)%y
   endif
endif
if (imark(1).gt.NTRACMAX) then
   if (print_node) then
      write(irefwr,*) 'PERROR(inittrac_cyl): ntot > NTRACMAX' 
      write(irefwr,*) 'ntot, NTRACMAX: ',imark(1),NTRACMAX
   endif
   call instop
endif
   
end subroutine inittrac_cyl
