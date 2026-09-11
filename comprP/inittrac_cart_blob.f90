subroutine inittrac_cart_blob()
use sepmodulecomio
use convparam
use blobs
use mpetrac
use tracers
use coeff
use mparallel
use geometry
use control
#ifdef MPI
use mpi
#endif
implicit none
integer :: ip,i,j,ntot,iy,ix,itrac,nxtrac,nytrac,ip1
real(kind=8) ::  dxtrac,dytrac,x,y,local_volume
real(kind=8) ::  piw,ym,xm,yoffset,ymin,ymax,row_intvl,col_intvl
real(kind=8) ::  dx,xstart,xend,xhere,xb,yb,theta,xbl,dz
character(len=80) :: tname,fdummy
integer :: nxtrac_here,nyhere,ierr=0
real(kind=8) :: density_ratio,factor_a,sumhi,scale_sumhi,htop,mass_now,massgoal,htot,funccf,stretchcorrect
integer :: nytractot
real(kind=8),dimension(:),allocatable :: htrac,ystretch,zstretch
logical my_turn

if (ratio_method) then

   ! call blobs_init_ratiomethod

else ! not ratio method; just discretize spheres.
   !! first compute number of tracers
   itrac=0
   do iblob=1,nblob
      nypix=int((2*r_blob(iblob)-dr)/dr+1) ! 
      do iy=nypix/2,-nypix/2,-1 ! start at the top
         yb=y_blob(iblob)+iy*dr ! global coordinates
         theta=acos(iy*dr/r_blob(iblob))
         xbl=r_blob(iblob)*sin(theta)
         nxpix=int(2*xbl/dr) + 1 ! one in center
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
      nypix=int((2*r_blob(iblob)-dr)/dr+1)  ! 
      do iy=nypix/2,-nypix/2,-1 ! start at the top
         yb=y_blob(iblob)+iy*dr ! global coordinates
         theta=acos(iy*dr/r_blob(iblob))
         xbl=r_blob(iblob)*sin(theta)
         nxpix=int(2*xbl/dr) + 1 ! one in center
!              write(irefwr,*) 'blobs: ',yb,xbl,nxpix,theta,ip+itrac
         do ix=-nxpix/2,nxpix/2
            xb=x_blob(iblob)+ix*dr
            itrac=itrac+1
            if (ip+itrac>ntrac_allocated) then
               if (print_node) write(irefwr,*) 'error: ip+itrac>ntrac_allocated'
               if (print_node) write(irefwr,*) iblob,ip+itrac,ntrac_allocated
               call instop
            endif
            tracer(ip+itrac)%x=xb
            tracer(ip+itrac)%y=yb
            if (ix==0.and.iy==0) icenterblob(iblob)=ip+itrac
         enddo
!              if (nxpix==1) write(irefwr,*) 'xb,yb: ',xb,yb
      enddo
      ntrac(iblob)=itrac
      volume_per_tracer(iblob)=pi*r_blob(iblob)*r_blob(iblob)/ ntrac(iblob) 
      if (print_node) write(irefwr,*) 'fill tracer density',iblob,ip+1,ip+ntrac(iblob)
      do i=1,ntrac(iblob)
         tracer(ip+i)%density=-Rb_local*delta_rhop(iblob)*volume_per_tracer(iblob)
      enddo
      ip=ip+ntrac(iblob)
   enddo
   if (print_node) then
      open(9,file='blobs.dat')
      do itrac=1,sum(ntrac(1:nblob))
          write(9,'(3e15.7)') tracer(itrac)%x,tracer(itrac)%y, tracer(itrac)%density
      enddo
      open(9,file='centerblobs.dat') 
      do iblob=1,nblob
        write(9,'(3i10,3e15.7)') iblob,icenterblob(iblob), ntrac(1),tracer(icenterblob(iblob))%x-x_blob(iblob), &
          &      tracer(icenterblob(iblob))%y-y_blob(iblob),volume_per_tracer(iblob)
      enddo
      close(9)
      ndist=1
      ntrac(1)=sum(ntrac(1:nblob)) ! collect all blobs as active tracers into first distribution
      ntot=ntrac(1)
    endif
 endif ! not.ratio_method

end subroutine inittrac_cart_blob
