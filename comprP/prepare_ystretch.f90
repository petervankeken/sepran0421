subroutine prepare_ystretch(htrac,zstretch,ystretch,nytractot,dytrac)
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
real(kind=8) :: htrac(*),zstretch(*),ystretch(*)


if (compress.and.stretch_tracers) then
   if (.not.fake_rho_bar.or.step_rho_background.or.exp_rho_background) then
      ! find distances between tracers that decrease as density increases
      htrac(1)=dytrac
      zstretch(1)=htrac(1)/2
      massgoal=dytrac ! assume rho=1 at top...
      ! make first estimate of layerthicknesses (should decrease with depth)
      htot=htrac(1)
      i=0
      do
         i=i+1
         if (i>nytractot) then
            if (print_node) write(irefwr,*) 'PERROR(prepare_ystretch): estimate for nytractot is too small: ',nytractot
            if (print_node) write(irefwr,*) 'i, zstretch(i): ',i-1,zstretch(i-1),htrac(i-1)
            call instop
         endif
         zstretch(i+1)=zstretch(i)+htrac(i)
         htrac(i+1)=htrac(i)
         mass_now=funccf(3,0.0_8,1.0_8-zstretch(i+1),0.0_8)*htrac(i+1)
         ! correct thickness of layer
         htrac(i+1)=htrac(i)*massgoal/mass_now
         ! reassign zstretch
         zstretch(i+1)=zstretch(i)+htrac(i+1)
         htot=htot+htrac(i+1)
         if (htot>1.0_8) exit
       enddo
       nyhere=i+1
       ! now restretch it a bit to make an integer number of tracers fit
       htrac(1:nytractot)=htrac(1:nytractot)/htot
       zstretch(1)=0.5*htrac(1)
       ystretch(nyhere)=1.0_8-zstretch(1)
       do i=2,nyhere
          zstretch(i)=zstretch(i-1)+htrac(i)
          ystretch(nyhere-i+1)=1.0_8-zstretch(i)
       enddo
       if (print_node) then
          open(9,file='ystretch.dat')
          write(9,*) 0,ystretch(1)
          do i=2,nyhere
             write(9,*) 0,ystretch(i),1.0_8/(ystretch(i)-ystretch(i-1))
          enddo
          close(9)
       endif
       !write(irefwr,*) 'stop in inittrac_cart'
       !call instop
       
   else ! .not.(.not.fake_rho_bar.or.step_rho_background.or.exp_rho_background)
       density_ratio=1+drho_background_dense
       ! first make a regular set of tracer centers in 1D that stretch out along the density profile
        ! we're assuming particle distance for top most particle h_n and that of bottom most particle h_1 are
        ! at a ratio dictated by the density ratio: h_n = density_ratio*h1. Sum h_i = 1. 

        ! if density=1+z*drho_background_dense then the ratio between two adjacent
        ! h's is a: h_i=a*h_i+1. Then h_n=a^n*h_1 or a=density_ratio^1/n. Or a=10^(1/n*log10(density_ratio))
        factor_a=10**(dytrac*log10(density_ratio))
        ! set trial h1 at top
        htrac(1)=1.0_8
        do i=2,nytractot
           htrac(i)=htrac(i-1)*factor_a
        enddo
        if (print_node) write(irefwr,*) 'htrac1,N: ',htrac(1),htrac(nytractot)
        if (print_node) write(irefwr,*) 'total length of htrac: ',sum(htrac(1:nytractot))
        scale_sumhi=1.0_8/sum(htrac(1:nytractot))
        htrac(1:nytractot)=htrac(1:nytractot)*scale_sumhi
        if (print_node) write(irefwr,*) 'total length of htrac: ',sum(htrac(1:nytractot))
        ! now form the array of center locations 
        ystretch(1)=0.5*htrac(1)
        htop=htrac(1)
        do i=2,nytractot
           ystretch(i)=htop+0.5*htrac(i-1)
           htop=htop+htrac(i)
        enddo
        if (print_node) then
           open(9,file='ystretch.dat')
           do i=1,nytractot
              write(9,*) 0.,ystretch(i)
           enddo
           nyhere=nytractot
           close(9)
        endif
   endif! loop creating array ystretch
else ! compress .and. stretch_tracers
   ! classical case of unstretch tracers
   htrac(1:nytractot)=dytrac
   do i=1,nytractot
     ystretch(i)=0.5*dytrac+(i-1)*dytrac
     zstretch(i)=1.0_8-ystretch(i)
   enddo
endif  ! compress .and. stretch_tracers


end subroutine prepare_ystretch
