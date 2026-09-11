!reference to `ieltogrid_'
!reference to `inittracheat_'
!reference to `inittrac_'
!reference to `tracerout_'
subroutine set_up_tracers()
use sepmodulecomio
use convparam
use coeff
use mpetrac
use tracers
use sepran_arrays
use mparallel
use geometry
use control
#ifdef MPI
use mpi
#endif
implicit none
integer :: iuserlc(100)
real(kind=8) :: userlc(1000)
!include 'SPcommon/ctimen'
integer :: ntot,ir,ierr
character(len=80) :: tfname
real(kind=8) :: ytop,ytop_d
real(kind=4) :: t01

if (print_node) write(irefwr,*) 'PINFO(set_up_tracers) itracoption = ',itracoption
! Create initial condition for markers/tracers
if (itracoption /= 0) then
   if (hannah) then
!     Special case: use one tracer to track top of plume
      itracoption=3 !!
      volume_per_tracer(1)=1e-3
      ntot=1
      ndist=1
      ntrac(1)=1
      if (restart) then
!        read information from restart file
         tfname='tracer.restart'
         inquire(file=tfname,exist=eos_exists)
         if (eos_exists) then
             open(9,file=tfname) 
!            The file also contains time - maybe should build in a consistency check
             read(9,*) tracer(1)%x,tracer(1)%y
             close(9)
         else
             if (print_node) then
              write(irefwr,*) 'PERROR(compr_start): file ',tfname
              write(irefwr,*) 'does not exist, but you should'
              write(irefwr,*) 'have a file with the tracer restart'
              write(irefwr,*) 'information'
             endif
             call instop
         endif
      else ! not restart
!        tracer on top of boundary layer with perturbation (see func)
         tracer(1)%x=0
         tracer(1)%y=plume_tracer_start
         if (print_node) then
           open(9,file='tracer.restart')
           write(9,*) tracer(1)%x,tracer(1)%y
           close(9)
         endif
      endif
   else ! .not.hannah
      if (cyl) then
        if (print_node) write(irefwr,'(''radius(min/max): '',2f10.7)') radius_min,radius_max
        do idist=1,ndist
          ! temporary reduction to JP
          if (x0tr(idist).lt.radius_min) x0tr(idist)=radius_min ! +1e-4
          if (x1tr(idist).gt.radius_max) x1tr(idist)=radius_max ! -1e-4
        enddo
      endif ! cyl
      ! always call inittrac as it allocates the arrays
      call inittrac()
      ! overwrite tracer coordinates upon restart
      if (restart) call readtrac(0)
   endif ! of .not.hannah  block
   call tracerout(0)
endif ! of itracoption /= 0 block

!#ifdef MPI
!call MPI_Barrier(MPI_COMM_WORLD,ierr)
!#endif
!if (print_node) write(6,*) 'stopping after inittrac/tracerout'
!call instop()

write(irefwr,*) 'itracoption=',itracoption
write(irefwr,*) 'ntrac(1): ',ntrac(1)
!write(irefwr,*) 'stop in set up tracers'
!call instop

! formerly in detrminmax: 
if (itracoption == 1 .or. itracoption == 3) then
!  keep center of tracers a bit away from the edge
!  This should be tied to the size of the tracer.
   if (volume_per_tracer(1) == 0d0) then
     write(irefwr,*) 'PERROR(detrminmax): volume_per_tracer is not set'
     write(irefwr,*) 'volume_per_tracer = ',volume_per_tracer
     call instop
   endif
   eps_tracer_bound(1)= sqrt(volume_per_tracer(1))
   eps_tracer_bound=0.0_8
else if (itracoption == 2) then
!  For the marker chain we should keep the boundary tracers
!  at the boundary
   eps_tracer_bound= 0d0
endif
eps_elem = eps_tracer_bound(1)



! tracerheat is an array that contains the heating per tracer
if (iqtype == 3 .or. iqtype == 5) then
    if (print_node) write(irefwr,*) 'PERROR(set_up_tracers):: iqtype=3 or 5 need work'
!   !call inittracheat(tracerheat,kmesh1,kprob1,isol2)
endif

if (itracoption == 2) then
   if (cyl) then
!    Initialize radial pixel grid
!    Number of pixels in radial direction
     nrpix = npix_radial
     drpix = (r2-r1)/nrpix
     dr_pixel = drpix
!    Number of pixels in middle of the cylinder r=(r1+r2)/2
     nthpix = nint(frac*pi*(r1+r2)*nrpix)
!    dxpix is the angular extent of each pixel
     if (print_node) then
      write(irefwr,*) '***********************************************'
      write(irefwr,*) 'Radial pixel grid specifications: '
      write(irefwr,*) 'Extent : ',nthpix,' by ',nrpix
      write(irefwr,*) 'Radial dimension of pixels: ',dr_pixel
     endif
     dthpix = frac*2*pi/nthpix
     do ir=1,nrpix
        dth_pixel(ir) = (r1+(ir-0.5)*drpix)*dthpix
     enddo
     if (print_node) then
      write(irefwr,*) 'Angular extent: ', dthpix
      write(irefwr,*) 'Which translates to an average length of: '
      do ir=1,nrpix,10
         write(irefwr,*) '     ',ir,dth_pixel(ir)
      enddo
      write(irefwr,*) '***********************************************'
     endif
 
   else 

     nypix = npix_radial
     nxpix = nint(nypix * rlampix)
     dxpix = rlampix/nxpix
     dypix = 1d0/nypix
     if (print_node) then
      write(irefwr,*) '**********************************************'
      write(irefwr,*) 'Cartesian pixel grid specifications: '
      write(irefwr,*) 'Extent : ',nxpix,' by ',nypix
      write(irefwr,*) 'Dimension of pixels : ',dxpix,dypix
      write(irefwr,*) '**********************************************'
     endif
     if (nxpix*nypix>NPIXMAX) then
        write(irefwr,*) 'specification of pixel grid is too large: '
        write(irefwr,*) 'nxpix*nypix=',nxpix*nypix
        write(irefwr,*) 'NPIXMAX=',NPIXMAX
        call instop
     endif
     ! output initial condition as pixel values
     call mardiv('markers.dat',-2)
     !write(irefwr,*) 'stop after pixel grid specification'
     !call instop

   endif

   if (ifollowchem.ne.0) then
      if (print_node) write(irefwr,*) 'PERROR(set_up_tracers): ifollowchem is obsolete'
      call instop
  
      nchem = 1
      if (restart) then
         !call readchem(coorreal,chemmark,nchem,'chem.start')
      else
         !call initchem(chemmark,kmesh1,kprob1,isol2,nchem,ncontzone,NTRACMAX)
      endif
   endif 
endif ! itracoption==2

if (cyl.and.itracoption/=0) then
!  Provide mapping of current mesh element numbers onto
!  an equidistant (r,theta) grid. This will make it easier
!  to update the position of the tracers.
   if (print_node) write(irefwr,'(''Create regular interpolation grid (ieltogrid): '')') 
   t00 = second()
   call ieltogrid()
   t01 = second()
   if (print_node) write(irefwr,*) 'ieltogrid took ',t01-t00,' seconds'
endif

! Test tracer interpolation scheme 
if (itracoption /= 0) call test_tracer_accuracy()
if (itracoption /= 0) tracers_are_set_up=.true.

end subroutine set_up_tracers
