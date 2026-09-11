subroutine determine_geometry(iuser,user)
use sepmodulecomio
use control
use coeff
use sepran_arrays
use geometry
use tracers
implicit none
integer ::  iuser(*) 
real(kind=8) :: user(*)
real(kind=8) :: volume_cyl,volint
integer :: i,ihelp

if (.not.print_node) pedebug=.false.
! Find volume and surface area for this geometry
iuser(2)=1
iuser(6)=7
iuser(7)=0
iuser(8)=icoorsystem
iuser(9)=0
iuser(10)=-6
 user(6)=1d0
if (pedebug) write(irefwr,*) 'determine_geometry: ',user(1)
volume = volint(0,1,1,kmesh1,kprob1,isol2,iuser,user,ihelp)
if (cyl) then
  if (axi) then
!   surf1 = r1, surf2=r2
    surf(1) = volume*3d0*r1*r1/(r2*r2*r2-r1*r1*r1)
    surf(2) = surf(1)*r2*r2/(r1*r1)
!   Use cylindrical volume to find frac (which represents an
!   relative angular fraction of the annulus, rather than a 
!   volume fraction).
    iuser(2)=1
    iuser(6)=7
    iuser(7)=0
!   Change coordinate system to Cartesian
    iuser(8)=0
    iuser(9)=0
    iuser(10)=-6
     user(6)=1d0
    volume_cyl=volint(0,1,1,kmesh1,kprob1,isol2,iuser,user,ihelp)
    frac = volume_cyl/(pi*r2*r2-pi*r1*r1)
!   This would be the volume fraction:
!   frac = volume/(4d0/3d0*pi*(r2*r2*r2-r1*r1*r1))
!   modify frac to be consistent with cylindrical definition
!   frac = 0.5*frac
!   (end volume fraction)

    if (frac.gt.0.5) then
!      Can't use more than half an annulus for axisymmetric coordinates
       if (print_node) then
          write(irefwr,*) 'PERROR(compr_start): frac > 0.5 and axi' 
          write(irefwr,*) '   cyl, axi   : ',cyl,axi
          write(irefwr,*) '   frac       : ',frac
        endif
       call instop
    endif
    if (print_node) then
      write(irefwr,'(''icoorsystem                    : '',2i10)') icoorsystem,icoor900
      write(irefwr,'(''volume                         : '',f10.4)') volume
      write(irefwr,'(''surface area                   : '',2f10.4)') surf(1),surf(2) 
      write(irefwr,'(''frac (0.5=full spherical shell): '',f10.4)') frac
    endif

  else

    surf(2) = volume*2d0*r2/(r2*r2-r1*r1)
    surf(1) = surf(2)*r1/r2
    frac = volume/(pi*r2*r2-pi*r1*r1)
    if (print_node) then
       write(irefwr,'(''icoorsystem            : '',2i10)') icoorsystem,icoor900
       write(irefwr,'(''volume                 : '',f8.4)') volume
       write(irefwr,'(''surface area           : '',2f8.4)') surf(1),surf(2)   
       write(irefwr,'(''frac (1=full cylinder) : '',f8.4)') frac
    endif
  endif
  dth_output = 2*frac*pi/nth_output
  dr_output  = (r2-r1)/nr_output

else

  ! Rectangular geometry; axisymmetric or Cartesian
  if (r2-r1.ne.1d0) then
     if (print_node) then
        write(irefwr,*) 'PERROR(compr_start): r2-r1 <> 1'
        write(irefwr,*) 'r2-r1: ',r2-r1
     endif
     call instop
  endif
  surf(1) = volume
  surf(2) = volume
  if (iclc.ne.4) then
     if (print_node) then
       write(irefwr,*) 'PERROR(compr_start): iclc should be 4'
       write(irefwr,*) 'for Cartesian models (cyl=.false.)'
     endif
     call instop
  endif
  ifilchoice=0
  do i=1,4
     ncurvn(i)=1
     icurv(i,1) = icloc(i)
  enddo
  call pefilxy(2)
  rlampix = (xcmax-xcmin)/(ycmax-ycmin)
  ! default perturbation has half wavelength of the width of the box
  if (wavel_perturb<0) wavel_perturb=rlampix
  frac=1d0
  if (print_node) then
     write(irefwr,'(''icoorsystem            : '',2i10)') icoorsystem, icoor900
     write(irefwr,'(''volume                 : '',f8.4)') volume
     write(irefwr,'(''surface area           : '',2f8.4)') surf(1),surf(2)   
     write(irefwr,'(''rlampix, dx_two_elem   : '',2f8.4)') rlampix,dx_two_elements
  endif
endif
surfbot=surf(1)
surftop=surf(2)
end subroutine determine_geometry
