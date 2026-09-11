!   Prepare lookup table /petracer/ f_tracrhsd
!   for use in the element rhsd construction
subroutine tracrhsd()
use sepmodulecomio
use sepmodulekmesh
use mparallel
use mpetrac
use msper01
use sepran_arrays
use sepran_interface
use geometry
use coeff
use control
use tracers
#ifdef MPI
use mpi
#endif
implicit none
integer :: ierr

if (.not.tracers_are_set_up) then
   if (print_node) write(irefwr,*) 'PWARN(tracrhsd): called while tracers are not yet set up'
   return
endif

if (mpi_partrac.and. numprocs>1) then
   i_f_tracrhsd = 2  ! use second part of f12_tracrhsd for partial buffer
   if (cyl) f1_tracrhsd(1:7,1:nelem,1)=0d0 ! clear part of f1_tracrhsd that will be used for reduction
   f2_tracrhsd(1:7,1:nelem,1)=0d0 ! clear part of f2_tracrhsd that will be used for reduction
else
   i_f_tracrhsd = 1 ! business as usual; ignore second part of f_tracrhsd
endif


if (cyl.and.itype_stokes==903) then
   if (print_node) then
     write(irefwr,*) 'PERROR(tracrhsd): need to update for cyl and 903'
   endif
   call instop()
endif
if (itracoption==1) then
   ! tracer Stokeslet method
   call tracrhsd01()
else
   ! markerchain method
   call tracrhsd02()
endif

#ifdef MPI
if (mpi_partrac.and. numprocs > 1) then
   ! add partial buffers of f_tracrhsd and reduce to first part
   if (cyl) then
     ! horizontal component for cyl only
     call MPI_REDUCE(f1_tracrhsd(1,1,2),f1_tracrhsd(1,1,1),7*nelem, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)
     if (ierr/=0) then
         write(irefwr,*) 'PERROR(trac_rhsd_stokes) MPI_REDUCE:',ierr
         write(irefwr,*) 'myid = ',myid
         call instop
     endif
   endif
   !write(6,'(''reduce: '',i5,2e15.7,i10,2i5)') myid,f2_tracrhsd(1,1,1:2),7*nelem,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD
   ! always reduce vertical component
   call MPI_REDUCE(f2_tracrhsd(1,1,2),f2_tracrhsd(1,1,1),7*nelem,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   if (ierr/=0) then
       write(irefwr,*) 'PERROR(trac_rhsd_stokes) MPI_REDUCE:',ierr
       write(irefwr,*) 'myid = ',myid
       call instop
   endif

   ! now broadcast f2 and f1 (if necessary)
   if (cyl)  call MPI_BCAST(f1_tracrhsd(1,1,1),7*nelem, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(f2_tracrhsd(1,1,1),7*nelem, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
   if (ierr/=0) then
      write(irefwr,*) 'PERROR(trac_rhsd_stokes) MPI_BCAST:',ierr
      write(irefwr,*) 'myid = ',myid
      call instop
   endif
endif
#endif

! check some values in f2_tracrhsd
!if (print_node) then
!   write(irefwr,'(''f2t: '',i5,3e15.7)') myid,f2_tracrhsd(1,1,1),minval(f2_tracrhsd(1:7,1:nelem,1)), &
!        & maxval(f2_tracrhsd(1:7,1:nelem,1))
!endif
!if (mpi_partrac.and.numprocs>1) then
!   write(99+myid,'(''f2t: '',i5,7e15.7)') myid,f2_tracrhsd(1:7,1,2)
!else
!   write(98,'(''f2t: '',i5,7e15.7)') myid,f2_tracrhsd(1:7,1,1)
!endif
!call instop

end subroutine tracrhsd

subroutine tracrhsd01()
use sepmodulecomio
use sepmodulekmesh
use mparallel
use mpetrac
use tracers
use coeff
use geometry
use sepran_arrays
use sepran_interface
use control
implicit none
    
integer :: ntot,i,j,k,ielh
real(kind=8) :: xm,ym,xn(6),yn(6),shapef(7),xi,eta,phiq(7)
real(kind=8) :: chbeta,funccf,rl(3)
integer :: nodno(6),nodlin(3),itrac_start,itrac_end
integer :: itrac,isub,iharz
logical :: fail,correct

if (nelem.gt.MAXELEM) then
   if (print_node) write(irefwr,*) 'PERROR(tracrhsd_cart): nelem > MAXELEM'
   if (print_node) write(irefwr,*) 'nelem, MAXELEM: ',nelem,MAXELEM
   call instop
endif
 
itrac_start = 1
itrac_end = ntrac(1) ! assume all active tracers are in first distribution

if (.not.cyl) call pefilxy(2)

! clear current buffer
if (cyl) f1_tracrhsd(1:7,1:nelem,i_f_tracrhsd)=0.0_8
f2_tracrhsd(1:7,1:nelem,i_f_tracrhsd)=0.0_8

if (ibuoy_trac/= 1 .and. ibuoy_trac/= 2) then
  if (print_node) write(irefwr,*) 'PERROR(tracrhsd_cartt: wrong ibuoy_trac= ',ibuoy_trac
  call instop
endif

if (cyl) then
  if (itype_stokes==900) then
     call tracrhsd01_900_cyl(kmesh1,itrac_start,itrac_end)
  else
    if (print_node) write(irefwr,*) 'tracrhsd01 needs updating for 903 and cyl'
   !   call tracrhsd01_903_cyl(kmesh1,itrac_start,itrac_end)
   !       & f2_tracrhsd(1,1,i_f_tracrhsd),itrac_start,itrac_end)
  endif

 
else ! not cylindrical

  if (itype_stokes==900) then
    call tracrhsd01_900_cart(kmesh1,itrac_start,itrac_end)
  else
    !write(irefwr,*) 'PERROR(tracrhsd01): needs updating for Cartesian 903'
    !call instop
    call tracrhsd01_903_cart(kmesh1,itrac_start,itrac_end)
  endif
endif !cyl
end subroutine tracrhsd01

subroutine tracrhsd02()
use sepmodulecomio
use sepmodulekmesh
use mparallel
use mpetrac
use tracers
use coeff
use geometry
use sepran_arrays
use sepran_interface
use control
implicit none
    
integer :: ntot,i,j,k,ielh
real(kind=8) :: xm,ym,xn(6),yn(6),shapef(7),xi,eta,phiq(7)
real(kind=8) :: chbeta,funccf,rl(3)
integer :: nodno(6),nodlin(3),itrac_start,itrac_end
integer :: itrac,isub,iharz
logical :: fail,correct

if (nelem.gt.MAXELEM) then
   if (print_node) write(irefwr,*) 'PERROR(tracrhsd_cart): nelem > MAXELEM'
   if (print_node) write(irefwr,*) 'nelem, MAXELEM: ',nelem,MAXELEM
   call instop
endif
 
if (.not.cyl) call pefilxy(2)

! clear current buffer
if (cyl) f1_tracrhsd(1:7,1:nelem,i_f_tracrhsd)=0.0_8
f2_tracrhsd(1:7,1:nelem,i_f_tracrhsd)=0.0_8

if (ibuoy_trac/= 1 .and. ibuoy_trac/= 2) then
  if (print_node) write(irefwr,*) 'PERROR(tracrhsd_cartt: wrong ibuoy_trac= ',ibuoy_trac
  call instop
endif

if (cyl) then
  if (print_node) write(irefwr,*) 'PERROR(tracrhsd02): not yet set up for cyl'
  call instop
  !if (itype_stokes==900) then
  !   call tracrhsd01_900_cyl(kmesh1,itrac_start,itrac_end)
  !else
  !  write(irefwr,*) 'tracrhsd01 needs updating for 903 and cyl'
  ! !   call tracrhsd01_903_cyl(kmesh1,itrac_start,itrac_end)
  ! !       & f2_tracrhsd(1,1,i_f_tracrhsd),itrac_start,itrac_end)
  !endif
 
else ! not cylindrical

  if (itype_stokes==900) then
    call tracrhsd02_900_cart()
  else
    if (print_node) write(irefwr,*) 'PERROR(tracrhsd01): needs updating for Cartesian 903'
    call instop
    !call tracrhsd02_903_cart(kmesh1,itrac_start,itrac_end)
  endif
endif !cyl
end subroutine tracrhsd02

! tracrhsd02_900_cart
! compute element contributions to load vector for markerchain method (Cartesian, element 900)
! PvK 20201013
subroutine tracrhsd02_900_cart()
use sepmodulekmesh
use coeff
use mpetrac
use tracers
use msper01
implicit none
integer :: i,ntot,ielh,ix,iy,ic,nodno(6)
real(kind=8) :: xm,ym,xn(7),yn(7),shapef(7),rho_here,weight,weight_fixed,funccf


rho_here=1.0_8
weight_fixed = -Rb_local*dxpix*dypix
!write(irefwr,*) 'weight_fixed: ',weight_fixed
!call instop
do ix=1,nxpix
   do iy=1,nypix
      ic=(iy-1)*nxpix+ix
      !write(irefwr,*) 'ix,iy,ic: ',ix,iy,ic
      if (ipix(ic)==1) then
         ! pixel is in the dense layer
         xm=(ix-0.5)*dxpix
         ym=(iy-0.5)*dypix
         call pedetel(1,xm,ym,ielh)
         call sper01(kmeshc,coor,nodno,xn,yn,ielh)
         !write(irefwr,*) 'xm,ym: ',xm,ym
         !write(irefwr,*) 'ielh: ',ielh
         !write(irefwr,*) 'xn: ',xn(1:6)
         !write(irefwr,*) 'yn: ',yn(1:6)
         call detshape7(xn,yn,xm,ym,shapef)
         if (compress.and..not.stretch_tracers) rho_here=funccf(3,xm,ym,ym)
         weight = rho_here*weight_fixed
         !write(irefwr,*) 'weight: ',weight,shapef(1:7)
         f2_tracrhsd(1:7,ielh,i_f_tracrhsd) = f2_tracrhsd(1:7,ielh,i_f_tracrhsd) +weight*shapef(1:7)
         !write(irefwr,'(''pixrhsd: '',3i5,7e15.7)') ix,iy,i_f_tracrhsd,weight*shapef(1:7)
         !call instop
      endif
   enddo
enddo
 
end subroutine tracrhsd02_900_cart
