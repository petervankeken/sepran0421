module tracers
  implicit none
  integer, parameter :: NDISTMAX=100
  logical :: chemwithrho
  real(kind=8) :: x0tr(NDISTMAX),x1tr(NDISTMAX),y0tr(NDISTMAX)
  real(kind=8) :: y1tr(NDISTMAX),dr,ainittr(NDISTMAX)
  real(kind=8) :: volume_dist(NDISTMAX),volume_per_tracer(NDISTMAX)
  integer :: itracoption,ntrac(NDISTMAX),ndist,idist,ifollowchem=0,ntracg(NDISTMAX)
  real(kind=8) :: radius_min,radius_max,eps_tracer_bound(NDISTMAX)

  ! cpix.inc
  integer, parameter :: NRPIXMAX=4000,NPIXMAX=4000*4000
  real(kind=8) :: rlampix,dxpix,dypix,drpix,dthpix
! integer,dimension(:),allocatable :: ipix,ielempix
! real(kind=8),dimension(:,:),allocatable :: xi_eta_pix
! real(kind=8),dimension(:),allocatable :: dth_pixel
  integer :: nxpix,nypix,nrpix,nthpix,npix_radial
  integer :: ielempix(NPIXMAX),ipix(NPIXMAX)
  real(kind=8) :: dr_pixel
  real(kind=8) :: dth_pixel(NRPIXMAX),pix(NPIXMAX),xi_eta_pix(3,NPIXMAX)

  ! c1mark.inc
  real(kind=8),dimension(10) :: wm,dm,ainit,y0m,dm_min
  integer :: nochain,ichain,imark(10)

end module tracers
