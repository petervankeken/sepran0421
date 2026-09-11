subroutine readbs_netcdf(fname,ivec)
use sepmodulecomio
use sepmodulesol
use control
use mtime
implicit none
character*(*) fname
integer :: ivec
include 'netcdf.inc'

integer :: ipvec,i,dimids(2),len
integer :: ncid,netcdf_err,ivec_id,nrdims,varid,extra_varid
integer :: iextra_id,nextra_length
integer :: nphys,nunkp,nusol,npoint=0,ndim
integer :: NEXTRA_MAX
parameter(NEXTRA_MAX=10)
real(kind=8) :: extra_array(NEXTRA_MAX),tstarth,t

call sepgetprobinfo(nphys,nunkp,nusol,ivec)
call sepgetprobinfo(nphys,nunkp,nusol,ivec)

!     *** open file
netcdf_err = nf_open(fname,NF_NOWRITE,ncid)
if (netcdf_err /= nf_noerr) then
   if (print_node) write(irefwr,*) 'file name: ',fname
   call handle_err(netcdf_err)
endif

!     *** inquire the id of length of vec
netcdf_err = nf_inq_dimid(ncid,'ndegfd',ivec_id)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
!     *** same for length of array extra
netcdf_err = nf_inq_dimid(ncid,'nextra',iextra_id)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

!     *** get the length of vec
netcdf_err = nf_inq_dimlen(ncid,ivec_id,len)
if (len /= nusol) then
   if (print_node) write(irefwr,*) 'PERROR(readbs_netcdf): length of vector'
   if (print_node) write(irefwr,*) ' on file is ',len,' which is different '
   if (print_node) write(irefwr,*) ' from the length of the defined vector: ',nusol
   call instop
endif
!     *** get the length of extra
netcdf_err = nf_inq_dimlen(ncid,iextra_id,nextra_length)
if (nextra_length > NEXTRA_MAX) then
  if (print_node) write(irefwr,*) 'PERROR(readbs_netcdf): length of extra array'
  if (print_node) write(irefwr,*) ' on file is ',nextra_length,' which is larger'
  if (print_node) write(irefwr,*) ' than the length of the defined vector: ',NEXTRA_MAX
  call instop
endif
!     *** get varid of extra array
netcdf_err = nf_inq_varid(ncid,'extra',extra_varid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
!     *** get array extra
netcdf_err = nf_get_var_double(ncid,extra_varid,extra_array)
if (restart) then
   tstarth = extra_array(1)
   time_now = tstarth
   if (print_node) write(irefwr,*) 'PINFO(readbs_netcdf): tstart_here = ',tstarth
endif
 

! get varid of the vector 
netcdf_err = nf_inq_varid(ncid,'vec',varid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

! read vector data 
netcdf_err = nf_get_var_double(ncid,varid,ks(ivec)%sol)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

! close file
netcdf_err = nf_close(ncid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

end subroutine readbs_netcdf
