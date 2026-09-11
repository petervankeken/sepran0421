subroutine writbs_netcdf(fname,ivec)
use sepmodulesol
use sepmodulecomio
use control
use convparam
use mtime
implicit none
character(len=*) :: fname
integer :: ivec
include 'netcdf.inc'

integer :: i,dimids(2),length_file_name
integer :: ncid,netcdf_err,ivec_id,nrdims,varid,extra_varid
integer :: iextra_id,nextra_length,nphys,nunkp,nusol,ndim,npoint
integer :: NEXTRA_MAX
parameter(NEXTRA_MAX=10)
real(kind=8) :: extra_array(NEXTRA_MAX),t
character(len=80) :: fname2,fname_in


call sepgetmeshinfo(ndim,npoint)
call sepgetprobinfo(nphys,nunkp,nusol,ivec)
!if (ivec==1.and.nusol/=ndim*npoint.and.print_node) then
!   write(irefwr,*) 'nous avons un petit problem: ',nusol,ivec,npoint
!   call instop
!endif

fname_in=fname(1:80)
length_file_name=len_trim(fname_in)
fname2=fname_in(1:length_file_name)

!write(irefwr,*) 'npoint = ',npoint
!write(irefwr,*) 'nusol  = ',nusol
!write(irefwr,*) 'writbs_netcdf: ',ivec,nusol,npoint

! create netcdf file
netcdf_err = nf_create(fname2,NF_CLOBBER,ncid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

! Metadata
! length of data in vector
netcdf_err = nf_def_dim(ncid,'ndegfd',nusol,ivec_id)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
! length of extra parameters (time etc.
nextra_length = 1
if (nextra_length > NEXTRA_MAX) then
   if (print_node) then
      write(irefwr,*) 'PERROR(netcdf_writbs): nextra_length is too'
      write(irefwr,*) 'large : ',nextra_length,NEXTRA_MAX
   endif
   call instop 
endif
netcdf_err = nf_def_dim(ncid,'nextra',nextra_length,iextra_id)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
nrdims=1
netcdf_err = nf_def_var(ncid,'vec',NF_DOUBLE, nrdims,ivec_id,varid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
nrdims=1
netcdf_err = nf_def_var(ncid,'extra',NF_DOUBLE, nrdims,iextra_id,extra_varid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
! End of definition mode

netcdf_err = nf_enddef(ncid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
    
! write data to file
! first write extra info (time etc.)
extra_array(1)=time_now
netcdf_err = nf_put_var_double(ncid,extra_varid,extra_array)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
! then write vector
netcdf_err = nf_put_var_double(ncid,varid,ks(ivec)%sol)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

! close file
netcdf_err = nf_close(ncid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

end subroutine writbs_netcdf


