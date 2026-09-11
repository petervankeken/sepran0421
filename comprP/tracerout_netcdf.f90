subroutine tracerout_netcdf(no)
use sepmodulecomio
use mtime
use tracers
use mparallel
use mpetrac
use control
implicit none
integer :: no
include 'netcdf.inc'
!include 'SPcommon/ctimen'
character(len=80) :: fname
integer :: ip,itrac,ntot

integer :: dimids(2)
integer :: ncid,netcdf_err,ivec_id,nrdims,varid,extra_varid
integer :: iextra_id,nextra_length,ntrac_id,ntrac_varid
integer :: NEXTRA_MAX
parameter(NEXTRA_MAX=10)
real(kind=8) :: extra_array(NEXTRA_MAX)

if (itracoption == 2.and.mpi_partrac) then
   if (print_node) then
      write(irefwr,*) 'PERROR(tracerout_netcdf): not yet available'
      write(irefwr,*) 'for itracoption == 2 and mpi_partrac ',itracoption,mpi_partrac
   endif
   call instop
endif



if (itracoption==1 .or. (itracoption==9 .and. mpi_partrac)) then
   ntot=sum(ntrac(1:ndist))
   if (.not.mpi_partrac) then
       write(fname,'(''Tracers/tracers.'',i4.4)') no
   else
       write(fname,'(''Tracers/'',i3.3,''/tracers.'',i4.4)') myid,no
   endif
else if (itracoption==2) then
   ntot=sum(imark(1:nochain))
   write(fname,'(''Tracers/markers.'',i4.4)') no
endif

do itrac=1,ntot
   coorreal(2*itrac-1)=real(tracer(itrac)%x)
   coorreal(2*itrac)=real(tracer(itrac)%y)
enddo

! create netcdf file
netcdf_err = nf_create(fname,NF_CLOBBER,ncid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

! Metadata
! length of data in vector
netcdf_err = nf_def_dim(ncid,'ndegfd',2*ntot,ivec_id)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
! length of ntrac array (ndist)
netcdf_err = nf_def_dim(ncid,'ndist',ndist,ntrac_id)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
! length of extra parameters (time etc.
nextra_length = 1
if (nextra_length > NEXTRA_MAX) then
   if (print_node) write(irefwr,*) 'PERROR(netcdf_writbs): nextra_length is too'
   if (print_node) write(irefwr,*) 'large : ',nextra_length,NEXTRA_MAX
   call instop 
endif
netcdf_err = nf_def_dim(ncid,'nextra',nextra_length,iextra_id)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
nrdims=1
netcdf_err = nf_def_var(ncid,'trac',NF_FLOAT, nrdims,ivec_id,varid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
nrdims=1
netcdf_err = nf_def_var(ncid,'extra',NF_DOUBLE, nrdims,iextra_id,extra_varid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
netcdf_err = nf_def_var(ncid,'ntrac',NF_INT, nrdims,ntrac_id,ntrac_varid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
! End of definition mode

netcdf_err = nf_enddef(ncid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
    
! write data to file
! *** first write extra info (time etc.)
extra_array(1)=time_now
netcdf_err = nf_put_var_double(ncid,extra_varid,extra_array)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
! then write tracer numbers
netcdf_err = nf_put_var_int(ncid,ntrac_varid,ntrac)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
! then write tracers
netcdf_err = nf_put_var_real(ncid,varid,coorreal)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)


! close file
netcdf_err = nf_close(ncid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

end subroutine tracerout_netcdf
