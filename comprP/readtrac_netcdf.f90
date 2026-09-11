!************************************************************
!   READTRAC_NETCDF
!
!   Not used for itracoption=3
!
!   PvK 970419
!************************************************************
subroutine readtrac_netcdf(nfout)
use sepmodulecomio
use mpetrac
use tracers
use coeff
use mparallel
use control
implicit none
character(len=80) :: nfname
integer :: no,nfout
!include 'SPcommon/ctimen'
include 'netcdf.inc'
integer :: ip,itrac,ntot,nmark
real(kind=8) :: x,y,global_volume
integer :: dimidds(2),len
integer :: ncid,netcdf_err,ivec_id,nrdims,varid,extra_varid
integer :: iextra_id,nextra_length,ntrac_id,ntrac_varid
integer :: NEXTRA_MAX
parameter(NEXTRA_MAX=10)
real(kind=8) :: extra_array(NEXTRA_MAX)
logical :: startexists

if (itracoption /= 1) then
   if (print_node) then
     write(irefwr,*) 'PERROR(readtrac_netcdf): not suited for'
     write(irefwr,*) '  itractoption /= 1'
   endif
   call instop
endif

if (numprocs == 1) then
    if (nfout == -1) then
       write(nfname,'(''tracers_start.nf'')') 
    else
       write(nfname,'(''Tracers/tracers.'',i4.4)') nfout
    endif
    inquire(file=nfname,exist=startexists)
    if (.not.startexists) then
       if (print_node) then
          write(irefwr,*) 'PERROR(readtrac_netcdf): tracer start'
          write(irefwr,*) ' file does not exist: ',nfname
       endif
       call instop
    endif
 else ! parallel run
    if (nfout == -1) then
       write(nfname,'(''tracers_'',i3.3,''.start.nf'')') myid
    else
       write(nfname,'(''tracers_'',i3.3,''.'',i3.3)') myid,nfout
    endif
    inquire(file=nfname,exist=startexists)
    if (.not.startexists) then
       if (print_node) then
          write(irefwr,*) 'PERROR(readtrac_netcdf): tracer start'
          write(irefwr,*) ' file does not exist: ',nfname
       endif
       call instop
    endif
    close(88)

endif ! parallel
if (print_node) write(irefwr,*) 'read tracer info from: ',nfname(1:len_trim(nfname))

!     *** open file
netcdf_err = nf_open(nfname,NF_NOWRITE,ncid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

!     *** inquire the id of vec
netcdf_err = nf_inq_dimid(ncid,'ndegfd',ivec_id)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

!     *** get the length of vec
netcdf_err = nf_inq_dimlen(ncid,ivec_id,len)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
if (print_node) write(irefwr,*) 'ndegfd=',len

!     *** inquire the id of ntrac
netcdf_err = nf_inq_dimid(ncid,'ndist',ntrac_id)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

!     *** get the length of ntrac array (ndist)
netcdf_err = nf_inq_dimlen(ncid,ntrac_id,ndist)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)
if (ndist > NDISTMAX) then
   if (print_node) then
      write(irefwr,*) 'PERROR(readtrac_netcdf): ndist is too large'
      write(irefwr,*) 'ndist = ',ndist
      write(irefwr,*) 'NDISTMAX = ',NDISTMAX
   endif
   call instop
endif
if (print_node) write(irefwr,*) 'ndist=',ndist

!     *** get varid of ntrac
netcdf_err = nf_inq_varid(ncid,'ntrac',ntrac_varid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

!     *** read ntrac data
netcdf_err = nf_get_var_int(ncid,ntrac_varid,ntrac)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

ntot = sum(ntrac(1:ndist))
if (print_node) write(irefwr,*) 'ntot: ',ntot
if (ntot /= len/2) then
   if (print_node) then
      write(irefwr,*) 'PERROR(readtracer_netcdf): inconsistent tracer '
      write(irefwr,*) ' numbers; length of tracer vector is ',len/2
      write(irefwr,*) ' but ntot computed from ntrac is     ',ntot
      write(irefwr,*) ' ndist = ',ndist
      write(irefwr,*) ' ntrac = ',(ntrac(idist),idist=1,ndist)
   endif
   call instop
endif

if (ntot>NTRACMAX) then
   if (print_node) then
     write(irefwr,*) 'PERROR(readtrac): ntot > NTRACMAX'
     write(irefwr,*) 'ntot    : ',ntot
     write(irefwr,*) 'NTRACMAX: ',NTRACMAX
     write(irefwr,*) 'ndist   : ',ndist
     write(irefwr,*) 'ntrac   : ',(ntrac(idist),idist=1,ndist)
   endif
   call instop
endif

!     *** get varid of the vector
netcdf_err = nf_inq_varid(ncid,'trac',varid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

!     *** read vector data
netcdf_err = nf_get_var_real(ncid,varid,coorreal)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

!     close file
netcdf_err = nf_close(ncid)
if (netcdf_err /= nf_noerr) call handle_err(netcdf_err)

do itrac=1,ntot
   tracer(itrac)%x=coorreal(2*itrac-1)
   tracer(itrac)%y=coorreal(2*itrac)
enddo
if (print_node) then
  write(irefwr,*) 'PINFO(readtrac_netcdf) '
  write(irefwr,*) 'ndist, ntrac : ',ndist, (ntrac(idist),idist=1,ndist)
  write(irefwr,*) 'ntot         : ',ntot
  write(irefwr,*) 'tracer(1)    : ',tracer(1)%x,tracer(1)%y
  write(irefwr,*) '               ',tracer(ntot)%x,tracer(ntot)%y
endif

 
return
end subroutine readtrac_netcdf
