!   TRACEROUT
!
!   Not used for itracoption=3
!
!   PvK 970413
subroutine tracerout(no)
use sepmodulecomio
use mtime
use tracers
use mpetrac
use mparallel
use control
implicit none
integer :: no
character(len=80) :: tname
integer :: ip,i,ntot
real(kind=8) :: x,y

if (netcdf) then
   call tracerout_netcdf(no)
   return
endif


if (itracoption == 1 .or. (itracoption==9.and.mpi_partrac)) then
  if (.not.mpi_partrac) then
    write(tname,'(''Tracers/tracers.'',i4.4)') no
    open(88,file=tname,form='unformatted')
    write(tname,'(''Tracers/tracers_header.'',i4.4)') no
    open(89,file=tname,form='formatted')
 else 
    write(tname,'(''Tracers/tracers_'',i4.4,''.'',i4.4)') myid,no
    open(88,file=tname,form='unformatted')
    write(tname,'(''Tracers/tracers_header_'',i4.4,''.'',i4.4)') myid,no
    open(89,file=tname,form='formatted')
 endif
 write(89,*) ndist,(ntrac(idist),idist=1,ndist)
 close(89)

  ntot=0
  do idist=1,ndist
     ntot=ntot+ntrac(idist)
  enddo


  ip=0
  do idist=1,ndist
    write(88) (real(tracer(ip+i)%x),real(tracer(ip+i)%y),i=1,ntrac(idist))
    !write(88) (coorreal(1,ip+i),coorreal(2,ip+i),i=1,ntrac(idist))
    ip = ip+ntrac(idist)
  enddo
  close(88)
  if (print_node) write(irefwr,*) 'Tracers stored at t = ',time_now,no
endif
 
if (itracoption == 2) then
  if (numprocs>1) then
!          *** this is likely an error: no need to use parallel with markerchain
     if (print_node) then
       write(irefwr,*) 'PERROR(tracerout): numprocs > 1 for' 
       write(irefwr,*) 'itracoption=2: ',numprocs
     endif
     call instop()
  endif
  write(tname,'(''solutions/markers.'',i4.4)') no
  open(88,file=tname)
  do i=1,imark(1)
     !write(88,*) coormark(1,i),coormark(2,i)
     write(88,*) tracer(i)%x,tracer(i)%y
  enddo
  close(88)
  if (print_node) write(irefwr,*) 'Markers stored at t = ',time_now,no
endif
  
end subroutine tracerout
