!   READTRAC
!
!   Not used for itracoption=3
!
!   PvK 970419
subroutine readtrac(jnout)
use control
implicit none
integer :: jnout

if (netcdf) then
   call readtrac_netcdf(jnout)
else
   call readtrac_old(jnout)
   return
endif

end subroutine readtrac

subroutine readtrac_old(jnout)
use sepmodulecomio
use mtime
use mpetrac
use tracers
use coeff
use mparallel
use control
implicit none
character(len=80) :: tname
integer :: no
!include 'SPcommon/ctimen'
integer :: ip,itrac,ntot,nmark,jnout
real(kind=8) :: x,y,global_volume

if (netcdf) then
   call readtrac_netcdf(jnout)
   return
endif
ip=0
if (itracoption == 1) then
   if (numprocs == 1) then
      if (jnout == -1) then
         write(tname,'(''tracers_header.start'')') 
      else
         write(tname,'(''tracers_header.'',i3.3)') jnout
      endif
      open(88,file=tname)
      read(88,*) ndist,(ntrac(idist),idist=1,ndist)
      ntot = sum(ntrac(1:ndist))
      write(irefwr,*) 'PINFO(readtrac): ndist = ',ndist
      write(irefwr,*) '                 ntrac = ', (ntrac(idist),idist=1,ndist)
      if (ntot>NTRACMAX) then
         write(irefwr,*) 'PERROR(readtrac): ntot > NTRACMAX'
         write(irefwr,*) 'ntot    : ',ntot
         write(irefwr,*) 'NTRACMAX: ',NTRACMAX
         write(irefwr,*) 'ndist   : ',ndist
         write(irefwr,*) 'ntrac   : ',(ntrac(idist),idist=1,ndist)
         call instop
      endif
      if (jnout == -1) then
         write(tname,'(''tracers.start'')') 
      else
         write(tname,'(''tracers.start'',i3.3)') jnout
      endif
      open(89,file=tname,form='unformatted') 
      do idist=1,ndist
         write(irefwr,*) 'read : ',ip,2*(ip+1)-1,2*(ntrac(idist)+1)-1
        read(89) (coorreal(2*(ip+itrac)-1),coorreal(2*(ip+itrac)),itrac=1,ntrac(idist))
          ip = ip+ntrac(idist)
        ! this should go into the tracer file
        volume_dist(idist) = (x1tr(idist)-x0tr(idist))*(y1tr(idist)-y0tr(idist))
        volume_per_tracer(idist)=volume_dist(idist)/ntrac(idist)
        if (volume_dist(idist) < 1e-9) then
           write(irefwr,*) 'PERROR(readtrac_netcdf): volume of ', ' distribution is very small: ',volume_dist(idist)
           write(irefwr,*) 'distribution : ',idist
           call instop
        endif
      enddo
      write(irefwr,*) 'done reading'
      close(88)
      close(89)
   else ! parallel run
      if (jnout == -1) then
         write(tname,'(''tracers_header_'',i3.3,''.start'')') myid
      else
         write(tname,'(''tracers_header_'',i3.3,''.'',i3.3)') myid,jnout
      endif
 
      open(88,file=tname)
      read(88,*) ndist,(ntrac(idist),idist=1,ndist)
      ntot = sum(ntrac(1:ndist))
      write(irefwr,*) 'PINFO(readtrac): myid, ndist,ntrac = ', myid,ndist,(ntrac(idist),idist=1,ndist)
      if (ntot>NTRACMAX) then
         write(irefwr,*) 'PERROR(readtrac): ntot > NTRACMAX'
         write(irefwr,*) 'ntot    : ',ntot
         write(irefwr,*) 'NTRACMAX: ',NTRACMAX
         write(irefwr,*) 'ndist   : ',ndist
         write(irefwr,*) 'ntrac   : ',(ntrac(idist),idist=1,ndist)
         call instop
      endif
      if (jnout == -1) then
         write(tname,'(''tracers_'',i3.3,''.start'')') myid
      else
         write(tname,'(''tracers_'',i3.3,''.'',i3.3)') myid,jnout
      endif
      open(89,file=tname,form='unformatted') 
      do idist=1,ndist
        read(89) (coorreal(2*(ip+itrac)-1),coorreal(2*(ip+itrac)), itrac=1,ntrac(idist))
          ip = ip+ntrac(idist)
      enddo
      close(88)
      close(89)
   endif


   write(irefwr,*) 'ntot = ',ntot
   do itrac=1,ntot
     tracer(itrac)%x = coorreal(2*itrac-1)
     tracer(itrac)%y = coorreal(2*itrac)
   enddo

endif  ! itracoption == 1
 
if (itracoption == 2) then
   open(88,file=tname)
   itrac=1
100      continue
     read(88,*,end=200) tracer(itrac)%x,tracer(itrac)%y
     itrac=itrac+1
     if (itrac.lt.NTRACMAX) goto 100
200      nmark=itrac
   imark(1)=nmark
endif
   

end subroutine readtrac_old
