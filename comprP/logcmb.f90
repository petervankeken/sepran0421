! log when particles get into or go out of the CMB
! PvK October 2019
subroutine logcmb(ichoice)
use sepmodulecomio
use coeff
use tracers
use mpetrac
use mparallel
use geometry
use brandenburg
use control
use convparam
use mtime
implicit none
integer :: ichoice
logical :: extract
integer :: itrac,ip1,ip2,ntot,idum
integer :: imelt,izone,lu,iingas,entered_cmb,left_cmb
real(kind=4) :: ran0
real(kind=8) :: basmin,basmax,harzmin,harzmax,ybefore
real(kind=8) :: r,theta,cost,rbefore,x,y,dtheta

if (.not.cyl) then
   write(irefwr,*) 'PERROR(logcmb): only suited for cyl right now'
   call instop
endif

if (mpi_parallel) then
   write(irefwr,*) 'PERROR(logcmb): logging of particles needs to be extended to parallel'
   call instop
endif


ntot=sum(ntrac(1:ndist))
if (ichoice==2) then
   ! final call. Flush remaining tracers in cmb
   do itrac=1,ntot
      if (tracer(itrac)%in_cmb) then
         write(LU_CMB_ENTRY,*) itrac,tracer(itrac)%CMBentertime,niter
      endif
   enddo
   call flush(LU_CMB_ENTRY)
   close(LU_CMB_ENTRY)
   return
endif

if (ichoice==0) then
   ! first call. Log all tracers in cmb
   entered_cmb=0
   do itrac=1,ntot
      x=tracer(itrac)%x
      y=tracer(itrac)%y
      r=sqrt(x*x+y*y)
      if (r<=r1+d_cmb) then 
         tracer(itrac)%in_cmb=.true.
         write(LU_CMB_ENTRY,*) itrac,0  
         tracer(itrac)%CMBentertime=0
         entered_cmb=entered_cmb+1
      endif  
   enddo
   flush(LU_CMB_ENTRY)
   write(irefwr,*) 'on first iteration ',entered_cmb,' tracers are in CMB'
   return
endif
 

entered_cmb=0
left_cmb=0
do itrac=1,ntot
   ! loop over all tracers
   x=tracer(itrac)%x
   y=tracer(itrac)%y
   r=sqrt(x*x+y*y)
   if (.not.tracer(itrac)%in_cmb.and.r<=r1+d_cmb) then
      ! log entry into the CMB
      tracer(itrac)%in_cmb=.true. ! note that tracer got to the CMB
      entered_cmb=entered_cmb+1
      tracer(itrac)%CMBentertime=niter
      write(LU_CMB_ENTRY,*) itrac,niter
   endif
   if (tracer(itrac)%in_cmb.and.r>r1+d_cmb) then
      ! log departure from CMB
      write(LU_CMB_EXIT,*) itrac,tracer(itrac)%CMBentertime,niter
      tracer(itrac)%CMBentertime=-1
      tracer(itrac)%in_cmb=.false. ! not that tracer left the CMB
      left_cmb=left_cmb+1
   endif
enddo
write(irefwr,'(''at time_now '',e14.7,'' tracers entering and leaving CMB: '',2i10)') time_now/tscale_dim,entered_cmb,left_cmb
call flush(LU_CMB_ENTRY)
call flush(LU_CMB_EXIT)

end subroutine logcmb
