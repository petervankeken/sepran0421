! ADVANCE_TRACERS_JGR97
! 
! Compute tracer advection in fixed velocity field
! to simulate benchmark from Van Keken et al., JGR, 1997
!  
!    Do for NOUT steps
!        tout = tout+dtout
!        While (t<tout)
!            compute CFL time step 
!            Integrate over new time step
!            output intermediate information
!        end
!        Output plots, write solution to file etc.
!    Done
!  
! PvK 220104
subroutine advance_tracers_JGR97()
use sepmodulecomio
use sepmodulemesh
use mtime
use geometry
use convparam
use mpetrac
use tracers
use coeff
use sepran_arrays
use control
use mparallel
use sepran_interface
#ifdef MPI
use mpi
#endif
implicit none

! LOCAL VARIABLES
real t01,t11,t10,cpu
real(kind=8) :: tnew,contln(10),format,yfaccn
integer :: jnout,iall,irhs2(5),jsmoot,numarr,ntot,myntot
integer :: inbetween,ibp,imoved,ielhlp,ipoint,ierr
real(kind=8) :: rmax,subdif,anorm,rmax2,p,q,dx,dy,func,x0,y0
real(kind=8) :: xlmin,xlmax,ylmin,ylmax,rmin,x,y,dxt,dyt
logical :: output
character(len=80) :: fimage,command

save irhs2

niter = 0
time_now=0
cpu=0
jnout = 1
inbetween = 0
imoved = 0
dxt=0
dyt=0
! ???
!if (cyl) then
!   dxt=0.4192
!   dyt=-0.5
!endif

!if (petest) call instop

if (print_node) t00 = second()
toutp = time_now + dtoutp
output = .false.

if (print_node.and..not.mpi_partrac) open(9,file='JGR97_path.dat')
! calculate time step; pefilcof stores velocity in user
! which then is used in pedtcf.
if (print_node) write(irefwr,*) 'call pefilcof2'
call pefilcof(2,iuser_here,user_here)
call pedtcf(kmesh1,iuser_here,user_here,dtcfl)
t00 = second()
if (print_node) write(irefwr,'(''   Initial condition, dtcfl: '',3f12.7)') tracer(1)%x,tracer(1)%y,dtcfl
x0 = tracer(1)%x-dxt
y0 = tracer(1)%y-dyt
if (print_node.and..not.mpi_partrac) write(9,*) x0,y0
xlmax = x0
xlmin = x0
ylmax = y0
ylmin = y0
rmin=r2
rmax=r1
!write(6,*) 'start loop',myid,ntrac(1)
do
   niter = niter+1
   
!  based time step on CFL criterion
   tstepp = dtcfl * tfac
     
!  make sure new time doesn't exceed output limits or maximum time
   tnew = time_now+tstepp
   if (tnew.ge.tmaxp) then
      tstepp = tmaxp-time_now
      time_now = tmaxp
      output=.true.
   else if (tnew.ge.toutp) then
      tstepp = toutp-time_now
      time_now = toutp
      output = .true.
   else
      time_now = tnew
      output = .false.
   endif
   !if (print_node) write(irefwr,*) niter,tnew

!  predict position of markers
   if (itracoption /= 0) then
!     interpolate velocity in old positions of markers
      call tracvel(kmesh1,kprob1,isol1,user_here)
      call move_tracers4(1,isol1,isol1,tstepp,tfac)
   endif
   ! no need to do corrector
   call copcoor
   xlmax = max ( xlmax , tracer(1)%x -dxt)
   xlmin = min ( xlmin , tracer(1)%x -dxt)
   ylmax = max ( ylmax , tracer(1)%y -dyt)
   ylmin = min ( ylmin , tracer(1)%y -dyt)
   x=tracer(1)%x-dxt
   y=tracer(1)%y-dyt
   rmin=min(rmin,sqrt(x*x+y*y))
   rmax=max(rmax,sqrt(x*x+y*y))
  
   if (print_node.and..not.mpi_partrac) write(9,*) tracer(1)%x-dxt,tracer(1)%y-dyt

   if (time_now>=tmaxp) exit
enddo
!write(6,*) 'done with loop: ',myid

if (print_node.and..not.mpi_partrac) then
  t01 = second()
  write(irefwr,'(''    Final time: '',f12.7)') time_now
  write(irefwr,'(''    Final x,y : '',2f12.7)') tracer(1)%x-dxt,tracer(1)%y-dyt
  dx = (tracer(1)%x-x0-dxt)
  dy = (tracer(1)%y-y0-dyt)
  write(irefwr,'(''    Error     : '',3f12.7)') dx,dy,sqrt(dx*dx+dy*dy)
  write(irefwr,'(''    CPU+niter: '',f12.2,2i10)') t01-t00,niter,ntrac(1)
  write(irefwr,'(''    Tracer path stored in JGR97_path.dat'')')
  close(9)
  write(irefwr,*) 'Path x min/max: ',xlmin,xlmax
  write(irefwr,*) 'Path y min/max: ',ylmin,ylmax
  write(irefwr,*) 'Path r min/max: ',rmin,rmax ! ,r1-rmin,r2-rmax
endif
if (mpi_partrac) then
   ntot=ntrac(1)
   myntot=ntrac(1)
!  write(6,*) 'before MPI_REDUCE: ',myid
#ifdef MPI
!  call MPI_REDUCE(myntot,ntot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#endif
  if (print_node) then
     t01=second()
     write(irefwr,'(''    CPU+niter: '',f12.2,2i10)') t01-t00,niter,ntot*numprocs
  endif
endif
call tracerout(0)

end subroutine advance_tracers_JGR97
