! PvK July 2024
subroutine intermediate_output()
use sepmodulekmesh
use sepmoduleoldrouts
use sepran_arrays
use geometry
use coeff
use control
implicit none
real(kind=8),allocatable,dimension(:) :: funcx,funcy
real(kind=8) :: smax,tmax,vrms,gnus,temp1,q1,q2
integer :: ielhlp=0,ihelp=0,irule,icurvs(3),nunkp,allocate_status

t2=second()
dcpu=t2-t1
cput=t2-cpustart
t1=t2

! check difference between temperature solution from last iteration
dif=anorm(0,3,0,kmesh,kprob,isol(2),islold(2),ielhlp)
tmax=anorm(1,3,0,kmesh,kprob,isol(2),isol(2),ielhlp)
dif=dif/tmax

if (.not.allocated(funcx)) then
   allocate(funcx(6+2*npoint),stat=allocate_status)
   if (allocate_status/=0) then
      if (print_node) write(irefwr,*) 'PERROR(intermediate_output): allocate on funcx failed: ',allocate_status   
      call instop
   endif
endif
if (.not.allocated(funcy)) then
   allocate(funcy(6+npoint),stat=allocate_status)
   if (allocate_status/=0) then
      if (print_node) write(irefwr,*) 'PERROR(intermediate_output): allocate on funcy failed: ',allocate_status   
      call instop
   endif
endif
funcx(1)=1.0_8*(6+2*npoint)
funcy(1)=1.0_8*(6+npoint)

! Vrms, Nu, and heatflow at top corners
call pevrms(vrms)
irule=1
call deriva(0,2,0,1,0,igradt,kmesh,kprob,isol(2),isol(2),iuser_here,user_here,ihelp)
temp1 = bounin(1,irule,1,1,kmesh,kprob,1,1,isol(2),iuser_here,user_here)
if (temp1 < 1e-7) then
   if (print_node) write(irefwr,*) 'PERROR(intermediate_output): temp along curve 1 is close to 0: ',temp1
   call instop
endif
! This should become part of a general nusselt() subroutine
gnus = -bounin(1,irule,1,2,kmesh,kprob,3,3,igradt,iuser_here,user_here)/temp1
icurvs(1)=-1
icurvs(2)=3
icurvs(3)=3
call compcr(0,kmesh,kprob,igradt,2,icurvs,funcx,funcy)
nunkp=nint(funcy(5))
q2 = -funcy(6)
q1 = -funcy(5+nunkp)


write(6,'(i5,f12.4,4f12.3,2f8.2)') niter,dif,vrms,gnus,q1,q2,dcpu,cput

if (allocated(funcx)) deallocate(funcx)
if (allocated(funcy)) deallocate(funcy)

end subroutine intermediate_output
