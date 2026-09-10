! create velocity and temperature solution vectors
! PvK July 2024
subroutine initialize_temp_vel()
use sepran_arrays
use control
implicit none
integer :: iprob_here,idum=0,ichvc,iu1lc(2)=0
real(kind=8) :: u1lc(2)=0.0_8

! create vectors, set up with conductive + perturbation for temperature and zero velocity
!  call create(0,kmesh,kprob,isol(1))
!  call create(0,kmesh,kprob,isol(2))
write(6,*) 'create temperature vector'
ichvc=1
iprob_here=2
iu1lc(1)=1 ! set to func(1,x,y,z)
call creavc(0,ichvc+(iprob_here-1)*1000,idum,isol(2),kmesh,kprob,iu1lc,u1lc)
write(6,*) 'create velocity vector'
iprob_here=1
iu1lc(1)=0 ! set u to zero
iu1lc(2)=0 ! set v to zero
call creavc(0,ichvc+(iprob_here-1)*1000,idum,isol(1),kmesh,kprob,iu1lc,u1lc)

! overwrite with stored solutions for irestart>=0 
if (irestart == 0) then
   call readbs_netcdf(Tstartfile,isol(2))
   if (read_velocity) call readbs_netcdf(UVstartfile,isol(1))
endif
! initialize islold(*)
call copyvc(isol(1),islold(1))
call copyvc(isol(2),islold(2))

end subroutine initialize_temp_vel
