!   MOVE_TRACERS4
!
!   Use 4th order Runge-Kutta to move tracers. Used both in
!   predictor and corrector step; for the corrector a temporal
!   average of isolold and isolnew is used. Time step for tracing
!   is 0.5*CFL timestep.
!
!   ichoice     i    1 = predictor step; use isolold only
!                    2 = corrector step; use average of isolold and isolnew
!   kmesh,kprob i    Sepran mesh and problem definition
!   isolold     i    Velocity vector at time t
!   isolnew     i    velocity vector at time t+tstep
!   tfac        i    Fraction of CFL timestep
!   tstep       i    time step
!
!   PvK 112505
!************************************************************
subroutine move_tracers4(ichoice,isolold,isolnew,tstep,tfach)
use sepmodulekmesh
use coeff
use control
use geometry
use mparallel
use sepran_arrays
implicit none
integer :: ichoice,isolold,isolnew
real(kind=8) :: tstep,tfach

if (cyl) then
   if (periodic) then
     if (print_node) then
        write(irefwr,*) 'PERROR(move_tracers): needs to be updated for cyl and periodic'
     endif
     call instop
!    call move_tracers4_cyl_periodic(ichoice,kmesh,kprob,isolold,isolnew,user_here,tstep,tfach)
   else 
     call move_tracers4_cyl(ichoice,isolold,isolnew,tstep,tfach)
   endif
else
   call move_tracers4_cart(ichoice,isolold,isolnew,tstep,tfach)
endif
call tracvel(kmesh1,kprob1,isolnew,user_here)

end subroutine move_tracers4


