!   INITTRACHEAT
!
!   Create initial distribution of tracer heating
!
!   062710 PvK: simple version to set up tracer-dependent heating
subroutine inittracheat(tracerheat,kmesh,kprob,isol)
use sepmodulecomio
use convparam
use control
use tracers
implicit none
integer :: kmesh(*),kprob(*),isol(*)
real(kind=8) :: tracerheat(*)
integer :: ntot
real(kind=8) :: x,y,r
integer :: ibuffr
common ibuffr(1)
real(kind=8) :: buffr(1)
equivalence(buffr(1),ibuffr(1))

if (print_node) write(irefwr,*) 'PERROR(inittracheat): subroutine is obsolete'
call instop

if (itracoption == 2) then
   if (print_node) then
      write(irefwr,*) 'PERROR(initchem): not suited for markerchain'
      write(irefwr,*) 'method'
   endif
   call instop
endif

ntot = 0
do idist=1,ndist
   ntot=ntot+ntrac(idist)
enddo

if (print_node) write(irefwr,*) 'update for inittracheat'
call instop

! Create array of length ntot that will contain the tracer heating
! call inittracheat01(tracerheat,ntot,coortrac)

end subroutine inittracheat

!!! subroutine inittracheat01(tracerheat,ntot,coortrac)
!!! implicit none
!!! integer :: ntot
!!! real(kind=8) :: tracerheat(ntot),coortrac(2,ntot)
!!! integer :: i
!!! real(kind=8) :: func_heat
!!! 
!!! do i=1,ntot
!!! !  for test comparison with other heating specifications
!!!       tracerheat(i) = func_heat(1,coortrac(1,i),coortrac(2,i))
!!! enddo
!!!   
!!! end subroutine inittracheat01
