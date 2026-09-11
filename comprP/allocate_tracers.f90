subroutine allocate_tracers(ntrac,nchemmax,fname)
use brandenburg
use mpetrac
use control
use coeff
implicit none
integer :: ntrac,nchemmax
character(len=*) :: fname
integer :: alloc_status

allocate(tracer(ntrac),stat=alloc_status); call trap_allocate(alloc_status,fname,'tracer')
allocate(coorreal(2*ntrac),stat=alloc_status); call trap_allocate(alloc_status,fname,'coorreal')
if (nchemmax>0) then
  allocate(chemmark(nchemmax*ntrac),stat=alloc_status); call trap_allocate(alloc_status,fname,'chemmark')
endif

tracers_allocated=.true.
ntrac_allocated=ntrac

end subroutine allocate_tracers

subroutine deallocate_tracers()
use coeff
use mpetrac
use brandenburg
implicit none

if (allocated(tracer)) deallocate(tracer)
if (allocated(tracerheat)) deallocate(tracerheat)
if (allocated(coorreal)) deallocate(coorreal)
if (allocated(chemmark)) deallocate(chemmark)


tracers_allocated=.false.

end subroutine deallocate_tracers
