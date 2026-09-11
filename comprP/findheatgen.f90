!   FINDHEATGEN
!
!   Determine radiogenic heatproduction in nodal points
!
!   kmesh,kprob  i    Mesh and problem definition
!   coormark     i    positions of the tracers
!   tracerheat   i    Array containing heating in the tracer
!   idens        o    vector of special structure containing 
!                     density of the tracers in the nodal points
!                     (only if iqtype=3)
!   iheat        o    vector of special structure containing
!                     heatprodcction in the nodal points
!
!   PvK 080600
subroutine findheatgen(idens,iheat)
use sepmodulecomio
use sepmodulekmesh
use sepmodulekprob
use mtime
use convparam
use coeff
use control
implicit none
integer ::  idens,iheat

if (compress.and.qwithrho) then
   if (print_node) write(irefwr,*) 'PERROR(findheatgen): not yet suited for compress and qwithrho'
   call instop
endif

!if (iqtype/=0) then
!   write(irefwr,*) 'PERROR(findheatgen): not yet suited for iqtype/=0'
!   call instop
!endif

call findheatgen01(npoint,ks(iheat)%sol)


end subroutine findheatgen

subroutine findheatgen01(npoint,internal_heating)
use sepmodulecomio
use control
use coeff
implicit none
integer :: npoint
real(kind=8) :: internal_heating(*)

if (iqtype/=0 .and. iqtype/=1) then
   if (print_node) then
     write(irefwr,*) 'PERROR(findheatgen): not yet suited for anything but constant heating'
     write(irefwr,*) 'iqtype / qwithrho = ',iqtype,qwithrho
   endif
   call instop
endif

if (iqtype==0) then
   internal_heating(1:npoint)=0.0_8
else if (iqtype==1) then
   internal_heating(1:npoint)=q_layer(1)
endif

end subroutine findheatgen01
