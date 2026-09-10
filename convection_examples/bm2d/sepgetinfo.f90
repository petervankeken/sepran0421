! PvK July 2024
subroutine sepgetprobinfo(nphys_here,nunkp_here,nusol_here,isol_here)
use sepmodulemain
use sepmodulekprob
implicit none
integer :: nphys_here, nunkp_here, isol_here, nusol_here

call sepactsolbf1(isol_here)
nphys_here = nphysi
nunkp_here = nunkpi
nusol_here = nusoli
end subroutine sepgetprobinfo

subroutine sepgetmeshinfo(ndim_here,npoint_here)
use sepmodulemain
use sepmodulekmesh
use sepmodulecpack
implicit none
integer :: ndim_here, npoint_here

!write(irefwr,*) 'ndim = ',ndim,npoint,inodenr
ndim_here = ndim
npoint_here= npoint
end subroutine sepgetmeshinfo


subroutine sepgetrhsdinfo(irhsd)
use sepmodulecomio
use sepmodulevecs
use control
implicit none
integer :: irhsd
if (print_node) write(irefwr,*) 'rhsd complexity: ',ks(irhsd)%compl

end subroutine sepgetrhsdinfo

