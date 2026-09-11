module brandenburg
! cjp.inc
real(kind=8) :: poffset,pfactor,x_nought=2.0_8,pb_param(3)
real(kind=8) :: b,c,chd,chs,clayer,chs_P,q_beta
real(kind=8) :: odrelax,trestart
real(kind=8) :: rcmbia(2)
integer :: iviscositylaw,isurfbound,ismove
integer :: ifreq,nbasdist(6)
integer :: gable_stokes_choice=0,iCHBT
integer :: modet=0,restartniter,restartout
logical :: Gzero,iweak,track_cmb_ingas=.false.
logical :: FBdebug,FBvisc,extractbas,duplicate,initialize_crust
logical :: restartchem,move_Hz_from_crust,old_gable_plates
integer, parameter :: NPLATEMAX=10
real(kind=8),dimension(NPLATEMAX) :: plate_vel
integer,dimension(NPLATEMAX) :: meltzone,ingaszone,plate_rhsd=0
logical,dimension(NPLATEMAX) :: diverging
integer :: nmeltzones,ningaszone
real(kind=8),allocatable :: GP(:,:),TP(:)

! strchain.inc
! Declare the max length of the marker chain
integer,parameter :: NCHAINMAX=20000

! real coordinates of the marker chain...this is set at t=0 and does not change
real(kind=8) :: xstrchain(NCHAINMAX),ystrchain(NCHAINMAX),thetachain(NCHAINMAX)
real(kind=8) :: xbpoint(20),ybpoint(20)

! Values of dphi/dx and dphi/dy in the actual element
real(kind=8) :: dphidx(6,NCHAINMAX),dphidy(6,NCHAINMAX)

! Topology: Link each tracer to its 6 corresponding nodes...this is set at t=0 and does not change
integer :: itopo_chain(6,NCHAINMAX)
real(kind=8) :: xnchain(6,NCHAINMAX),ynchain(6,NCHAINMAX),shapefchain(6,NCHAINMAX)

! Strain Rate Tensor... only need to store 3 components b/c of symmetry; this will change
! strchain_tensor(1,*) = 2(du/dx)
! strchain_tensor(2,*) = (du/dy + dv/dx)
! strchain_tensor(3,*) = 2(dv/dy)
real(kind=8) :: stten(3,NCHAINMAX)

! Shear Stress... this vector holds the shear stress computed from (tau = eta*gammadot)
real(kind=8) :: sttau(NCHAINMAX),bttau(20)

! Other marker chain information
real(kind=8) :: scy,scxmin,scxmax,strc_r
real(kind=8) :: plate_boundaries(NPLATEMAX+1),whalf
real(kind=8) :: zmelt,xmelt,xmelt_left(NPLATEMAX+1),xmelt_right(NPLATEMAX+1),zcrust,zharz,d_cmb
! imeltnow(1,i): number of basalt tracers melted in melt zone i
! imeltnow(2,i): number of harzburgite tracers melted in melt zone i (when using ! ratio method)
integer :: nsmark,intlink(2,NPLATEMAX+1),nplates,istress_curve,nmeltzone,imeltnow(2,NPLATEMAX+1)
integer :: plate_sol(NPLATEMAX)=0,irotation=0 ! solution arrays pointing towards each plate solution
logical :: gable_plates=.false.,recomputeG=.false.

end module brandenburg
