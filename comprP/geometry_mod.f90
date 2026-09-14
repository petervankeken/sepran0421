module geometry
real(kind=8),parameter :: pi = 3.14159265358979323846_8, twopi=2.0_8*pi
real(kind=8) :: r1,r2,frac,rmid
real(kind=8) :: T_top,T_bot
real(kind=8) :: volume,surf(10),wavel_perturb,surfbot,surftop
integer :: ibottom,itop,imid,nbp,iboundpoints(10),icloc(10),iclc,ntheta
logical :: cyl,axi,half,quart,periodic,full,insulbot,eighth
integer :: isideboundary

! pexcyc.inc
integer, parameter :: NXCYCMAX=10000
real(kind=8),dimension(NXCYCMAX) :: xc,yc
real(kind=8) ::xcmin,xcmax,ycmin,ycmax,dx_two_elements
integer :: nx,ny

integer :: ifilchoice,ncurvn(4),icurv(4,10)
real(kind=8) :: dth_output,dr_output
integer :: nth_output,nr_output
logical :: rav_output

integer :: nlay
real(kind=8) :: r_zint(10),zint_d(10)

integer,parameter :: NELEM_TOPO=120000
real(kind=8) :: rtop,rbot
real(kind=8) :: rtop_max,rtop_min,rtop_threshold,eps_elem
real(kind=8) :: rbot_max,rbot_min,rbot_threshold
integer, dimension(NELEM_TOPO) :: ipxmin,ipxmax,ipymin,ipymax
integer, dimension(NELEM_TOPO) ::  ipix_rmin,ipix_rmax,ipix_thmin,ipix_thmax
logical, dimension(NELEM_TOPO) ::  boundary_elem,curved_elem


integer,parameter :: NRDIVMAX=1001,NDIVMAX=10000
real(kind=8) :: drgrid,dthgrid(NRDIVMAX),rgrid(NRDIVMAX)
integer :: nrdiv,nthgrid(NRDIVMAX),ithgrid(NRDIVMAX),ngridtot,idiv(NDIVMAX)

integer :: ncontzone,imovecont
real(kind=8) :: thmincont(10),thmaxcont(10),revol_cont,rcont

integer,parameter :: NXELM=3001,NYELM=3001
real(kind=8) :: dxel,dyel,xoff,yoff,grid_height,grid_width
integer :: nxel,nyel,ielgrid(NXELM,NYELM)


end module geometry
