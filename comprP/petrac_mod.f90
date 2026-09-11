module mpetrac

  integer NTRACMAX,MAXELEM
  parameter(NTRACMAX= 16000000,MAXELEM=120000)
  real(kind=8),dimension(:,:,:),allocatable :: f1_tracrhsd,f2_tracrhsd
  integer :: i_f_tracrhsd
  real(kind=8),dimension(:),allocatable :: tracerheat
  real(kind=4),dimension(:),allocatable :: coorreal
  real(kind=8),dimension(:),allocatable ::  coormark,velmark,densmark,coornewm,velnewm
  integer :: ntrac_allocated=0
  integer :: tracer_number
  logical :: xi_eta_stored=.false.

  real(kind=4),allocatable,dimension(:) :: chemmark
  integer :: nchem
  logical :: duplicate_tracers,tracers_allocated,pre_allocate_tracers=.true.,predict_with_RK4=.true.

  type, public :: type_tracer
     real(kind=8) :: x,y,xnew,ynew,xhalf,yhalf,u,v,unew,vnew,xi,eta,density,CMBentertime
     integer :: ielem
     logical :: meltable,ingasable,in_cmb
  end type type_tracer
  type(type_tracer), allocatable :: tracer(:)

end module mpetrac
