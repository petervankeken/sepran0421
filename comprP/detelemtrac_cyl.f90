subroutine detelemtrac_cyl(ic,calling_routine)
use sepmodulecomio
use sepmodulemain
use sepmodulekmesh
use control
use msper01
use mpetrac
use tracers
use mparallel
use geometry
use sepran_arrays
!use sepran_interface
implicit none
     
integer :: ic
character(len=*) :: calling_routine

real(kind=4) :: t01
integer :: j,ntot,ielh
integer :: nodno(6),nodlin(3),imissed,itrac,isub,ifound
real(kind=8) ::  xm,ym,xn(6),yn(6),un(6),vn(6),xi,eta
real(kind=8) ::  rl(3),phiq(6),u,du,func,zm,dutot,dumin,dumax
integer :: inaccurate,nmark,ichoice=0,ioffset=0
logical :: guess_first,correct,first
data first/ .true./

!if (print_node) verbose=.true.
t00 = second()

if (ic/=1.and.ic/=10) then
   if (print_node) write(irefwr,*) 'PERROR(detelemtrac_cyl): ic should be 1 or 10'
   call instop
endif

if (.not.print_node) verbose=.false.
inaccurate=0
ntot=0
if (itracoption == 1 .or. itracoption == 3 .or. itracoption == 9) then
   ntot = 0
   do idist=1,ndist
      ntot=ntot+ntrac(idist)
   enddo
else if (itracoption == 2) then
   nmark = 0
   do ichain=1,nochain
      nmark = nmark + imark(ichain)
   enddo
   ntot = nmark
endif
imissed=0
ifound=0
if ((verbose.or.petest).and.print_node) write(irefwr,*) 'ntot: ',ntot

if (ic.eq.1) then
!  *** Get element number for each tracer stored in tracer%(x,y)
   guess_first=.true.
   if (first) then 
      guess_first=.false.
      first=.false.
   endif
!  loop over tracers
   do itrac=1,ntot
      xm = tracer(itrac)%x
      ym = tracer(itrac)%y
!     try previously determined element first
      ielh = tracer(itrac)%ielem
      if (verbose) write(irefwr,'('' itrac,x,y,ielh: '',i7,2f12.3,i10,L5)') itrac,xm,ym,tracer(itrac)%ielem,guess_first
      ichoice=0
      ioffset=0
      call pedeteltrac(ichoice,ioffset,xm,ym,nodno,nodlin,rl,xn,yn,un,vn,ielh,imissed, &
     &          ifound,phiq,xi,eta,guess_first,'detelemtrac',itrac)
      tracer(itrac)%ielem=ielh
      tracer(itrac)%xi = xi ! xi_eta(1,itrac)=xi
      tracer(itrac)%eta = eta !xi_eta(2,itrac)=eta
      !if (verbose) write(irefwr,*) 'xi,eta,iel: ',xi,eta,ielh
      !write(irefwr,*) 'after first tracer determination in detelemtrac_cyl'
      !call instop
   enddo
   xi_eta_stored = .true.
   t01=second()
   if (print_node.and.petest) then
      write(irefwr,*) 'neighbor: ',needed_neighbors,needed_neighbors_neighbors,missed,xi_eta_stored
   endif

else if (ic.eq.10) then

!  test the element determination using func(10) as analytical
!  function; the nodal point values should have been stored
!  in user(10) before calling this routine
   dutot=0
   dumin=1
   dumax=0
   do itrac=1,ntot
      xm = tracer(itrac)%x
      ym = tracer(itrac)%y
      ielh = tracer(itrac)%ielem
      xi = tracer(itrac)%xi ! xi_eta(1,itrac)
      eta = tracer(itrac)%eta ! xi_eta(2,itrac)
      if (verbose) write(irefwr,*) 'before sper01: ',itrac,ielh,xi,eta
      call sper01(kmeshc,coor,nodno,xn,yn,ielh)
      call getshape6_xi_eta(xi,eta,phiq)
      do j=1,6
         un(j) = func(10,xn(j),yn(j),zm)
      enddo
      u = 0
      do j=1,6
         u = u+un(j)*phiq(j)
      enddo
      du = abs(u - func(10,xm,ym,zm))
      if (verbose) write(irefwr,*) 'du: ',du
      dumax = max(dumax,du)
      dumin = min(dumin,du)
      dutot = du+dutot
      if (du.gt.1e-5) then
         inaccurate = inaccurate+1
         if (print_node) then
           write(irefwr,*) 'PWARN(detelemtrac): ', 'interpolation is inaccurate ',itrac,u,func(10,xm,ym,zm), &
     &              (phiq(j),j=1,6),xi,eta, phiq(1)+phiq(2)+phiq(3)+phiq(4)+phiq(5)+phiq(6)
         endif
      endif
    enddo
    if (print_node) then
      write(irefwr,*) 'PINFO(detelemtrac) ', 'Av/min/max error in interpolation' ,' using stored xi,eta: '
      write(irefwr,'(15x,3e15.7,i10)') dutot/ntot,dumin,dumax, inaccurate
    endif

endif

end subroutine detelemtrac_cyl
