subroutine detelemtrac_cart(ic,calling_routine)
use sepmodulecomio
use sepmodulemain
use sepmodulekmesh
use mparallel
use control
use mpetrac
use msper01
use tracers
use mparallel
use sepran_arrays
use sepran_interface
implicit none
     
integer :: ic
character(len=*) :: calling_routine

real(kind=8) ::  un(6),vn(6),func,du
real(kind=8) ::  xi,eta,rl(3),phiq(6),u,zm,dutot,dumin,dumax
real(kind=8) ::  xn(6),yn(6),xm,ym
integer :: ielh,itrac,nodno(6),ntot,nmark,nodlin(3)
integer :: isub
integer :: inaccurate,j
logical :: correct
   
inaccurate=0
call pefilxy(2,kmesh1,kprob1,isol1)

ntot=1
if (itracoption == 1 .or. itracoption == 3 .or. (itracoption==9.and.mpi_partrac)) then
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

!write(lu_out,*) 'in detelemtrac: ',myid,ntot

!if (print_node) pedebug=.true.
if (ic.eq.1.or.ic.eq.2) then 
  do itrac=1,ntot
     if (ic.eq.1) then
        xm = tracer(itrac)%x
        ym = tracer(itrac)%y
     else if (ic.eq.2) then
        xm = tracer(itrac)%xnew
        ym = tracer(itrac)%ynew
     endif
     call pedetel(1,xm,ym,ielh)
     if (pedebug) write(irefwr,*) 'ielh, xm,ym: ',ielh,xm,ym
     correct=.false.
     if (ielh>0) then
        call sper01(kmeshc,coor,nodno,xn,yn,ielh)
        correct=checkinelem(3,xn,yn,xm,ym,nodno,nodlin,rl, xi,eta,phiq,isub,ielh)
     endif
     if (.not.correct) then
        if (print_node) then
           write(irefwr,*) 'PERROR(velintmark): pedetel missed this one'
           write(irefwr,*) 'xm, ym, ielh: ',xm,ym,ielh,xi,eta,correct
        endif
        call instop
     endif
     tracer(itrac)%ielem=ielh
     tracer(itrac)%xi = xi ! xi_eta(1,itrac)=xi
     tracer(itrac)%eta = eta !xi_eta(2,itrac)=eta
  enddo
  xi_eta_stored =.true.

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
      !write(irefwr,'(''detelemtrac: '',4e15.7,i5,2e15.7)') xm,ym,xi,eta,ielh,u,du
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
   write(irefwr,*) 'PINFO(deteltrac): error in interpolation ', 'using stored xi,eta: '
   write(irefwr,'(15x,3e15.7)') dutot/ntot,dumin,dumax
 endif

endif

end subroutine detelemtrac_cart


