subroutine remarker
use geometry
implicit none
if (cyl) then
   call remarker_cyl
else
   call remarker_cart
endif
return
end

!   REMARKER_CYL
!
!   Regrid the markerchain. Check if the interval between
!   markers has become larger than the tolerance DMAX. If so,
!   add a marker in the center of this interval using linear
!   interpolation.
!
!   PvK 960830
!
!   PvK 072604
subroutine remarker_cyl
use sepmodulecomio
use tracers
use mpetrac
use coeff
use geometry
use control
implicit none
integer :: nmnew(10)
integer :: ipm,ipn,nmark,inew,iold,ip,ip1,inewtot
real(kind=8) :: dx,dy,rl,x1,x2,y1,y2,distance,xnew,ynew,rnew,thetanew
real(kind=8) :: rm1,rm2,theta1,theta2

if (.not.quart.and..not.half.and..not.eighth) then
! we need to deal with the case when the old markers
! straddle x=0.
  if (print_node) then
    write(irefwr,*) 'PERROR(remarker_cyl): not yet suited for'
    write(irefwr,*) 'full cylinder geometry'
  endif
  call instop
endif

ipm = 0
ipn = 0
if (nochain.ne.1) stop 'remarker_cyl nochain<>1'
inewtot=1
do ichain=1,nochain
   nmark = imark(ichain)
   ip = 2*(ipn+1)-1
   tracer(ipn+1)%x = tracer(ipn+1)%xnew
   tracer(ipn+1)%y = tracer(ipn+1)%ynew
   inew=1
   iold=1
10       continue
     inew = inew+1
     inewtot = inewtot + 1
     ip = 2*(ipm+iold)-1
     x1 = tracer(ipm+iold)%xnew
     y1 = tracer(ipm+iold)%ynew
     iold = iold+1
100        continue
      ip = 2*(ipm+iold)-1
      x2 = tracer(ipm+iold)%xnew
      y2 = tracer(ipm+iold)%ynew
      distance  = sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1) )
      if (distance.lt.dm_min(ichain)) then
!       tracer can be removed. Continue along the chain
!       until distance has become large enough
        if (iold.lt.nmark-1) then
           iold = iold+1
           goto 100
        endif
      endif
        
     if (distance.gt.dm(ichain)) then
!       add a marker
        xnew = 0.5*(x1+x2)
        ynew = 0.5*(y1+y2)
        rnew = sqrt(xnew*xnew+ynew*ynew)
        if (rnew.ge.rtop_threshold.or.rnew.le.rbot_threshold) then
!          use cylindrical coordinates for interpolation
           rm1 = sqrt(x1*x1+y1*y1)
           rm2 = sqrt(x2*x2+y2*y2)
           theta1 = acos(y1/rm1)
           theta2 = acos(y2/rm2)
           thetanew = 0.5*(theta1+theta2)
           rnew = 0.5*(rm1+rm2)
           xnew = rnew*sin(thetanew)
           ynew = rnew*cos(thetanew)
        else
          call checkbounds(0,xnew,ynew)
        endif
        ip = 2*(ipn+inew)-1
        tracer(ipn+inew)%x   = xnew
        tracer(ipn+inew)%y   = ynew
        inewtot=inewtot+1
        inew=inew+1
     endif

     if (inewtot.gt.NTRACMAX) then
!       markerchain has become to big: routine failed
        if (print_node) then
           write(irefwr,*) 'PERROR(remarker_cyl): markerchain too big'
           write(irefwr,*) 'inewtot = ',inewtot
           write(irefwr,*) 'NTRACMAX = ',NTRACMAX
        endif
        call instop
     endif

     ip = 2*(ipn+inew)-1
     tracer(ipn+inew)%x   = x2
     tracer(ipn+inew)%y = y2
   if (iold.le.nmark-1) goto 10

   nmnew(ichain) = inew
   ipm = ipm+nmark
   ipn = ipn+inew
 enddo

 do ichain=1,nochain
    imark(ichain) = nmnew(ichain)
 enddo

 nmark = imark(1)
 xi_eta_stored=.false.
 end subroutine remarker_cyl


!   REMARKER_CART
!
!   Regrid the markerchain. Check if the interval between
!   markers has become larger than the tolerance DMAX. If so,
!   add a marker in the center of this interval using linear
!   interpolation.
!
!   PvK 960830
!
!   PvK 072604
subroutine remarker_cart
use sepmodulecomio
use mpetrac
use tracers
use control
implicit none
integer :: nmnew(10)
integer :: ipm,ipn,nmark,inew,iold,ip,inewtot
real(kind=8) :: dx,dy,rl,x1,x2,y1,y2,distance,xnew,ynew

ipm = 0
ipn = 0
inewtot=1
do ichain=1,nochain
   nmark = imark(ichain)
   ip = 2*(ipn+1)-1
   tracer(ipn+1)%x=tracer(ipn+1)%xnew
   tracer(ipn+1)%y=tracer(ipn+1)%ynew
   inew=1
   do iold=1,nmark-1
     inew = inew+1
     inewtot = inewtot+1
     ip = 2*(ipm+iold)-1
     x1 = tracer(ipm+iold)%xnew
     y1 = tracer(ipm+iold)%ynew
     x2 = tracer(ipm+iold+1)%xnew
     y2 = tracer(ipm+iold+1)%ynew
     distance  = sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1) )
     if (distance.gt.dm(ichain)) then
        xnew = 0.5*(x1+x2)
        ynew = 0.5*(y1+y2)
        ip = 2*(ipn+inew)-1
        tracer(ipn+inew)%x = xnew
        tracer(ipn+inew)%y = ynew
        inew=inew+1
        inewtot = inewtot+1
     endif
     if (inewtot.gt.NTRACMAX) then
!       markerchain has become to big: routine failed
        if (print_node) then
           write(irefwr,*) 'PERROR(remarker_cart): markerchain too big'
           write(irefwr,*) 'inewtot, NTRACMAX: ',inewtot,NTRACMAX
        endif
        call instop
     endif
     ip = 2*(ipn+inew)-1
     tracer(ipn+inew)%x = x2
     tracer(ipn+inew)%y = y2
   enddo
   nmnew(ichain) = inew
   ipm = ipm+nmark
   ipn = ipn+inew
 enddo

 do ichain=1,nochain
    imark(ichain) = nmnew(ichain)
 enddo
! write(irefwr,*) 'remarker_cart: ', imark(1),tracer(1)%x,tracer(imark(1))%x
 xi_eta_stored=.false.

 end subroutine remarker_cart



