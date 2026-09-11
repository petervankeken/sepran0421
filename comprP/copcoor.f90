! COPCOOR
! PvK 990413
subroutine copcoor
use mpetrac
use tracers
implicit none

integer :: nmark,i,j,ip1,ip2,ntot

if (itracoption == 2) then
   ! markerchain method
   nmark=0
   do ichain=1,nochain
      nmark=nmark+imark(ichain)
   enddo
   ntot = nmark
else ! includes itracoption==3
   ntot=0
   do idist=1,ndist
      ntot=ntot+ntrac(idist)
   enddo
endif


do i = 1,ntot
   ip1 = 2*i-1
   ip2 = 2*i
   tracer(i)%x = tracer(i)%xnew
   tracer(i)%y = tracer(i)%ynew
enddo

return
end

