! Single Stokes solve with a fixed blob of buoyancy +/- phase change
subroutine blobby
use sepmodulecomio
use sepmoduleoldrouts
use control
use geometry
use coeff
use sepran_arrays
implicit none
!integer,parameter :: NITERMAX=20
integer :: iprob_here,ichvc,idum,ncntln=0,ipoint=0,i1=0,iinvec(2)
real(kind=8) :: contln(11),p=0.0_8,q=0.0_8,tmin,tmax,vrms,vrms_old,T_bot_here,rinvec(2)

! Create a double Gaussian blob for T'
iprob_here=2
ichvc=1
iu1lc(1)=14 
call creavc(0,ichvc+(iprob_here-1)*1000,idum,isol2,kmesh1,kprob1,iu1lc,u1lc)
! form T=T'+Tbar
call algebr(6,1,iadia,i1,i1,kmesh1,kprob1,tmin,tmax,p,q,ipoint)
write(6,*) 'blobby: adiabatic min/max: ',tmin,tmax
iinvec(1)=2
iinvec(2)=27
rinvec(1)=1.0
rinvec(2)=1.0
call manvec(iinvec,rinvec,isol2,iadia,isol2,kmesh1,kprob1)
! reset bottom boundary condition just in case
T_bot_here=eos_T(eos_np)
call bvalue(0,1002,kmesh1,kprob1,isol2,T_bot_here,ibottom,ibottom,1,0)
call algebr(6,1,isol2,i1,i1,kmesh1,kprob1,tmin,tmax,p,q,ipoint)
write(6,*) 'blob T min/max, T_bot_here: ',tmin,tmax,T_bot_here

! solve Stokes equation
vrms_old=0.0_8
niter=0
do 
   niter=niter+1
   call stokesP_new(1)
   call pevrms(vrms,isol1)
   if (compress.and..not.TALA) then
      dif=abs(vrms_old-vrms)/vrms
      vrms_old=vrms
      if (print_node) write(irefwr,'(i5,f12.4,e15.7)') niter,vrms,dif
      if (dif>1e-8.and.niter<NITERMAX) then
        cycle
      else
        exit
      endif
   else 
      exit
   endif
enddo
write(6,*) 'vrms = ',vrms
inout=1
call bmout()



if (print_node) then
   ncntln=-100
   write(ourplotname,'(''PLOTS/BLOB.0000'')')
   call plotc1(1,kmesh1,kprob1,isol2,contln,ncntln,15.0_8,1.0_8,1)
   write(irefwr,*) 'stopping in blobby'
endif
call instop
end subroutine blobby
