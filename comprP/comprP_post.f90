subroutine comprP_post
use sepmoduleoldrouts
use sepmodulecomio  ! provides lu's for reading & writing
use sepmodulecpack  ! parallel computing interface
use sepran_arrays ! access to kmesh, kprob, isol, etc. and interfaces to subroutines
use sepran_interface
use control
use mtime
use dtm_elem  ! controls building of pressure mass matrix
use coeff
use tracers
implicit none
integer :: ndim,nphys,npoint,nunkp
integer :: modelv,mconv,nparm,nsup1,idum,iread,ihelp,icurvs(2)
integer :: idegfd,i,ichoice,isubr,ncntln
real(kind=8) :: rayleigh,umax,outputtime,contln(10)
character(len=120) :: namedof(4)
integer :: ichvc,iprob_here,ip,ipbuoy2,ipuser,map2(5),jchoice=0
real(kind=8) :: q_a(20),vrms,pee,quu,solmin,solmax
logical :: first=.true.,plot_vel
real(kind=4) :: cput
integer :: ih,im,is,i1=0,iplots(4),ncoorc(200),jinput(200),ivec
real(kind=8) :: rplots(20),coorcn(200),C_int
character(len=80) :: plottext
integer :: ipoststart,ipostend,ipoststep


open(99,file='postprocessing')
read(99,*) ipoststart,ipostend,ipoststep
close(99)

open(99,file='tC.dat')
jchoice=10 ! first time calling; do averages
do inout=ipoststart,ipostend,ipoststep
   write(ourplotname,'(''solutions/T.'',i4.4)') inout
   write(irefwr,*) 'ourplotname: ',ourplotname
   call readbs_netcdf(ourplotname,isol2)
   write(ourplotname,'(''T.'',i4.4)') inout
   call toGMT(jchoice,rlampix,nres_GMT,ourplotname)
   jchoice=11
!  if (output_velocity_solution) then
!     write(ourplotname,'(''solutions/UV.'',i4.4)') inout
!     call readbs_netcdf(ourplotname,isol1,kmesh1,kprob1)
!  endif
   write(irefwr,*) itracoption
   if (itracoption>0 .and. Coutput) then
      call readtrac_netcdf(inout)
      call tracdens()
      call find_C_int(0,iuser_here,user_here,C_int)
      write(99,*) inout*dtoutp,C_int
   endif
enddo
close(99)


end subroutine comprP_post
