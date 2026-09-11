subroutine peplotit
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
real(kind=8) :: rplots(20),coorcn(200)
character(len=80) :: plottext

save ncoorc,iplots,coorcn,rplots,jinput,map2

if (.not.print_node) return


if (print_node) then
   call algebr(6,1,isol1,i1,i1,kmesh1,kprob1,solmin,solmax,pee,quu,ip)
   write(irefwr,'(''umin/max       : '',2e15.7)') solmin,solmax
   call algebr(6,2,isol1,i1,i1,kmesh1,kprob1,solmin,solmax,pee,quu,ip)
   write(irefwr,'(''vmin/max       : '',2e15.7)') solmin,solmax
   if (itype_stokes==903) then
      call algebr(6,3,isol1,i1,i1,kmesh1,kprob1,solmin,solmax,pee,quu,ip)
      write(irefwr,'(''pmin/max       : '',2e15.7)') solmin,solmax
   endif
   call algebr(6,1,isol2,i1,i1,kmesh1,kprob1,solmin,solmax,pee,quu,ip)
   write(irefwr,'(''Tmin/max       : '',2e15.7)') solmin,solmax

   if (inout<0) inout=1000+inout
   iplots(1)=2
   iplots(2)=1
   rplots(1)=15.0_8
   rplots(2)=1.0_8
   write(ourplotname,'(''PLOTS/TEMP.'',i4.4)') inout
   call plotcn(kmesh1,kprob1,isol2,iplots,rplots,contln,jinput,plottext,ncoorc,coorcn)

   inquire(file='plot_velocity',exist=plot_vel)
   if (plot_vel) then
      write(ourplotname,'(''PLOTS/VEL.'',i4.4)') inout
      write(plottext,'(''V: '',i4)') inout
      call plotvc(1,2,isol1,isol1,kmesh1,kprob1,15.0_8,1.0_8,0.0_8)
   endif
endif
!write(irefwr,*) 'compute and plot stream'
call stream(1,ivec,istream,0,0,kmesh1,kprob1,isol1)
write(ourplotname,'(''PLOTS/STREAM.'',i4.4)') inout
write(plottext,'(''STREAM: '',i4)') inout
call plotcn(kmesh1,kprob1,istream,iplots,rplots,contln,jinput,plottext,ncoorc,coorcn)

if (fieldC) then
   write(ourplotname,'(''PLOTS/C.'',i4.4)') inout
   write(plottext,'(''C: '',i4)') inout
   call plotcn(kmesh1,kprob1,idens,iplots,rplots,contln,jinput,plottext,ncoorc,coorcn)
endif
!write(irefwr,*) fieldC,ourplotname
!call instop

if (idivv>0) then
   write(ourplotname,'(''PLOTS/DIVV.'',i4.4)') inout
   write(plottext,'(''DIV U: '',i4)') inout
   call plotcn(kmesh1,kprob1,idivv,iplots,rplots,contln,jinput,plottext,ncoorc,coorcn)
endif
if (ipress>0) then
   write(ourplotname,'(''PLOTS/P.'',i4.4)') inout
   write(plottext,'(''P: '',i4)') inout
   call plotcn(kmesh1,kprob1,ipress,iplots,rplots,contln,jinput,plottext,ncoorc,coorcn)
endif
if (icurlv>0) then
   write(ourplotname,'(''PLOTS/CURLV.'',i4.4)') inout
   write(plottext,'(''CURL V: '',i4)') inout
   call plotcn(kmesh1,kprob1,icurlv,iplots,rplots,contln,jinput,plottext,ncoorc,coorcn)
endif

! something odd in sepran 0220: when using intcoor matrix info is overwritten
! something odd in sepran 1020: when using intcoor kpropb info is overwritten
if (.not.noGMT) then
  write(ourplotname,'(''T.'',i4.4)') inout
  jchoice=10+jchoice
  call toGMT(jchoice,rlampix,nres_GMT,ourplotname)
endif


if (.not.no_temperature_solution) then
  write(ourplotname,'(''VTK/T_'',i4.4)') inout
  !write(irefwr,*) 'create ',ourplotname
  ichoice=1
  idegfd=0
  namedof(1)='temperature'
  call sol2vtu(1,isol2,ourplotname,0,namedof)
endif
namedof(1)='velocity'
write(ourplotname,'(''VTK/UV_'',i4.4)') inout
!write(irefwr,*) 'create ',ourplotname
call sol2vtu(2,isol1,ourplotname,0,namedof)
if (itypv>0) then
  namedof(1)='viscosity'
  write(ourplotname,'(''VTK/eta_'',i4.4)') inout
  !write(irefwr,*) 'create ',ourplotname
  call sol2vtu(1,ivisc,ourplotname,1,namedof)
endif
if (Di>0) then
  write(ourplotname,'(''VTK/PHI_'',i4.4)') inout
  ichoice=1
  idegfd=0
  namedof(1)='phi'
  write(6,*) 'write iphi to vtu',ourplotname
  write(6,*) 'iphi/iphi2/ivisdip: ',iphi,iphi2,ivisdip
  call sol2vtu(1,iphi,ourplotname,0,namedof)
  call system('ls VTK')
endif

if (fieldC) then
  write(ourplotname,'(''VTK/C_'',i4.4)') inout
  !write(irefwr,*) 'create ',ourplotname
  ichoice=1
  idegfd=0
  namedof(1)='C'
  call sol2vtu(1,idens,ourplotname,0,namedof)
endif

write(ourplotname,'(''solutions/T.'',i4.4)') inout
call writbs_netcdf(ourplotname,isol2,kmesh1,kprob1)
if (output_velocity_solution) then
   write(ourplotname,'(''solutions/UV.'',i4.4)') inout
   call writbs_netcdf(ourplotname,isol1,kmesh1,kprob1)
endif

end subroutine peplotit
