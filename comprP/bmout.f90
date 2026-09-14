subroutine bmout()
use sepmoduleoldrouts
use sepmodulecomio  ! provides lu's for reading & writing
use sepmodulecpack  ! parallel computing interface
use sepmodulekmesh
use sepmodulevecs
use sepran_arrays ! access to kmesh, kprob, isol, etc. and interfaces to subroutines
use sepran_interface
use geometry
use control
use coeff
use convparam
use dtm_elem  ! controls building of pressure mass matrix
use mdebugTbars
implicit none
integer :: nusol
integer :: modelv,mconv,nparm,nsup1,idum,iread,ihelp,icurvs(2)
integer :: idegfd,i,ichoice,isubr,ncntln,noutje,iinder(8),ielhlp
real(kind=8) :: rayleigh,umax,outputtime,contln(10),veloc(10),veloc_d(10),vismin,vismax
real(kind=8) :: temp00,temp1,gnus1,gnus2,q1,q2,q3,q4,avT(10),avwork(10),avphi(10),avvisdip(10)
real(kind=8) :: pee,quu,f_cyl,vrmsa(4),pressure_max,curlv_max,tempph
character(len=120) :: namedof(4)
integer :: ichvc,iprob_here,ip,ipoint,iu1(2)
real(kind=8) :: q_a(4),vrms,dNu,dq,dw,q_top,q_bot,Rah
real(kind=8) :: plate_mobility,plateness,plasticity_integral,davphi,rinvec(20)
logical :: output_pressure

call pevrms(vrms,isol1)
if (tackley.or.(tosi15.and.(tosi15_case==2.or.tosi15_case==4))) then
   call create_plasticity
   call find_intP(plasticity_integral)
endif
if (ibench_type==300) then
   ! just Stokes solve with cpu output
   if (print_node) write(irefwr,'(2i15,f15.7,f12.2)') nstokes_solves,2*npoint,vrms,cpu_stokes
   call instop
endif

icurvs(1)=0
icurvs(2)=3
funcx(1)=20005
funcy(1)=20005
ichoice=0
call compcr(ichoice,kmesh1,kprob1,isol1,1,icurvs,funcx,funcy)
umax=maxval(funcy(6:5+int(funcy(5))))
!if (print_node) write(irefwr,'(''umax: '',e15.7)') umax
call nusseltP(q_a)
temp1=q_a(3)
temp00=q_a(4)
gnus1=q_a(1)/(temp1-temp00)  ! update for delta1K and delta1C
gnus2=q_a(2)/(temp1-temp00)
if (print_node) write(irefwr,*) 'temp1, temp0 = ',temp1,temp00
if (compute_curl_v) then
   call sepactsolbf1(isol1)
   call pefilcof(1,iuser_here,user_here)
   iinder(1)=8
   iinder(2)=1
   iinder(3)=0
   iinder(4)=5  ! curl v 
   iinder(5)=0
   iinder(6)=0
   iinder(7)=0
   iinder(8)=2
   ! ipresse is special vector containing value of pressure in the elements (centroids)
   call deriv(iinder,icurlv,kmesh1,kprob1,isol1,iuser_here,user_here)
   curlv_max = anorm(1,3,1,kmesh1,kprob1,icurlv,icurlv,ielhlp)
   if (print_node) write(irefwr,*) 'curlv_max: ',curlv_max
endif


!call check_stress_tensor(iuser_here,user_here)
output_pressure=.false.
if (itype_stokes==900.and.output_pressure) then
   call sepactsolbf1(isol1)
   call pefilcof(1,iuser_here,user_here)
   iinder(1)=8
   iinder(2)=1
   iinder(3)=0
   iinder(4)=20  ! pressure in the elements
   iinder(5)=0
   iinder(6)=0
   iinder(7)=0
   iinder(8)=2
   ! ipresse is special vector containing value of pressure in the elements (centroids)
   call deriv(iinder,ipresse,kmesh1,kprob1,isol1,iuser_here,user_here)
   pressure_max= anorm(1,1,1,kmesh1,kprob1,ipresse,ipresse,ielhlp)
   ! shift so that P is all positive (to match 903 with essential b.c.)
   call print_and_shift_P(ks(ipresse)%sol,nelem)
   if (print_node) write(irefwr,*) 'maximum pressure 900: ',pressure_max
   ! distribute values to nodal points
   if (ipress==0) then
      iu1=0
       u1=0
      call creavc(0,1001,1,ipress,kmesh1,kprob1,iu1,u1,iu1,u1)
   endif
   call convert_elem_to_coor(npoint,nelem,ndim,npelm,coor,kmeshc,ks(ipresse)%sol,ks(ipress)%sol)
   !call prinrv(ipress,kmesh1,kprob1,1,1,'work')
    
else if (itype_stokes==903.and.output_pressure) then
   call coefpress
   call print_and_shift_P(ks(ipress)%sol,npoint)
endif
if (print_node) write(irefwr,*) 'before peplotit'
call peplotit()

if (print_node) then
   open(9,file='comprP.diagnostics')
   write(9,'(2f12.3)') q_a(1),umax
   close(9)
endif

call writbs_netcdf('T_restart.nf',isol2)
call writbs_netcdf('UV_restart.nf',isol1)

noutje=1
call surfacevel(kmesh1,kprob1,isol1,veloc,veloc_d,noutje)
if (tosi15.or.tackley) then
   call find_mobility_plateness_cart(vrms,veloc(3),plate_mobility,plateness)
   if (print_node) then
     open(9,file='plates.dat')
     write(irefwr,*)
     write(irefwr,'(5a15)') 'Ra','sigma_y','<plasticity>','mobility','plateness'
     write(9,'(5e15.7)') Ra,sigma_y,plasticity_integral,plate_mobility,plateness
     write(irefwr,'(5e15.7)') Ra,sigma_y,plasticity_integral,plate_mobility,plateness
     write(9,'(5a15)') 'Ra','sigma_y','<plasticity>','mobility','plateness'
     close(9)
   endif
endif
call averageT(1,iuser_here,user_here,avT)
if (Di > 0.or.compress) then
   call getvisdip
   !call prinrv(iphi2,kmesh1,kprob1,1,1,'work')
   call averageT(2,iuser_here,user_here,avphi)
   !call compressionwork(kmesh1,kprob1,isol1,isol2,icompwork)
   call compressionwork()
   !call prinrv(icompwork,kmesh1,kprob1,1,1,'work')
   call averageT(3,iuser_here,user_here,avwork)
   call phifromgrad(kmesh1,kprob1,isol1,ivisdip)
   call averageT(2,iuser_here,user_here,avvisdip)
   avvisdip(1)=avvisdip(1)*DiRa
   avphi(1) = avphi(1)*DiRa  ! already divided by volume
   avwork(1) = avwork(1)*DiRa ! same; note division with Ra in Di/Ra 
   if (print_node) then
      write(irefwr,'(''<phi>     : '',1e15.7)') avphi(1)! avphi(1)*DiRa,DiRa
      write(irefwr,'(''<work>    : '',1e15.7,f12.2)') avwork(1),100*(avphi(1)-avwork(1))/avphi(1)!,avwork(1)*Di
      write(irefwr,'(''<avvisdip>: '',1e15.7)') avvisdip(1)!,avvisdip(1)*Di
   endif
endif

if (tosi15) then
   ! get minimum and maximum viscosity
   call algebr(6,1,ivisc,ivisc,ivisc,kmesh1,kprob1,vismin,vismax,pee,quu,ipoint)
endif

if (print_node) then
   write(irefwr,*)
   write(irefwr,'(''Niter     '',i15,f15.2)') niter,cpu_total
   write(irefwr,'(''Nu, <T>   '',4f15.7,2i5)') gnus1,gnus2,temp1,temp00,itop,ibottom
   write(irefwr,'(''Vrms      '',f15.7,f15.7)') vrms ! ,vrms_d
   write(irefwr,'(''avT       '',f15.7,f15.3,f15.7)') avT(1),avT(2),avT(1)/deltaT_dim

   write(irefwr,'(''Maximum surface velocity '',2f12.7)') veloc(1),veloc_d(1)
   write(irefwr,'(''Average surface velocity '',2f12.7)') veloc(2),veloc_d(2)
   write(irefwr,'(''RMS surface velocity     '',f12.7)') veloc(3)

   write(irefwr,*) 'Difference in heat flow through top and bottom: '
   write(irefwr,'(2f12.3,$)') gnus1*surf(2),gnus2*surf(1)
   if (gnus1 > 0) then
      dNu = (gnus1*surf(2)-gnus2*surf(1))/(gnus1*surf(2))*100
      write(irefwr,'(f12.3)') dNu
   else
      write(irefwr,*)
   endif

   ! final output for various benchmarks
   write(irefwr,*) 'ibench_type=',ibench_type,tosi15
   if (ibench_type>0.or.tosi15) then
      if (did_not_converge) then
        write(irefwr,'(e12.5,f10.4,L5,'' did not converge '')') Ra,Di,compress
      else if (Di.gt.0) then
!       if (delta1K) then
!         write(irefwr,100) Ra,Di,compress,gnus1,gnus2,vrms,veloc(1),veloc(2),avT(1)/deltaT_dim,avphi(1), &
!      &     avwork(1),(avphi(1)-avwork(1))/avwork(1)*100,dNu
!         write(irefwr,101) Ra,Di,compress,gnus1,gnus2,vrms,veloc(1),veloc(2),avT(1)/deltaT_dim,avphi(1), &
!      &     avwork(1),(avphi(1)-avwork(1))/avwork(1)*100,dNu
!       else
        ! make delta1K adjustment in averagetemp() instead except for <T>
          if (delta1K) then
             avT(1)=avT(1)/deltaT_dim
             Rah=Ra*deltaT_dim
          else
             Rah=Ra
          endif
          write(irefwr,100) Rah,Di,compress,gnus1,gnus2,vrms,veloc(1),veloc(2),avT(1),avphi(1), &
       &     avwork(1),(avphi(1)-avwork(1))/avwork(1)*100,dNu
          write(irefwr,101) Rah,Di,compress,gnus1,gnus2,vrms,veloc(1),veloc(2),avT(1),avphi(1), &
       &     avwork(1),(avphi(1)-avwork(1))/avwork(1)*100,dNu
!       endif
     else if (tosi15) then
        write(irefwr,102) avT(1)/deltaT_dim,gnus1,gnus2,vrms,veloc(3),veloc(1),vismin,vismax,avwork(1), &
     &     avphi(1),(avphi(1)-avwork(1))/(max(avwork(1),avphi(1)))*100
     else
        write(irefwr,100) Ra,Di,compress,gnus1,gnus2,vrms,veloc(1),veloc(2),avT(1),0.,0.,0.,dNu
        write(irefwr,101) Ra,Di,compress,gnus1,gnus2,vrms,veloc(1),veloc(2),avT(1),0.,0.,0.,dNu
     endif
   endif
endif
100   format(e12.5,f10.6,L5,2(f14.6),3(f14.6),3(f12.5),f10.3,f10.3)
101   format(e12.5,',',f10.6,',',L5,',',2(f14.6,','),3(f14.6,','),3(f12.5,','),f10.3,',',f10.3)
102   format(6f12.4,2g15.7,2f12.4,f5.2,'%')


if (ibench_type==7 .and. print_node) then
  ! cylindrical benchmark
  f_cyl=r1/r2
  gnus1=abs(gnus1*log(f_cyl)/(1-f_cyl))
  gnus2=abs(gnus2*f_cyl*log(f_cyl)/(1-f_cyl))
  call vrms_r_theta(isol1,vrmsa)
  write(irefwr,'(''vrms       '',4f15.7)') vrmsa(1:4)
  write(irefwr,'(3f8.3,5f12.3,5f8.3)') avT(1),gnus1,gnus2,vrmsa(1),vrmsa(3),vrmsa(2),veloc(1),veloc(4),avphi(1),avwork(1)
endif

if (ibench_type==17 .and. print_node) then
   ! this computes vrms based on polar coordinates
   ! use standard Cartesian instead
   !call vrms_r_theta(isol1,vrmsa)
   q_top=gnus1*surf(2)
   q_bot=gnus2*surf(1)
   dq=(q_top-q_bot-volume*q_layer(1))/q_top
   if (Di>0) then
      dw=(avphi(1)-avwork(1))/avwork(1)
   else
      dw=0
   endif
   write(irefwr,'(e15.4,f12.3,8f12.5,2f6.3)') Ra,Di,q_layer(1),q_top,q_bot,avT(1),vrms,veloc(1),avphi(1),avwork(1),dq*100,dw*100
   write(irefwr,'(e15.4,'','',f12.3,'','',f12.5,'',''i5,'','',7(f12.5,'',''),2(f6.3,'',''))') Ra,Di,q_layer(1),&
      & eos_type,q_top,q_bot,avT(1), &
      &  vrms,veloc(1),avphi(1),avwork(1),dq*100,dw*100
endif

end subroutine bmout

subroutine print_and_shift_P(usol,N)
use geometry
implicit none
integer :: N
real(kind=8) :: usol(*)
integer :: i
real(kind=8) :: L2_P=0.0_8,Pmax=0.0_8,Pmin=0.0_8

!Pmax=maxval(usol(1:N))
!Pmin=minval(usol(1:N))
!L2_P=0.0_8
do i=1,N
   usol(i)=usol(i)-Pmin
   !write(irefwr,*) 'Pshift: ',usol(i)
   !L2_P=L2_P+usol(i)*usol(i)
enddo
!Pmax=maxval(usol(1:N))
!Pmin=minval(usol(1:N))
!L2_P=sqrt(L2_P)/N
!write(irefwr,'(''min, max, L2: '',3e15.7)') Pmin,Pmax,L2_P

end subroutine print_and_shift_P

subroutine convert_elem_to_coor(npoint,nelem,ndim,npelm,coor,kmeshc,usole,usol)
use sepmodulecomio
use control
implicit none
integer :: npoint,nelem,ndim,npelm,kmeshc(*)
real(kind=8) :: coor(ndim,npoint),usole(nelem),usol(npoint)
real(kind=8),dimension(:),allocatable :: weights
integer :: ielem,allocate_status,ip,i,inod,inoda(6)
real(kind=8) :: b1,b3,c1,c3,xn(6),yn(6),weight,Pmin,Pmax

if (npelm /= 6) then
   if (print_node) write(irefwr,*) 'PERROR(convert_elem_to_coor): expecting 6 node triangle but npelm=' ,npelm
   call instop
endif

if (.not.allocated(weights)) then
   allocate(weights(npoint),stat=allocate_status)
   if (allocate_status/=0) then
      if (print_node) write(irefwr,*) 'PERROR(convert_elem_to_coor): allocate failed: ',allocate_status
      call instop
   endif
endif
weights=0.0_8
usol=0.0_8

!open(766,file='geom.dat')
!open(767,file='centroids.dat')
!ip=0
!do ielem=1,nelem
!  do i=1,npelm
!     inod=kmeshc(ip+i)
!     xn(i)=coor(1,inod)
!     yn(i)=coor(2,inod)
!     write(766,*) xn(i),yn(i)
!     !write(irefwr,*) inod,coor(1,inod),coor(2,inod)
!  enddo
!  ip=ip+npelm
!  !write(766,*) xn(1),yn(1)
!  !write(766,'(''>'')') 
!  !usole(ielem)=7000*(sum(yn(1:6))/6)
!  !write(767,*) sum(xn(1:6))/6,sum(yn(1:6))/6
!enddo
!close(766)
!close(767)
   

ip=0
! for test usole=1.0_8
!open(766,file='presse.dat')
do ielem=1,nelem
  !write(irefwr,*) 'element: ',ielem
  do i=1,npelm
     inod=kmeshc(ip+i)
     xn(i)=coor(1,inod)
     yn(i)=coor(2,inod)
     !write(irefwr,*) inod,coor(1,inod),coor(2,inod)
  enddo
  !write(766,*) sum(xn(1:6))/6,sum(yn(1:6))/6,usole(ielem)
  ! find measure of size of element 
  b1 = yn(3) - yn(5)
  b3 = yn(1) - yn(3)
  c1 = xn(5) - xn(3)
  c3 = xn(3) - xn(1)
  weight=-c3*b1+b3*c1
  ! write(irefwr,*) '    weight= ',weight
  ! distribute centroid value to vertices in element weighted by size
  do i=1,npelm
     inod=kmeshc(ip+i)
     usol(inod)=usol(inod)+usole(ielem)*weight
     weights(inod)=weights(inod)+weight
  enddo
  ip=ip+npelm
enddo
!close(766)


! now fill rest of nodal points
!ip=0
!do ielem=1,nelem
!   do i=1,npelm
!     inoda(i)=kmeshc(ip+i)   
!   enddo
!   usol(inoda(2))=0.5*( usol(inoda(1))+usol(inoda(3)) )
!   usol(inoda(4))=0.5*( usol(inoda(3))+usol(inoda(5)) )
!   usol(inoda(6))=0.5*( usol(inoda(5))+usol(inoda(1)) )
!   weights(inoda(2))=0.5*( weights(inoda(1)) + weights(inoda(3)) )
!   weights(inoda(4))=0.5*( weights(inoda(3)) + weights(inoda(5)) )
!   weights(inoda(6))=0.5*( weights(inoda(5)) + weights(inoda(1)) )
!   ip=ip+npelm
!enddo

Pmin=0.0_8
Pmax=0.0_8
!open(766,file='press.dat')
!open(767,file='weights.dat')
do i=1,npoint
   ! write(767,*) coor(1,i),coor(2,i),weights(i),usol(i)
   usol(i)=usol(i)/weights(i)
   ! write(766,*) coor(1,i),coor(2,i),usol(i)
   Pmin=min(Pmin,usol(i))
   Pmax=max(Pmax,usol(i))
enddo
!close(766)
!close(767)
!write(irefwr,*) 'min/max usol(3): ',Pmin,Pmax
   
if (allocated(weights)) deallocate(weights)

end subroutine convert_elem_to_coor
