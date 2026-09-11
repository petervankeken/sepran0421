!    CORRH
! 
!    Solves    .
!            M T  + S T  = f
!    using a CN scheme. The system is non-linear, S depends on time.
! 
!    Parameters:
!       kmesh     i   Standard Sepran array
!       kprob     i   Standard Sepran array
!       intmat    i   Information of the large matrix for this problem
!       isol      o   Solution vector T(n+1)
!       islold    i   Old solution vector T(n): output of predictorstep
!       (i)user   i   Contains information on coefficients
!       numold    i   number of solutionvectors in ivcold
!       ivcold    i   'old solution array' used to transport 
!                     information between problems
!       matrs     i   Stiffness matrix S(n)
!       matrm     i   Mass matrix (assumed to be constant in time)
!       irhsd     i   Righthandside vector f(n)
!  
!    290889 PvK
!    Adapted to 0619 PvK 070519
subroutine corrh(kmesh,kprob,intmat,isol,isolold,user,iuser,matrs,matrm,irhsd)
use sepmoduleoldrouts
use sepmodulecomio
use coeff
use control
use mtime
implicit none
integer :: kmesh,kprob,intmat,isol,isolold,matrs,matrm,irhsd
integer :: iuser(*)
real(kind=8) :: user(*)
integer :: matrs1=0,ivec1=0,ivec2=0,ivec3=0,ifirst
integer :: inpsol(20),iread
real(kind=8) :: rinsol(10),bee,t,DONE=1.0_8,get_resmem
character(len=10) :: tname
real(kind=4) :: t11,t12

integer :: iinbld(20),ichoice,i
real(kind=8) :: t2step,p,q,alpha2,beta2
save ifirst,iinbld,matrs1,inpsol,rinsol
save ivec1,ivec2,ivec3


t11=second()
t2step = tstepp*0.5d0

call pefilcof(2,iuser,user)

! Built matrices and rhs vector at time t
! Stiffness matrix and
iinbld=0
if (idia==2) then
  iinbld(1) = 4
  iinbld(2) = 1
  iinbld(3) = 1
  iinbld(4) = idia
else
  iinbld(1) = 3
  iinbld(2) = 1
  iinbld(3) = 0
endif
if (print_node.and.red_flag) write(irefwr,*) '     build: ',get_resmem()
call build(iinbld,matrs1,intmat,kmesh,kprob,irhsd,matrm,isol,isolold,iuser,user)
! then add, if necessary, the contribution of the tracer heating
if (iqtype == 5) then
   if (print_node) write(irefwr,*) 'PERROR(corrh): not yet ready for iqtype==5'
   call instop
!  call tracer_heat_rhsd(kmesh1,kprob1,kmesh,kprob,irhsd,2)
endif

     
! ivec1 = M*uold  add it to irhsd*dt
! write(irefwr,*) 'matrm: ',(matrm(i),i=1,5),idia
ichoice=idia+3
if (print_node.and.red_flag) write(irefwr,*) '     maver: ',get_resmem()
call maver(matrm,isolold,ivec1,intmat,kprob,ichoice)
if (print_node.and.red_flag) write(irefwr,*) '     algebr: ',get_resmem()
call algebr(3,0,irhsd,ivec1,irhsd,kmesh,kprob,tstepp,DONE,p,q,i)

! S1 = M + dt/2*S1
if (print_node.and.red_flag) write(irefwr,*) '     sepaddmat: ',get_resmem()
!call addmat(kprob,matrs1,matrm,intmat,t2step,alpha2,1.0_8,beta2)
call sepaddmat(matrs1,matrm,t2step,alpha2,DONE,beta2)

! ivec1 = S1*usol(p) ; add it to irshd
! write(irefwr,*) 'matrs1: ',(matrs1(i),i=1,5),idia
if (print_node.and.red_flag) write(irefwr,*) '    maver: ',get_resmem()
call maver(matrs1,isol,ivec1,intmat,kprob,6)
bee = -1.0_8
if (print_node.and.red_flag) write(irefwr,*) '    algebr: ',get_resmem()
call algebr(3,0,irhsd,ivec1,irhsd,kmesh,kprob,DONE,bee,p,q,i)

! ivec1 = S*uold ; add it to irhsd
if (print_node.and.red_flag) write(irefwr,*) '    maver: ',get_resmem()
call maver(matrs,isolold,ivec1,intmat,kprob,5)
t2step=-t2step
if (print_node.and.red_flag) write(irefwr,*) '    algebr: ',get_resmem()
call algebr(3,0,irhsd,ivec1,irhsd,kmesh,kprob,DONE,t2step,p,q,i)

if (isolmethod8 <= 0) then
   inpsol(1) = 3
!  Direct solution method with profile storage; not positive definite
   inpsol(2) = 0
!  solution method
   inpsol(3) = 0
else
   inpsol(1) = 14
!  ipos
   inpsol(2) = 0
!  solution_method 1=CG 2=CGS 3=GMRES 4=GMRESR
   inpsol(3) = isolmethod8
!  ipreco 1=diag 2=eisenstat 3=ILU
   inpsol(4) = ipreco8
   inpsol(5) = maxiter8
   inpsol(6) = iprint8
!  write(irefwr,*) 'ipreco8 etc: ',inpsol(3:6),cgeps8
!  matrix is symmetric(1) or not (0)?
   inpsol(7) = 0
!  dimension of Krylov space for GMRES
   inpsol(8) = 20
!  ISTART: 1=starts with given vector
   inpsol(9) = 1 
!  KEEP: 0=destroy preconditioning matrix; 
   inpsol(10)= 0
   inpsol(11)= ireler8  ! ireler
   inpsol(12)= 0
!  NTRUNC: for GMRESR: maximum number of search directions
   inpsol(13)= 5
!  NINNER: for GMRESR: maximum number of inner loop iters
   inpsol(14)= 5
!  EPS: required accuracy for iteration process
   rinsol(1) = cgeps8
endif
iread=-1
if (print_node.and.red_flag) write(irefwr,*) '    solve: ',get_resmem()
call solvel(inpsol,rinsol,matrs1,isol,irhsd,intmat,kmesh,kprob,iread)
t12=second()
if (petiming .and. print_node) write(irefwr,*) 'time in corrh: ',t12-t11

end subroutine corrh
