!   IEULH
!
!   The temperature equation is solved using implicit euler
!   Information concerning the coefficients of the equation
!   are stored in (i)user and ivcold.
! 
!
!   Pvk  160589/290889
subroutine ieulh(kmesh,kprob,intmat,isol,isolold,user,iuser,matrs,matrm,irhsd)
use sepmoduleoldrouts
use sepmodulecomio
use coeff
use control
use mtime
implicit none
integer :: kmesh,kprob,intmat,isol,matrs,matrm,irhsd,isolold
integer :: iuser(*)
real(kind=8) :: user(*)
real(kind=4) :: t11,t12

integer :: ifirst,ivec1=0,ivec2=0,ivec3=0,matr1=0,ipser(100),mrhsd=0
integer :: iinbld(20),ichoice,ip
real(kind=8) :: pser(100),p,q,alpha2=0.0_8,beta2=0.0_8,bee,cee
real(kind=8) :: DONE,get_resmem
character(len=10) :: tname
integer :: inpsol(20),iread
real(kind=8) :: rinsol(20)
save ifirst,ivec1,ivec2,ivec3,matr1,ipser,pser,mrhsd
save iinbld
data ifirst/0/
data ipser(1),pser(1)/100,100/


t11=second()
call copyvc(isol,isolold)
call pefilcof(2,iuser,user)
!     *** Built matrices and rhs vector at time t
!     *** Stiffness matrix and 
iinbld=0
if (ifirst.eq.0.or.idia.eq.2) then
  iinbld(1) = 4
  iinbld(2) = 1
  iinbld(3) = 1
  iinbld(4) = idia
  ifirst = 1
else
  iinbld(1) = 3
  iinbld(2) = 1
  iinbld(3) = 0
endif
!write(irefwr,*) 'build with: ',iinbld(1:4)
if (print_node.and.red_flag) write(irefwr,*) '    build:',get_resmem()
call build(iinbld,matrs,intmat,kmesh,kprob,irhsd,matrm,isol,isolold,iuser,user)

! then add, if necessary, the contribution of the tracer heating
if (iqtype == 5) then
   if (print_node) write(irefwr,*) 'PERROR(ieulh): not yet ready for iqtype==5'
   call instop
!  write(irefwr,*) 'call tracer_heat_rhsd from ieulh'
!  call tracer_heat_rhsd(kmesh,kprob,kmesh2,kprob2,irhsd,1)
endif

! call getsepraninfomatr(matrs)
!call getsepraninfomatr(matrm)

if (printmatrix) then
   call prinmt(intmat,matrs,kprob)
   call instop
endif

! Built the system of equations
! M*uold
if (red_flag.and.print_node) write(irefwr,*) 'M*uold',get_resmem()
ichoice=idia+3
call maver(matrm,isolold,ivec1,intmat,kprob,ichoice)

! write(irefwr,*) 'ieulh: ',tstepp
! tstep*f + M*uold
if (red_flag.and.print_node) write(irefwr,*) '    algebr:',get_resmem()
DONE=1.0_8
call algebr(3,0,ivec1,irhsd,ivec3,kmesh,kprob,DONE,tstepp,p,q,ip)

! M + tstep*S
if (red_flag.and.print_node) write(irefwr,*) '    copymt:',get_resmem()
call copymt(matrs,matr1,kprob)
if (red_flag.and.print_node) write(irefwr,*) '    addmat:',get_resmem()
call addmat(kprob,matr1,matrm,intmat,tstepp,alpha2,DONE,beta2)
!call sepaddmat(matr1,matrm,tstepp,alpha2,DONE,beta2)

! ( M + tstep*S ) * usol(p)
if (print_node.and.red_flag) write(irefwr,*) '    maver:',get_resmem()
call maver(matr1,isol,ivec2,intmat,kprob,6)

! tstep*f + M*uold - (M+tstep*S) * usol(p)
bee = -1.0_8
if (print_node.and.red_flag) write(irefwr,*) '    algebr:',get_resmem()
call algebr(3,0,ivec3,ivec2,ivec1,kmesh,kprob,DONE,bee,p,q,ip)


inpsol=0
if (isolmethod8.le.0) then
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
   inpsol(3) = 1
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
if (pedebug.and.print_node) write(irefwr,*) 'call solvel with ',inpsol(1:10)
if (print_node.and.red_flag) write(irefwr,*) '    solvel:',get_resmem()
call solvel(inpsol,rinsol,matr1,isol,ivec1,intmat,kmesh,kprob,iread)
   
t12=second()
if (petiming .and. print_node) write(irefwr,*) 'time in ieulh: ',t12-t11

!if (petest) then
!   call prinrv(isol,kmesh,kprob,5,1,'temperature')
!endif
end subroutine ieulh

