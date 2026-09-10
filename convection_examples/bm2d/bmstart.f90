! PvK July 2024
subroutine bmstart()
use sepmodulekmesh
use sepmoduleoldrouts
use sepmodulevecs
use sepran_arrays
use control
use coeff
use geometry
implicit none
integer :: allocate_status,num_user_here,inputcr(10)=0
real(kind=8) :: rinputcr(10)=0.0_8,u1lc(10)=0.0_8
integer :: ielhlp=0,iu1lc(10)=0

cpustart=second()
t1=cpustart

! set defaults for parameters in fillcoef8/900 and update by namelist input
call readname_list()

! start sepran, read mesh & problem information, set up structure of matrices
call sepstr(kmesh,kprob,intmat)


call initialize_temp_vel()

call determine_geometry()

! create viscosity vector after temperature has been initialized
rinputcr=0.0_8
inputcr(1)=7
inputcr(2)=1 ! number of vectors to be created
inputcr(3)=2 ! ICHVC type of vector (2=special structure)
inputcr(4)=1 ! IPROB
inputcr(5)=1 ! IVEC type of special structure (1= one dof per nodal point)
inputcr(6)=0 ! ICOMPL
inputcr(7)=0 ! IFILL (0=create vector and set values to zero)
call creatv(kmesh,kprob,ivisc,inputcr,rinputcr)

! set up iadia, idens, iviscplast, ivisclin just for compatibility with future extension / merger with comprP
!inputcr(4)=2 ! iprob
 call creatv(kmesh,kprob,isecinv,inputcr,rinputcr)
 call creatv(kmesh,kprob,iviscplast,inputcr,rinputcr)
 call creatv(kmesh,kprob,ivisclin,inputcr,rinputcr)
 call creatv(kmesh,kprob,iplasticity,inputcr,rinputcr)
u1lc(1) = 0
 u1lc(1) = 0d0
! ichcrv = ichcr + (iprob-1)*1000
call creavc(0,1001,1,iphi,kmesh,kprob,iu1lc,u1lc,iu1lc,u1lc)
call creavc(0,1001,1,iphi2,kmesh,kprob,iu1lc,u1lc,iu1lc,u1lc)
call creavc(0,1001,1,ivisdip2,kmesh,kprob,iu1lc,u1lc,iu1lc,u1lc)
call creavc(0,1001,1,ivisdip,kmesh,kprob,iu1lc,u1lc,iu1lc,u1lc)
call creavc(0,1001,1,iadia,kmesh,kprob,iu1lc,u1lc,iu1lc,u1lc)
call creavc(0,1001,1,ipress,kmesh,kprob,iu1lc,u1lc,iu1lc,u1lc)
call creavc(0,1001,1,icompwork,kmesh,kprob,iu1lc,u1lc,iu1lc,u1lc)



! clearly I don't understand the creatv manpage because this doesn't work:
!inputcr(7)=3 ! IFILL, use next two positions to indicate what
!inputcr(8)=0 ! prescribe all degrees of freedom
!inputcr(9)=-1 ! fill with constant value of rinputcr(-1)
!inputcr(10)=0 ! use all nodes
!rinputcr(1)=1.0_8
!so just fill with zeros and then adjust
call creatv(kmesh,kprob,idens,inputcr,rinputcr)
call creatv(kmesh,kprob,irho,inputcr,rinputcr)
call creatv(kmesh,kprob,ialpha,inputcr,rinputcr)
call set_array_to_one(npoint,ks(idens)%sol)
call set_array_to_one(npoint,ks(irho)%sol)
call set_array_to_one(npoint,ks(ialpha)%sol)

! set up space for coefficients in user_here
num_user_here=10+7*npoint
allocate(user_here(num_user_here),stat=allocate_status)
if (allocate_status /= 0) then
   write(6,*) 'PERROR(simple): error in allocating user_here'
   call instop
endif
iuser_here(1)=NUSERMAX
user_here(1)=num_user_here

! set b.c. for both problems
call presdf(kmesh,kprob,isol)

end subroutine bmstart

subroutine determine_geometry()
use sepmodulecomio
use control
use coeff
use sepran_arrays
use geometry
implicit none
integer :: iuserh(100),ihelp
real(kind=8) :: userh(100),volint

if (cyl) then
   if (print_node) write(irefwr,*) 'PERROR(determine_geometry): not yet suited for cyl'
   call instop
endif
iuserh=0
 userh=0.0_8
iuserh(1)=100
 userh(1)=1.0_8*100

iuserh(2)=1
iuserh(6)=7
iuserh(8)=icoor800
iuserh(9)=0
iuserh(10)=-6
 userh(6)=1.0_8
volume=volint(0,1,1,kmesh,kprob,isol(2),iuserh,userh,ihelp)

if (print_node) then
   call system('mkdir -p GMT PLOTS VTK solutions Tracers')
endif

end subroutine determine_geometry

subroutine set_array_to_one(npoint,sepranvector)
implicit none
integer,intent(in) :: npoint
real(kind=8),intent(inout) :: sepranvector(*)
integer :: i

sepranvector(1:npoint) = 1e0_8

end subroutine set_array_to_one
