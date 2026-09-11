subroutine force_balance_setup()
use sepmodulecomio
use control
use sepran_arrays
use brandenburg
implicit none

real(kind=8) :: testtau,strchain_integrate
real(kind=8) :: rincre(10),format,factor,yfaccn
integer :: inpcre(10),isoltryfb
character(len=80) :: fname

!     include 'mysepar.inc'
!     include 'strchain.inc'
!     include 'cjp.inc'

! Make sure to use exact numbers from JP
plate_boundaries(1)=0.0_8
plate_boundaries(2)=0.7853975
plate_boundaries(3)=1.570795
plate_boundaries(4)=2.3561925
plate_boundaries(5)=3.14159
plate_boundaries(6)=3.9269875
plate_boundaries(7)=4.712385
plate_boundaries(8)=5.4977825
plate_boundaries(9)=6.28318


if (FBdebug.and.print_node) then 
  write(irefwr,*) 'force_balance_debug needs checking on func(4)'
  call instop
  ! Use velocity field for partrac as stress benchmark
  ! Synthetic velocity; create velocity vector
  ! Total number of entries
  inpcre(1)=10
  ! Number of vectors to be created
  inpcre(2)=1
  ! Type of vector  (1=solution vector)
  inpcre(3)=1
  ! problem number
  inpcre(4)=1
  ! sequence number of array of special structure
  inpcre(5)=0
  ! 0=real vector
  inpcre(6)=0
  ! IFILL: filling is defined by IFILL commands
  inpcre(7)=3
  ! fill first degree of freedom
  inpcre(8)=1
  ! fill using FUNC with ichoice=4
  inpcre(9)=4
  ! fill all nodes
  inpcre(10)=0
  call creatv(kmesh1,kprob1,isoltryfb,inpcre,rincre)
  ! fill second degree of freedom
  inpcre(8)=2
  ! fill using FUNC with ichoice=5
  inpcre(9)=5
  ! fill all nodes
  inpcre(10)=0
  call creatv(kmesh1,kprob1,isoltryfb,inpcre,rincre)
endif
call strchain_start()
call strchain_analyze()

if (FBdebug.and.print_node) then
   call strchain_compute(user_here,iuser_here,0)
   testtau = strchain_integrate(0)
   write(irefwr,*)'integrated surface shear stress = ',testtau
endif

end subroutine force_balance_setup
