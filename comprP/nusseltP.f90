subroutine nusseltP(q_a)
use sepmoduleoldrouts
use sepmodulecomio  ! provides lu's for reading & writing
use sepmodulecpack  ! parallel computing interface
use sepran_arrays ! access to kmesh, kprob, isol, etc. and interfaces to subroutines
use sepran_interface
use geometry
use convparam
use control
use coeff
implicit none
real(kind=8) :: q_a(*)
integer :: ihelp=0
real(kind=8) :: done=1.0_8,dzero=0.0_8
real(kind=8) :: rcondtop,rcondbot,funccf
save ihelp

if (cyl.and.axi) then
   if (print_node) then
      write(irefwr,*) 'PERROR(nusseltP): not yet suited for cyl & axi: ',cyl,axi
   endif
   call instop
endif

if (compress.and.solve_for_Tperturb) then
   call compute_total_T()
else
   call copyvc(isol2,isol3)
endif

rcondbot=1.0_8
rcondtop=1.0_8
if (cyl) then
   if (icondtype>=1) then
      rcondbot=funccf(4,r1,dzero,dzero)
      rcondtop=funccf(4,r2,dzero,dzero)
   endif
   call deriva(2,2,0,1,2,igradT,kmesh1,kprob1,isol3,isol3,iuser_here,user_here,ihelp)
   ! write(6,*) 'igradT, isol3: ',igradT,isol3
   q_a(1) = bounin(2,3,1,1,kmesh1,kprob1,itop,itop,igradt,iuser_here,user_here)
   q_a(2) = bounin(2,3,1,1,kmesh1,kprob1,ibottom,ibottom,igradt,iuser_here,user_here)
   q_a(3) = bounin(1,3,1,1,kmesh1,kprob1,ibottom,ibottom,isol3,iuser_here,user_here)
   q_a(4) = bounin(1,3,1,1,kmesh1,kprob1,itop,itop,isol3,iuser_here,user_here)
   ! temperature per unit surface area
else
   if (icondtype>=1) then
      rcondbot=funccf(4,dzero,done,dzero)
      rcondtop=funccf(4,dzero,dzero,dzero)
   endif
   call deriva(0,2,0,1,0,igradT,kmesh1,kprob1,isol3,isol3,iuser_here,user_here,ihelp)
   q_a(1)=-bounin(2,3,1,1,kmesh1,kprob1,itop,itop,igradT,iuser_here,user_here)
   q_a(2)=-bounin(2,3,1,1,kmesh1,kprob1,ibottom,ibottom,igradT,iuser_here,user_here)
   q_a(3)=bounin(1,3,1,1,kmesh1,kprob1,ibottom,ibottom,isol3,iuser_here,user_here)
   q_a(4)=bounin(1,3,1,1,kmesh1,kprob1,itop,itop,isol3,iuser_here,user_here)
endif
q_a(1) = -rcondtop*q_a(1)/surftop
q_a(2) = -rcondbot*q_a(2)/surfbot
q_a(3) = q_a(3)/surfbot
q_a(4) = q_a(4)/surftop
q_a(1:4) = abs(q_a(1:4))
q_a(1:4)=abs(q_a(1:4))


end subroutine nusseltP

subroutine compute_total_T()
use sepran_arrays
implicit none
integer :: iinvec(20),i,ipoint,i1
real(kind=8) :: rinvec(20),p,q,tmin,tmax

iinvec(1)=5
! *** With option 27 create linear combination: isol3  = isol2 + iadia
iinvec(2)=27
iinvec(3)=0
iinvec(4)=0
iinvec(5)=1
! Store option for manipulation in rinvec(1) (see funvec)
rinvec(1)=1d0
rinvec(2)=1d0
call manvec(iinvec,rinvec,isol2,iadia,isol3,kmesh1,kprob1)
end subroutine compute_total_T

