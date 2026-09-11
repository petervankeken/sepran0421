program trymc
use sepmoduleoldrouts
use control
implicit none
integer :: kmesh=0,kprob=0,intmat=0,isol=0,irhsd=0,isolold=0,matrm=0,matrs=0
integer, parameter :: NIUSER=1000,NUSER=1000
integer :: iuser(NIUSER)
real(kind=8) :: user(NUSER),resmem,resmem_old
integer :: commat_in(14),iinbld(10),i


! set up system command to check resident memory of the main process
call get_checkmem_command

iuser(1)=NIUSER
user(1)=NUSER*1.0_8
call sepstr(kmesh,kprob,intmat)
commat_in=0
commat_in(1)=5  ! number of entries in commat_in
commat_in(2)=2 ! 
commat_in(3)=0
commat_in(4)=0
commat_in(5)=1  ! iprob = problem number 
call matstruc(commat_in,kmesh,kprob,intmat)

iuser(2:NIUSER)=0
user(2:NUSER)=0.0_8
call filcof(iuser,user,kprob,kmesh,1)
iinbld=0
iinbld(1)=4
iinbld(2) = 1
resmem=get_resmem()
do i=1,20
   resmem_old=resmem
   resmem=get_resmem()
   write(6,'(i7,f12.3)') i,1024*1024*(resmem-resmem_old)
   call build(iinbld,matrs,intmat,kmesh,kprob,irhsd,matrm,isol,isolold,iuser,user)
enddo


end program trymc

