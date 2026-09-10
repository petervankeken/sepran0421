! PvK July 2024
module sepran_arrays
implicit none
integer :: kmesh=0,kprob=0,isol(3)=0,intmat(2)=0,istream=0 
integer :: matr(2)=0,matrm(2)=0,irhsd(2)=0,islold(2),igradt=0
integer :: ivisc=0,idens=0,iadia=0,ivisclin=0,iviscplast=0
integer :: iheat=0,isecinv=0,iphi=0,iplasticity=0,iphi2=0,iwork=0,iwork2=0
integer :: ivisdip=0,ivisdip2=0,icompwork=0,icompwork2=0,icond=0
integer :: ialpha=0,ipress=0,irho=0
integer, parameter :: NUSERMAX=10000
integer :: iuser_here(NUSERMAX)
real(kind=8), dimension(:), allocatable :: user_here
end module sepran_arrays
