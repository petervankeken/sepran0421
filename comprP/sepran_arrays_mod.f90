module sepran_arrays
implicit none
integer :: isol1=0,isol2=0,isol3=0
integer :: kmesh1=0,kprob1=0
integer :: intmat1=0,matr1,irhsd1=0,isolold1(3)=0,matrm1=0
integer :: intmat2=0,matr2=0,irhsd2=0,isolold2=0,matrm2=0
integer :: isecinv=0,ivisc=0,iphi=0,iheat=0,ipress=0,ialpha=0,ibulkmod=0,icp=0,iphi2=0
integer :: igradT=0,icompwork=0,icond=0,irho=0,iadia=0,ibigGamma=0,latent_heat=0,iarea=0
integer :: ivisclin=0,iviscplast=0,idens=0,icurl=0,istream=0,ivisdip=0,idivv=0,ipresse=0,icurlv=0,idens2=0
integer,parameter :: NUSERMAX=10000
integer :: iinput(100),iuser_here(NUSERMAX)
real(kind=8),dimension(:),allocatable :: user_here,funcx,funcy,coorint,solint
real(kind=8) :: u1lc(3)
integer :: iu1lc(3)
end module sepran_arrays
