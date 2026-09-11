program reference_rho
implicit none
real(kind=8), parameter :: Di=0.50,rho_r=3330,delta=2.0_8,gamma_0=1.0
real(kind=8) :: wdim(2),zdim,zdim_ph(2)
real(kind=8) :: z,drhodz,rho2(2886),dz,deltarhorel(2),rho2_1(2885),rho2_2(2885),rho3(2885),deltarho(2)
real(kind=8) :: bigGamma,w,dgdz(2,2886),bigG
integer :: i,N,nph=1,iph

deltarhorel(1)=0.05
deltarhorel(2)=0.1
deltarho(1)=3300*deltarhorel(1) ! 3619*deltarhorel(1)
deltarho(2)=4000*deltarhorel(2) ! 3998
zdim_ph(1)=400
zdim_ph(2)=670
wdim=5
rho2(1)=rho_r
N=2885
dz=1.0/N
! first construct background density
do i=2,N
   z=(i-1)*dz
   drhodz=rho2(i-1)*Di/gamma_0
   rho2(i)=rho2(i-1)+drhodz*dz
enddo
rho2_1(1)=rho_r
do i=2,2885
   z=(i-1)*dz
   zdim=z*2885
   w=wdim(1)/2885
   bigG=bigGamma(zdim,zdim_ph(1),wdim(1))
   dgdz(1,i)=2*deltarho(1)/w*bigG*(1-bigG)
   drhodz=rho2_1(i-1)*Di/gamma_0 + dgdz(1,i)
   rho2_1(i)=rho2_1(i-1)+drhodz*dz
enddo

rho2_2(1)=rho_r
do i=2,2885
   z=(i-1)*dz
   zdim=z*2885
   do iph=1,2
      w=wdim(iph)/2885
      bigG=bigGamma(zdim,zdim_ph(iph),wdim(iph))
      dgdz(iph,i)=2*deltarho(iph)/w*bigG*(1-bigG)
   enddo
   drhodz=rho2_2(i-1)*Di/gamma_0 +  dgdz(1,i) + dgdz(2,i)
   rho2_2(i)=rho2_2(i-1)+drhodz*dz
enddo

do i=1,2885
   z=(i-1)*dz
   zdim=z*2885
   write(6,'(6f12.3)') zdim,rho2(i),rho2_1(i),rho2_2(i),dgdz(1,i),dgdz(2,i)
enddo

end program reference_rho

real(kind=8) function bigGamma(z,z0,w)
implicit none
real(kind=8) :: z,z0,w
bigGamma = 0.5*(1+tanh((z-z0)/w))
end function bigGamma
