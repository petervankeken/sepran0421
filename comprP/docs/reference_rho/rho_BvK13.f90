program reference_rho
implicit none
real(kind=8) :: Di=0.65,rho_r=3330,delta=2.2_8,gamma_0=1.0
real(kind=8) :: wdim(2),zdim,zdim_ph(2)
real(kind=8) :: z,drhodz,rho2(2886),dz,deltarhorel(2),rho2_1(2885),rho2_2(2885),rho3(2885),deltarho(2)
real(kind=8) :: bigGamma,w,dgdz(2,2886),bigG,alpha
integer :: i,N,nph=1,iph

read(5,*) Di,delta
deltarhorel(1)=0.05
deltarhorel(2)=0.1
deltarho(1)=3300*deltarhorel(1) ! 3619*deltarhorel(1)
deltarho(2)=4000*deltarhorel(2) ! 3998*deltarhorel(2)
zdim_ph(1)=400
zdim_ph(2)=670
wdim=5
N=2885
dz=1.0/N
! first construct background density
rho2(1)=rho_r
do i=2,N
   z=(i-1)*dz
   drhodz=Di*alpha(rho2(i-1)/rho_r,delta)*rho2(i-1)
   rho2(i)=rho2(i-1)+drhodz*dz
   !write(6,*) alpha(rho2(i-1)/rho_r,delta),rho2(i),drhodz
enddo
rho2_1(1)=rho_r
do i=2,2885
   z=(i-1)*dz
   zdim=z*2885
   w=wdim(1)/2885
   bigG=bigGamma(zdim,zdim_ph(1),wdim(1))
   dgdz(1,i)=2*deltarho(1)/w*bigG*(1-bigG)
   drhodz=Di*alpha(rho2_1(i-1)/rho_r,delta)*rho2_1(i-1) + &
      & dgdz(1,i)
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
   drhodz=Di*alpha(rho2_2(i-1)/rho_r,delta)*rho2_2(i-1) +&
      & dgdz(1,i) + dgdz(2,i)
   rho2_2(i)=rho2_2(i-1)+drhodz*dz
enddo

do i=1,2885
   z=(i-1)*dz
   zdim=z*2885
   write(6,'(7f12.3)') zdim,rho2(i),rho2_1(i),rho2_2(i),dgdz(1,i),dgdz(2,i),alpha(rho2_2(i)/rho_r,delta)
enddo

end program reference_rho

real(kind=8) function bigGamma(z,z0,w)
implicit none
real(kind=8) :: z,z0,w
bigGamma = 0.5*(1+tanh((z-z0)/w))
end function bigGamma

real(kind=8) function alpha(rho,delta)
implicit none
real(kind=8) :: rho,delta
alpha=rho**(-delta)
end function alpha
