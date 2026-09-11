program rho_Gamma
implicit none

real(kind=8) :: T0=1500,z0=400.0_8/2885,gamma=0.1_8/1500.0_8,z,T,zi,dz=0.005_8

real(kind=8) :: bigGammabar,bigGammabar0

z=z0
T=T0-200
zi=z-z0-gamma*(T-T0)
bigGammabar0=0.5*(1+tanh((z-z0)/dz))
bigGammabar=0.5*(1+tanh(zi/dz))
write(6,'(8f12.3,:)') z-z0,zi,bigGammabar0,bigGammabar,bigGammabar-bigGammabar0



end program rho_Gamma
