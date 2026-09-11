subroutine getvisdip
use sepmodulekmesh
use sepran_arrays
use coeff
implicit none
integer :: iuserlc(100),ip,jtime
real(kind=8) :: userlc(100),phi(100000)
integer :: iinder(8),icheld

phi=0.0_8
userlc=0.0_8
iuserlc=0
userlc(1)=1.0_8*100
iuserlc(1)=100
iuserlc(2)=1

! First find secinv:secinv
! Use of icheld=11 (viscous dissipation) seems inaccurate. 
! With 0220 version it seems better to compute edot::edot which
! should be equivalent to tau::edot but apparently is not
icheld=10
iinder=0
iinder(1)=8
iinder(2)=1
iinder(3)=0
iinder(4)=icheld
iinder(8)=2

iuserlc(6) = 8
! define integer :: information:    jtime, modelv intrule icoor mcontv
! ( 1) Type of Navier-Stokes equations
jtime = 0
ip=7
iuserlc(ip+1) = jtime
! ( 2) - type of constitutive equations
! specify eta through coefficient 12
iuserlc(ip+2) = 1
! ( 3) - type of interpolation
iuserlc(ip+3) = intrule900 + 100*interpol900
! ( 4) - type of coordinate system (0=Cartesian; 1=Axisymmetric)
iuserlc(ip+4) = icoor900
! ( 5) - compressibility formulation
iuserlc(ip+5) = mcontv

! real information
iuserlc(ip+6)=-7
 userlc(7)=0.0_8 ! eps 
iuserlc(ip+7)=-8
 userlc(8)=1.0_8 ! rho assumed constant here
iuserlc(ip+12)=-9
 userlc(9)=1d0 ! eta assumed constant here
!write(irefwr,*) 'iphi, isol: ',iphi,isol1

call deriv(iinder,iphi,kmesh1,kprob1,isol1,iuserlc,userlc)

!call prinrv(iphi,kmesh1,kprob1,4,1,'visdip')
!call copyvc(isol2,iphi2)
!write(irefwr,*) 'iphi, isol: ',iphi,isol1,iphi2
!call change_problem_for_phi(iphi,iphi2)
!call pecopy(0,phi,iphi)
!write(irefwr,*) 'minval/maxval phi: ',minval(phi(1:npoint)),maxval(phi(1:npoint))

! multiply with viscosity
call coefvis()
call pemultphivisc(icheld) !kmesh1,kprob1,iphi,ivisc)

end subroutine getvisdip
