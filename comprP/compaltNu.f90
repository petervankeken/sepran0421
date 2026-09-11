subroutine compaltNu(kmesh, kprob, itemp,isol1,iphiy,irho,ipress,icond,dTdzint,icurve,Qz,Nfluxz)
use sepmodulekmesh
use sepmodulekprob
use sepmoduleprobinf
use sepmodulecomio
use control
use tracers ! for rlampix
use geometry
implicit none
integer, intent(in) :: kmesh,kprob,itemp,isol1,iphiy,irho,Nfluxz,icurve,ipress,icond
real(kind=8), intent(out) :: dTdzint,Qz(*)
interface 
   subroutine compaltNu001(ndim,npoint,coor,Tsol,isol1,phi,nusol,TN,wN,phiyN,N,dx,dz,icurve)
   integer, intent(in) :: ndim,npoint,nusol,N,isol1,icurve
   real(kind=8), intent(in) :: Tsol(:),coor(:,:),phi(:)
   real(kind=8), intent(out) :: TN(:),wN(:),phiyN(:),dx,dz
   end subroutine compaltNu001
   subroutine compaltNu003(ndim,npoint,coor,Tsol,u,w,phiy,rho,press,cond,nusol,TN,uN,wN,phiyN,rhoN,pressN,kN,N,irow)
   integer, intent(in) :: ndim,npoint,nusol,N,irow
   real(kind=8), intent(in) :: Tsol(:),coor(:,:),phiy(:),u(:),w(:),rho(:),press(:),cond(:)
   real(kind=8), intent(out) :: TN(:),uN(:),wN(:),phiyN(:),rhoN(:),pressN(:),kN(:)
   end subroutine compaltNu003
end interface
integer :: i,irow,N
real(kind=8), allocatable, dimension(:) :: TN,uN,wN,phiyN,dTdz,u,w,rhoN,pressN,kN
real(kind=8) :: dx,dz
logical :: FD

Qz(1:Nfluxz)=0.0_8
call sepactsolbf1(itemp)
if (nusoli /= npoint) then
   if (print_node) write(irefwr,*) 'PERROR(altNu): nusol /= npoint: ',nusoli,npoint
   call instop
endif
N=nint(sqrt(1.0_8*nusoli))
if (N*N /= nusoli) then
   if (print_node) write(irefwr,*) 'PERROR(compaltNu): N*N /= nusol: ',N*N,nusoli
   call instop
endif

if (HoLiu87.and.Nfluxz==-1) then
   if (print_node) then 
      write(irefwr,*) 'PERROR(compaltNu): incompatible options: HoLiu87 but Nfluxz has not been set'
      write(irefwr,*) '                   Nfluxz=',Nfluxz,' but should be the number of nodal points in vertical'
   endif
   call instop
endif
    

FD=.false.
!write(6,*) 'after compaltNu001: ',FD
if (FD.or..not.HoLiu87) then
   allocate(TN(3*N),wN(3*N),phiyN(3*N),dTdz(N))
   ! Retrieve temperature values and upward velocity component at top and at the two rows below it 
   write(6,*) 'call to compaltNu001: ',FD,.not.HoLiu87,FD.or..not.HoLiu87
   call compaltNu001(ndim,npoint,coor,ks(itemp)%sol,isol1,ks(iphiy)%sol,nusoli,TN,wN,phiyN,N,dx,dz,icurve)
   if (FD) then
     ! first order
     dTdz=TN(1:N)-TN(N+1:2*N)
     dTdz=dTdz/dz
     ! second order central
     dTdz=TN(1:N)-TN(2*N+1:3*N)
     dTdz=0.5*dTdz/dz
     ! second order backward
     dTdz=3*TN(1:N)-4*TN(N+1:2*N)+TN(2*N+1:3*N)
     dTdz=0.5*dTdz/dz
     ! trapezoid rule
     dTdzint=0.5*dTdz(1)+0.5*dTdz(N)+sum(dTdz(2:N-1))
     dTdzint=dTdzint*dx
   else
     ! classic finite element integration of dT/dx across the facets of the elements at the boundary
     call compaltNu002(TN,wN,N,dx,dz,dTdzint)
   endif
   deallocate(TN,wN,phiyN,dTdz)
endif

! N is number of nodal points in vertical
! Nfluxz is the number of bilinear elements spanned by the nodalpoints

if (HoLiu87) then
   if (nint(rlampix*1000.0_8) /= 1000) then
      if (print_node) then 
         write(irefwr,*) 'PERROR(computealtNu): HoLiu87 approach is only implemented for rlam=1 but rlam=',rlampix
      endif
      call instop
   endif
   dz=1.0_8/(Nfluxz-1)
   dx=dz ! clearly we assume here that mesh is uniform with same resolution in x and y...
   !write(6,'(''xcmin,xcmax: '',4e15.7,3i5)') xcmin,xcmax,dx,dz,nx,ny,Nfluxz
   allocate(TN(2*N),wN(2*N),phiyN(2*N),rhoN(2*N),uN(2*N),pressN(2*N),kN(2*N))
   allocate(u(npoint),w(npoint))
   call pecopy(2,u(1:npoint),isol1)
   call pecopy(3,w(1:npoint),isol1)
   !write(6,*) 'min/max w: ',minval(w(1:npoint)),maxval(w(1:npoint))
   !!!!! Assume ipress is filled if (add_P_to_Nu) call coefpress  ! form vector ipress
   if (petest) open(9,file='Qz_components.dat')
   do irow=1,N-1 
      ! step through the grid to compute Q(z) 
      ! extract two rows of grid points representing the corner points of bilinear elements
  !subroutine compaltNu003(ndim,npoint,coor,Tsol,u,w,phiy,rho,press,cond,nusol,TN,uN,wN,phiyN,rhoN,pressN,kN,N,irow)
       call compaltNu003(ndim,npoint,coor,ks(itemp)%sol,u,w,ks(iphiy)%sol,ks(irho)%sol,ks(ipress)%sol,ks(icond)%sol, & 
          & nusoli,TN,uN,wN,phiyN,rhoN,pressN,kN,N,irow)
      ! clearly we assume here that mesh is uniform with same resolution in x and y...
      call compaltNuHoLiu(Tn,uN,wN,phiyN,rhoN,pressN,kN,N,dx,dx,dTdzint,irow)
      Qz(irow)=dTdzint
   enddo
   if (petest) close(9)
   deallocate(TN,wN,uN,phiyN,u,w,rhoN,pressN,kN)
endif
   
end subroutine compaltNu

! Implement flux computation of Ho-Liu et al., 1987
! Construct constant flux in bilinear triangles that combines conductive and advective flux
! TN,wN,phiyN arrays contain three rows with top row first.
! Orientation needs to be taken into account when constructing top or bottom row of the bilinear elements
! PvK October 2021
subroutine compaltNuHoLiu(TN,uN,wN,phiyN,rhoN,pressN,kN,N,dx,dz,dTdzint,irow)
use sepmodulecomio
use control
use geometry
use coeff
implicit none
integer,intent(in) :: N,irow
real(kind=8), dimension(2*N),intent(in) :: TN,uN,wN,phiyN,rhoN,pressN,kN
real(kind=8),intent(in) :: dx,dz ! spacing between nodal points
real(kind=8),intent(out) :: dTdzint
integer :: i,iNe,j,k,ip,Nelem,orientation
integer, parameter :: Nk=4,Ni=4
real(kind=8) :: phiy_part,adv_part,cond_part,Th(Ni),qzh,u(Ni),w(Ni),Phiy1,Phiy2,Phiy3,Phiy4,kh(Ni),kav
real(kind=8) :: rho(4),funccf,yh,tau21,tau22,cond_int,adv_int,phiyint,Pw1,Pw2,Pw3,Pw4,Pint
logical :: first=.true.
real(kind=8) :: Ne(4,4),dNedxi(4,4),dNedeta(4,4),xg(2,4),xi,eta,upwardvelocity,P_part ! Shapefunctions and Guass quadrature points
real(kind=8),parameter :: invsqrt3=1.0_8/sqrt(3.0_8)
save first

   ! set up shapefunction values in Gauss points
   xg(1,1)= -invsqrt3
   xg(2,1)= -invsqrt3
   xg(1,2)=  invsqrt3
   xg(2,2)= -invsqrt3
   xg(1,3)=  invsqrt3
   xg(2,3)=  invsqrt3
   xg(1,4)= -invsqrt3
   xg(2,4)=  invsqrt3
   do k=1,4
      xi=xg(1,k)
      eta=xg(2,k)
      Ne(1,k)=(1.0_8-eta)*(1.0_8-xi)
      Ne(2,k)=(1.0_8+eta)*(1.0_8-xi)
      Ne(3,k)=(1.0_8+eta)*(1.0_8+xi)
      Ne(4,k)=(1.0_8-eta)*(1.0_8+xi)
      dNedxi(1,k)=-1+xi
      dNedxi(2,k)= 1-xi
      dNedxi(3,k)= 1+xi
      dNedxi(4,k)=-1-xi
      dNedeta(1,k)= eta-1
      dNedeta(2,k)=-eta-1
      dNedeta(3,k)= eta+1
      dNedeta(4,k)=-eta+1
   enddo
   Ne=Ne/4
   dNedxi=dNedxi/4
   dNedeta=dNedeta/4
   do k=1,4
      if (abs(sum(Ne(1:4,k))-1.0)>1e-7 ) then 
         if (print_node) write(irefwr,*) 'PERROR(compaltNuHoLiu): problem with shapefunction ',k,' sum /= 1: ',sum(Ne(1:4,k))
         call instop
      endif
   enddo

phiyint=0
adv_int=0.0_8
cond_int=0.0_8
Pint=0.0_8
dTdzint=0.0_8
! Loop over row of elements
!write(6,*) 'TN(x=0,1),Tn/wN(x=0,2): ',TN(1),TN(N+1),wN(N+1)
!if (irow==22) write(6,*) 'N=',N
do i=1,N-1
   Th(1)=TN(i)
   Th(2)=TN(i+1)
   Th(3)=TN(i+N+1)
   Th(4)=TN(i+N)
   kh(1)=kN(i)
   kh(2)=kN(i+1)
   kh(3)=kN(i+N+1)
   kh(4)=kN(i+N)
   kav=sum(kh(1:4))/4
   !if (add_wT_Nu_quad) then
   !   if (abs(Di)>1e-7.and.mcontv>=1) then
   !      ! define rho in quadrature points
   !      do k=1,4
   !         yh=(irow-1)*dz+dz*(1-xg(2,k))/2.0_8
   !         rho(k)=exp(DiG*(1-yh))
   !         !rho(k)=funccf(3,0.0_8,yh,0.0_8)
   !      enddo
   !   else
   !      rho(1:4)=1.0_8
   !   endif
   !else

   ! look up from nodal point info
   rho(1)=rhoN(i)
   rho(2)=rhoN(i+1)
   rho(3)=rhoN(i+N+1)
   rho(4)=rhoN(i+N)
   !endif
   w(1)=wN(i)
   w(2)=wN(i+1)
   w(3)=wN(i+N+1)
   w(4)=wN(i+N)
   u(1)=uN(i)
   u(2)=uN(i+1)
   u(3)=uN(i+N+1)
   u(4)=uN(i+N)
   Phiy1=phiyN(i)
   Phiy2=phiyN(i+1)
   Phiy3=phiyN(i+N+1)
   Phiy4=phiyN(i+N)
   Pw1=pressN(i)*w(1)
   Pw2=pressN(i+1)*w(2)
   Pw3=pressN(i+N+1)*w(3)
   Pw4=pressN(i+N)*w(4)
   
   ! qzh in element
   !qzh=1.0_8/(2*dz)*(kh(1)*Th(1)+kh(1)*Th(2)-kh(3)*Th(3)-kh(4)*Th(4))
   cond_part = kav/(2*dz)*(Th(1)+Th(2)-Th(3)-Th(4))
   qzh = cond_part
   adv_part=0.0_8   
   phiy_part=0.0_8
   !if (add_wT_Nu.and..not.add_wT_Nu_quad) then
      ! add advective part
      adv_part=adv_part+rho(1)*w(1)*(4*Th(1)+2*(Th(2)+Th(4))+Th(3))
      adv_part=adv_part+rho(2)*w(2)*(4*Th(2)+2*(Th(1)+Th(3))+Th(4))
      adv_part=adv_part+rho(3)*w(3)*(4*Th(3)+2*(Th(2)+Th(4))+Th(1))
      adv_part=adv_part+rho(4)*w(4)*(4*Th(4)+2*(Th(1)+Th(3))+Th(2))
      adv_part=adv_part/(4.0_8*9.0_8)
      qzh=qzh+adv_part
   !else if (add_wT_Nu) then
   !   ! use explicit Gaussian quadrature to improve on rho(y)
   !   adv_part=0.0_8
   !   do k=1,Nk   ! loop over Gauss points
   !      do iNe=1,Ni ! loop over shapefunctions for temperature
   !         upwardvelocity=0.0_8
   !         do j=1,Ni ! loop over shapefunctions for velocity
   !            upwardvelocity=upwardvelocity+w(j)*Ne(j,k)
   !         enddo
   !         adv_part=adv_part+rho(k)*Th(iNe)*upwardvelocity*Ne(iNe,k)
   !      enddo
   !   enddo
   !   adv_part=adv_part/4 ! scale for volume of reference element
   !   qzh=qzh+adv_part
   !endif
   !if (add_phi_Nu.and..not.add_phi_Nu_quad) then
      ! add partial viscous dissipation term u*tau_21+v*tau_22
      !phiy_part=phiy_part+4*Phiy1+2*(Phiy2+Phiy4)+Phiy3
      !phiy_part=phiy_part+4*Phiy2+2*(Phiy1+Phiy3)+Phiy4
      !phiy_part=phiy_part+4*Phiy3+2*(Phiy2+Phiy4)+Phiy1
      !phiy_part=phiy_part+4*Phiy4+2*(Phiy1+Phiy3)+Phiy2
      !phiy_part=phiy_part/(4.0_8*9.0_8)
      phiy_part=1.0_8/4*(Phiy1+Phiy2+Phiy3+Phiy4)
   !else if (add_phi_Nu) then
   ! numerical quadrature of viscous dissipation term v_j tau_2j
   !  phiy_part=0.0_8
   !  do k=1,4
   !     do iNe=1,4
   !        ! in each Gauss point form tau_21 and tau_22 from velocity gradients
   !        tau21=0.0_8
   !        tau22=0.0_8
   !        do j=1,4
   !           ! dv1/dx2+dv2/dx1
   !           tau21 = tau21 + u(j)*dNedxi(j,k)+w(j)*dNedeta(j,k)
   !           if (mcontv>=1) then
   !              ! compressible stress tensor - use gradients
   !              tau22 = tau22 + 4.0_8/3*w(j)*dNedxi(j,k) - 2.0_8/3*u(j)*dNedeta(j,k)
   !           else
   !              tau22 = tau22 + 2.0_8*u(j)*dNedeta(j,k) + 2.0_8*w(j)*dNedxi(j,k)
   !           endif
   !        enddo ! j
   !        ! now form v_j tau_2j = u*tau_21+v*tau_22 = sum u(i)*Ne(i)*tau21 + v(i)*Ne(i)*tau22
   !        phiy_part = phiy_part + u(iNe)*Ne(iNe,k)*tau21 + w(iNe)*Ne(iNe,k)*tau22
   !     enddo ! iNe
   !  enddo ! k
   !  phiy_part=phiy_part/4 ! rescale for volume of the reference element
   !  phiy_part=phiy_part/(2*dz) ! rescale for chain rule (dx==dz, of course)
   !endif
   if (add_P_to_Nu) then
      P_part=1.0_8/4*(Pw1+Pw2+Pw3+Pw4)
   endif
   ! Note: negative sign
   qzh=qzh-DiRa*phiy_part+DiRa*P_part
   if (petest.and.irow==22) write(6,'(i5,f12.7,5e12.5)') i,qzh,P_part,Pw1,Pw2,Pw3,Pw4
   cond_int=cond_int+cond_part*dx
   adv_int=adv_int+adv_part*dx
   phiyint=phiyint+phiy_part*dx
   Pint=Pint+P_part*dx
   dTdzint=dTdzint+qzh*dx
enddo
if (petest) write(9,'(8e15.7)') (irow-0.5)*dx,cond_int,adv_int,rho(1),-DiRa*phiyint,DiRa*Pint,Pint,dTdzint

!if (irow==22) then
!   write(6,*)
!   write(6,'(''row 22: '',2i5,8e15.7 )') irow,N,adv_part,cond_int,adv_int,wN(1),wn(2),Th(1),Th(2),dx
!endif

end subroutine compaltNuHoLiu

! Find temperature, velocity, and viscous dissipation in two rows of nodal points 
! Store row info in TN,wN,phiyN.
! PvK October 2021
!      call compaltNu003(ndim,npoint,coor,ks(itemp)%sol,u,w,ks(iphiy)%sol,ks(irho)%sol,ks(ipress)%sol,ks(icond)%sol, & 
!         & nusoli,TN,uN,wN,phiyN,rhoN,pressN,kN,N,irow)
subroutine compaltNu003(ndim,npoint,coor,Tsol,          u,w,phiy,         rho,         press,         cond,  & 
          & nusol, TN,uN,wN,phiyN,rhoN,pressN,kN,N,irow)
use sepmodulecomio
use geometry
implicit none
integer, intent(in) :: ndim,npoint,nusol,N,irow
real(kind=8), intent(in) :: Tsol(:),coor(:,:),phiy(:),u(:),w(:),rho(:),press(:),cond(:)
real(kind=8), intent(out) :: TN(:),uN(:),wN(:),phiyN(:),rhoN(:),pressN(:),kN(:)
integer :: i,j

! quick solution to find upward velocity component in Cartesian
!if (irow==1) then
!   write(6,*) 'minmaxval T: ',minval(Tsol(1:npoint)),maxval(Tsol(1:npoint))
!   write(6,*) 'minmaxval k: ',minval(cond(1:npoint)),maxval(cond(1:npoint))
!   write(6,*) 'minmaxval r: ',minval(rho(1:npoint)),maxval(rho(1:npoint))
!   write(6,*) 'minmaxval P: ',minval(press(1:npoint)),maxval(press(1:npoint))
!   write(6,*) 'minmaxval w: ',minval(w(1:npoint)),maxval(w(1:npoint))
!   write(6,*) 'minmaxval f: ',minval(phiy(1:npoint)),maxval(phiy(1:npoint))
!endif

TN=0.0_8
wN=0.0_8
uN=0.0_8
phiyN=0.0_8
rhoN=0.0_8
pressN=0.0_8
kN=0.0_8
! simpler copy of bottom rows
do i=1,2*N  ! copies row irow and irow+1
   TN(i)=Tsol(i+(irow-1)*N)
   wN(i)=w(i+(irow-1)*N)
   phiyN(i)=phiy(i+(irow-1)*N)
   rhoN(i)=rho(i+(irow-1)*N)
   pressN(i)=press(i+(irow-1)*N)
   kN(i)=cond(i+(irow-1)*N)
enddo

end subroutine compaltNu003



























! Find temperature, velocity, and viscous dissipation in the top and bottom rows for alternative Nu computation
! Store row info in TN,wN,phiyN.
! PvK October 2021
subroutine compaltNu001(ndim,npoint,coor,Tsol,isol1,phi,nusol,TN,wN,phiyN,N,dx,dz,icurve)
use sepmodulecomio
use geometry
use control
implicit none
integer, intent(in) :: ndim,npoint,nusol,N,isol1,icurve
real(kind=8), intent(in) :: Tsol(:),coor(:,:),phi(:)
real(kind=8), intent(out) :: TN(:),wN(:),phiyN(:),dx,dz
real(kind=8), allocatable, dimension(:) :: w
integer :: i,j
real(kind=8) :: dz2

if (print_node) then
   write(irefwr,*) 'PERROR(compaltNu001): this part of altNu is obsolete until improved'
endif
call instop

! quick solution to find upward velocity component in Cartesian
allocate(w(npoint))
call pecopy(3,w(1:npoint),isol1)
write(6,*) 'minmax w001: ',minval(w(1:npoint)),maxval(w(1:npoint))

wN=0.0_8
phiyN=0.0_8
if (icurve==itop) then
   j=0
   do i=npoint-N+1,npoint ! ,N-1
      j=j+1
      TN(j)=Tsol(i)
      wN(j)=w(i)
   enddo 
   do i=npoint-2*N+1,npoint-N ! ,N-1
      j=j+1
      TN(j)=Tsol(i)
      wN(j)=w(i)
   enddo 
   do i=npoint-3*N+1,npoint-2*N ! ,N-1
      j=j+1
      TN(j)=Tsol(i)
      wN(j)=w(i)
   enddo 
else if (icurve==ibottom) then
   ! simpler copy of bottom rows
   do i=1,3*N
      TN(i)=Tsol(i)
      wN(i)=w(i)
   enddo
endif
!write(6,'(''wN,TN: '',2f12.5)') sum(TN(1:N-1)),sum(wN(1:N-1))
!write(6,'(''wN,TN: '',2f12.5)') sum(TN(N:2*N-1)),sum(wN(N:2*N-1))
!write(6,'(''wN,TN: '',2f12.5)') sum(TN(2*N:3*N-1)),sum(wN(2*N:3*N-1))
deallocate(w)
dz=coor(2,npoint)-coor(2,npoint-N)
dz2=coor(2,npoint-N)-coor(2,npoint-2*N)
dx=1.0_8/(N-1)
!TN=0
end subroutine compaltNu001

! Compute Nu for the top row (hack for comparison with Scott King October 2021)
! Compute effect of dT/dy+wT in the top row on the nodal points of the top row 
subroutine compaltNu002(TN,wN,N,dx,dz,dTdzint)
use sepmodulecomio
use control
implicit none
integer,intent(in) :: N
real(kind=8), dimension(3*N),intent(in) :: TN,wN
real(kind=8),intent(in) :: dx,dz ! spacing between nodal points
real(kind=8),intent(out) :: dTdzint
integer :: i,j,ip,Nelem
real(kind=8) :: xn(6),yn(6),xl(3),yl(3),a(3),b(3),c(3),area_element,adv_part,cond_part
real(kind=8) :: phil(3),phiq(6),dphidksi(6),dphideta(6),Telem(6),Th(3),dThdy(3),dTelemdy(6),wh(6),welem(6)
! reference triangle in barycentric coordinates
xn(1:6)=(/ 0.0_8,0.5_8,1.0_8,0.5_8,0.0_8,0.0_8/)
yn(1:6)=(/ 0.0_8,0.0_8,0.0_8,0.5_8,1.0_8,0.5_8/)
xl(1:3)=xn(1:5:2)
yl(1:3)=yn(1:5:2)

!write(irefwr,*) 'in compaltNu002'
Nelem=(N-1)/2
dTdzint=0.0_8
! loop over quadratic element facets
do i=1,N-1,2 ! ,N-1,2
   ! fill temperature ; reorder 
   Telem(1)=TN(i)
   Telem(2)=TN(i+1)
   Telem(3)=TN(i+2)
   Telem(4)=TN(i+2+N) 
   Telem(5)=TN(i+2+2*N)
   Telem(6)=TN(i+1+N)
   welem(1)=wN(i)
   welem(2)=wN(i+1)
   welem(3)=wN(i+2)
   welem(4)=wN(i+2+N) 
   welem(5)=wN(i+2+2*N)
   welem(6)=wN(i+1+N)
   !write(6,'(''sumT, sumwT: '',2e15.7)') sum(Telem(1:6)),sum(welem(1:6)*Telem(1:6))
   !write(6,'(''Telem: '',36x,36x,6f12.5)') (Telem(j),j=1,6)
   dThdy=0
   wh=0
   ! find gradient in top row of nodal points (defined by eta=0)
   do ip=1,3
      call getshape_quad_xi_eta(xn(ip),yn(ip),phil,phiq,dphidksi,dphideta)
      !write(6,'(2f8.3,3f12.5)') xn,yn,(dphideta(j),j=1,3)
      do j=1,6
         dThdy(ip)=dThdy(ip)+Telem(j)*dphideta(j)
      enddo
      dThdy(ip)=dThdy(ip)/(2*dz)
      cond_part=dThdy(ip)
      ! PvK 10/14/2021: this doesn't work as phi(1:3) can be 1, but wT(1:3)=0 by definition....
      !if (add_wT_Nu) then
      !   ! add advective component; just lump into dThdy
      !   adv_part=0
      !   do j=1,6
      !      adv_part=adv_part +welem(j)*Telem(j)*phiq(j)
      !      write(6,'(10x,5e15.7)') welem(j),Telem(j),phiq(j),welem(j)*Telem(j)*phiq(j)
      !   enddo
      !   dThdy(ip)=dThdy(ip)+adv_part
      !   write(6,'(2i5,5e12.5)') i,ip,cond_part,adv_part,sum(phiq(1:6)),sum(welem(1:6)),sum(Telem(1:6))
      !endif
   enddo
   ! integrate over top element facet. Element width is 2*dx
   dTdzint=dTdzint+dx*(0.5*(dThdy(1)+dThdy(3))+dThdy(2))
enddo
! Correct for height of element (assumed to have reference value of 1 above)
!write(6,*) 'dTdzint=',dTdzint
!call instop
end subroutine compaltNu002

