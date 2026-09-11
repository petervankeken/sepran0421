subroutine tracrhsd01_903_cart(kmesh1,itrac_start,itrac_end)
use sepmodulecomio
use sepmodulekmesh
use coeff
use msper01
use mpetrac
use sepran_interface
implicit none
integer :: kmesh1,itrac_start,itrac_end
integer :: itrac,ielh
integer :: nodno(7),nodlin(3),isub,i
real(kind=8) :: rl(3),phiq(7),shapef(7),xi,eta,xn(6),yn(6)
real(kind=8) :: xm,ym,weight,funccf,z,rho_here
logical :: correct

rho_here=1.0
if (ibuoy_trac==1) then
  do itrac=itrac_start,itrac_end
     xm = tracer(itrac)%x
     ym = tracer(itrac)%y
     if (compress.and..not.stretch_tracers) then
        rho_here=funccf(3,xm,ym,z)
     endif
     call pedetel(1,xm,ym,ielh)
     call sper01(kmeshc,coor,nodno,xn,yn,ielh)
     correct=checkinelem(3,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,ielh)
     weight = rho_here*tracer(itrac)%density*funccf(7,xm,ym,0d0)
     f2_tracrhsd(1:6,ielh,i_f_tracrhsd) = f2_tracrhsd(1:6,ielh,i_f_tracrhsd) +weight*phiq(1:6)
     !write(irefwr,'(''f2: '',16e12.3)') xm,ym,rho_here,funccf(7,xm,ym,0d0),phiq(1:6),f2_tracrhsd(1:6,ielh,i_f_tracrhsd)
  enddo
else 
  do itrac=itrac_start,itrac_end
     xm = tracer(itrac)%xnew
     ym = tracer(itrac)%ynew
     call pedetel(1,xm,ym,ielh)
     call sper01(kmeshc,coor,nodno,xn,yn,ielh)
     correct=checkinelem(3,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,ielh)
     if (compress) then
        rho_here=funccf(3,xm,ym,z)
     endif
     weight = rho_here*tracer(itrac)%density*funccf(7,xm,ym,0d0)
     f2_tracrhsd(1:6,ielh,i_f_tracrhsd) = f2_tracrhsd(1:6,ielh,i_f_tracrhsd) +weight*phiq(1:6)
     !write(irefwr,'(''f2: '',16e12.3)') xm,ym,rho_here,funccf(7,xm,ym,0d0),phiq(1:6),f2_tracrhsd(1:6,ielh,i_f_tracrhsd)
  enddo
endif
!itrac=itrac-1
!write(irefwr,'(L2,3i6,15e12.3)') correct,ibuoy_trac,ielh,itrac,xm,ym,f2_tracrhsd(1:6,ielh,i_f_tracrhsd),phiq(1:6),tracer(itrac)%density

end subroutine tracrhsd01_903_cart

subroutine tracrhsd01_900_cart(kmesh1,itrac_start,itrac_end)
use sepmodulekmesh
use coeff
use msper01
use mpetrac
use sepran_interface
implicit none
integer :: kmesh1,itrac_start,itrac_end
integer :: itrac,ielh
integer :: nodno(7),nodlin(3),isub,i
real(kind=8) :: rl(3),phiq(7),shapef(7),xi,eta,xn(6),yn(6)
real(kind=8) :: xm,ym,weight,funccf,z,rho_here
logical :: correct

rho_here=1.0
if (ibuoy_trac==1) then
  do itrac=itrac_start,itrac_end
     xm = tracer(itrac)%x
     ym = tracer(itrac)%y
     if (compress.and..not.stretch_tracers) then
        rho_here=funccf(3,xm,ym,z)
     endif
     call pedetel(1,xm,ym,ielh)
     call sper01(kmeshc,coor,nodno,xn,yn,ielh)
     correct=checkinelem(3,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,ielh)
     call getshape7_xi_eta(xi,eta,shapef)
     
     weight = rho_here*tracer(itrac)%density*funccf(7,xm,ym,0d0)
     f2_tracrhsd(1:7,ielh,i_f_tracrhsd) = f2_tracrhsd(1:7,ielh,i_f_tracrhsd) +weight*shapef(1:7)
     !write(irefwr,'(''f2_1: '',16e12.3)') xm,ym,rho_here,funccf(7,xm,ym,0d0),phiq(1:6),f2_tracrhsd(1:6,ielh,i_f_tracrhsd)
  enddo
else 
  do itrac=itrac_start,itrac_end
     xm = tracer(itrac)%xnew
     ym = tracer(itrac)%ynew
     call pedetel(1,xm,ym,ielh)
     call sper01(kmeshc,coor,nodno,xn,yn,ielh)
     correct=checkinelem(3,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,ielh)
     if (compress.and..not.stretch_tracers) then
        rho_here=funccf(3,xm,ym,z)
     endif
     call getshape7_xi_eta(xi,eta,shapef)
     weight = rho_here*tracer(itrac)%density*funccf(7,xm,ym,0d0)
     f2_tracrhsd(1:7,ielh,i_f_tracrhsd) = f2_tracrhsd(1:7,ielh,i_f_tracrhsd) +weight*shapef(1:7)
     !write(irefwr,'(''f2_2: '',16e12.3)') xm,ym,rho_here,funccf(7,xm,ym,0d0),phiq(1:6),f2_tracrhsd(1:6,ielh,i_f_tracrhsd)
  enddo
endif

end subroutine tracrhsd01_900_cart


subroutine tracrhsd01_900_cyl(kmesh1,itrac_start,itrac_end)
use coeff
use msper01
use geometry
use mpetrac
use eos
use convparam
use sepran_interface
implicit none
integer :: kmesh1,itrac_start,itrac_end

integer :: itrac,ielh,ioffset
integer :: nodno(7),nodlin(3),isub,i,ikelmo
integer :: imissed,ifound
real(kind=8) :: rl(3),phiq(7),shapef(7),xi,eta,xn(7),yn(7),zdim
real(kind=8) :: xm,ym,weight,funccf,un(6),vn(6),r,rho_here,z
logical :: correct

ioffset=0

rho_here=1.0
if (ibuoy_trac==1) then
   do itrac=itrac_start,itrac_end
      xm = tracer(itrac)%x
      ym = tracer(itrac)%y
      call pedeteltrac(1,ioffset,xm,ym,nodno,nodlin,rl,xn,yn,un,vn,ielh,imissed,ifound, &
          &  shapef,xi,eta,.false.,'stokes',itrac)
      xn(7)=(xn(1)+xn(3)+xn(5))/3
      yn(7)=(yn(1)+yn(3)+yn(5))/3
      call getshape7_xi_eta(xi,eta,shapef)
      weight = tracer(itrac)%density*funccf(7,xm,ym,0d0)
      do i=1,7
         r=sqrt(xn(i)*xn(i)+yn(i)*yn(i))
         if (compress) then
            z=r2-r
            zdim=z*height_dim
            rho_here=get_rhobar_dim(zdim)/rho_dim
         endif
        
         f1_tracrhsd(i,ielh,i_f_tracrhsd) = f1_tracrhsd(i,ielh,i_f_tracrhsd) + rho_here*weight*shapef(i)*xn(i)/r
         f2_tracrhsd(i,ielh,i_f_tracrhsd) = f2_tracrhsd(i,ielh,i_f_tracrhsd) + rho_here*weight*shapef(i)*yn(i)/r
      enddo
   enddo
else
   do itrac=itrac_start,itrac_end
      xm = tracer(itrac)%xnew
      ym = tracer(itrac)%ynew
      call pedeteltrac(1,ioffset,xm,ym,nodno,nodlin,rl,xn,yn,un,vn,ielh,imissed,ifound, &
          &  shapef,xi,eta,.false.,'stokes',itrac)
      xn(7)=(xn(1)+xn(3)+xn(5))/3
      yn(7)=(yn(1)+yn(3)+yn(5))/3
      call getshape7_xi_eta(xi,eta,shapef)
      weight = tracer(itrac)%density*funccf(7,xm,ym,0d0)
      do i=1,7
         r=sqrt(xn(i)*xn(i)+yn(i)*yn(i))
         f1_tracrhsd(i,ielh,i_f_tracrhsd) = f1_tracrhsd(i,ielh,i_f_tracrhsd) + weight*shapef(i)*xn(i)/r
         f2_tracrhsd(i,ielh,i_f_tracrhsd) = f2_tracrhsd(i,ielh,i_f_tracrhsd) + weight*shapef(i)*yn(i)/r
      enddo
   enddo
endif
   
end subroutine tracrhsd01_900_cyl
