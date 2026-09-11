module eos
contains

  real(kind=8) function get_rhobar_dim(zdim)
  use convparam
  use coeff
  implicit none
  real(kind=8) :: zdim,z,rho_n

  z=zdim/height_dim
  if (iadiabat<=0.or.EBA) then
     get_rhobar_dim=rho_dim
  else if (iadiabat==2 .or. iadiabat==1) then
     ! WA with constant alpha
     get_rhobar_dim=rho_dim*exp(DiG*z)
  else if (iadiabat==3) then
     ! WA with variable alpha following BvK13
     get_rhobar_dim=rho_dim*(-DiG*palpha*z+1)**(-1.0_8/palpha)
  endif
  !write(irefwr,'(i5,8e15.7)') iadiabat,z,zdim,get_rhobar_dim,DiG,rho_dim,-DiG*palpha*z+1
  
  ! catch special test case with simplified rhobar
  if (fake_rho_bar) then
    if (step_rho_background) then
      rho_n=1.0_8
      if (z>0.5_8) rho_n=1.0_8+drho_background_dense
    else  if (exp_rho_background) then
      rho_n=exp(0.5*z)
    else
      rho_n=1+drho_background_dense*z
    endif
    get_rhobar_dim=rho_dim*rho_n
  endif

  
  end function get_rhobar_dim
   
  ! Specify adiabat in K
  real(kind=8) function get_Tbar_dim(zdim)
  use convparam
  use coeff
  use control
  implicit none
  real(kind=8) :: zdim,y

  y=1.0_8-zdim/height_dim
  if (EBA) then
     get_Tbar_dim  = 0d0
  else 
     get_Tbar_dim = Ts_dimK
  endif

  end function get_Tbar_dim 

  ! get_alpha_dim: specify non-dimensional value for thermal expansivity
  ! function will become obsolete since we'll move to lookup from vector ialpha
  real(kind=8) function get_alpha_dim(zdim,rhodim)
  use convparam
  implicit none
  real(kind=8) :: zdim,rho_d,rhodim
  real(kind=8) :: y,rho_nd,z

! y=1-zdim/height_dim
! z=zdim/height_dim
! rho_nd=rhodim/rho_dim ! nice names
! get_alpha_dim=alpha_dim/(rho_nd*rho_nd)
! return
  
  rho_d=get_rhobar_dim(zdim)
  rho_nd=rho_d/rho_dim
  if (ialphatype==0) then
     get_alpha_dim=1
  else if (ialphatype==3) then
     get_alpha_dim=0.2+0.8*y
  else if (ialphatype==1) then
      get_alpha_dim=1.0_8/(rho_nd*rho_nd)
  endif
  get_alpha_dim=get_alpha_dim*alpha_dim

  end function get_alpha_dim

  real(kind=8) function d2adiabatdz2(y)
  use convparam
  use coeff
  implicit none
  real(kind=8) :: y

  if (iadiabat.eq.1.or.EBA) then
     d2adiabatdz2 = 0d0
  else if (iadiabat.eq.2) then
     ! Adiabat is assumed to be Tbar=Ta*(exp(Di*z)-1)
     d2adiabatdz2 = Tbars_nondimK*Di*Di*exp(Di*(1-y))! / deltaT_dim
  endif

  end function d2adiabatdz2

end module eos
