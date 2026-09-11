!   INIT_EOS
!
!   Read EoS from file and store in Sepran vector ivec
!   ichoice  = 1    density
!              2    alpha
!              3    cp
!              4    Ks
!              5    conductivity (derived from rho)
!              6    adiabatic temperature
!
!   1D profiles of rho, alpha etc. are stored in pecof900.inc
!   eos_alpha_d etc. contains dimensional values
!   eos_alpha etc. contains non-dimensional values
subroutine init_eos(ichoice,ivec,kmesh1,kprob1)
use sepmoduleoldrouts
use convparam
use coeff
use geometry
use control
!use sepmodulecpack
!use sepmodulemain
use sepmodulekmesh
!use sepmodulemeshinf
!use sepmodulemesh
use sepmodulekprob
!use sepmoduleprobinf
implicit none
interface
    subroutine init_eos01(ichoice,coor,npoint,ndim,vec)
      integer,intent(in) :: ichoice,npoint,ndim
      real(kind=8),intent(in) :: coor(ndim,npoint)
      real(kind=8),intent(inout) :: vec(:)
    end subroutine init_eos01
    subroutine init_eos03(vec,pcond,eos_rhoav,npoint)
      integer,intent(in) :: npoint
      real(kind=8),intent(in) :: pcond,eos_rhoav
      real(kind=8),intent(inout) :: vec(:)
    end subroutine init_eos03
end interface
integer :: ichoice,ivec,kmesh1,kprob1
real(kind=8) :: vec_av,u1(3),y0
integer :: iu1(3),i
integer :: ipvec,ipcoor,ncntln,iuser(100),ihelp
real(kind=8) :: contln(6),format,user(100),volint,p,q,tmin,tmax,adiabat
integer :: i1,ipoint

if (ichoice==5) then
   if (print_node) then
      write(irefwr,*) 'PERROR(init_eos): requires update for conductivity in look up tables for ichoice=',ichoice
   endif
   call instop
endif

iuser(1)=100
 user(1)=100
if (ivec == 0 ) then
   iu1(1)=0
    u1(1)=0d0
   call creavc(0,2,1,ivec,kmesh1,kprob1,iu1,u1,iu1,u1)         
endif
!if (ivec /= 2) then
!    write(irefwr,*) 'PERROR(init_eos): ivec /= 2'
!    call instop
!endif

if (ichoice < 0 .or. ichoice > 6) then
   if (print_node) then
      write(irefwr,*) 'PERROR(init_eos): incompatible ichoice'
      write(irefwr,*) 'ichoice = ',ichoice
   endif
   call instop
endif
 
call sepactsolbf1(ivec)
call init_eos01(ichoice,coor,npoint,ndim,ks(ivec)%sol)

do i=2,iuser(1)
   iuser(i)=0
enddo
do i=2,nint(user(1))
   user(i)=0d0
enddo
iuser(2)=1
iuser(6)=7
iuser(7)=0
iuser(8)=icoor900
iuser(9)=0
iuser(10)=-6
 user(6)=1d0
! vec_av is dimensional average
!!!!!PvK: this needs work for confusion in integration (same as iphi->iphi2)
vec_av = volint(0,2,1,kmesh1,kprob1,ivec,iuser,user,ihelp)
vec_av = vec_av / volume
! write(irefwr,*) 'volume, npoint: ',volume,npoint,vec_av
if (ichoice == 1) then
   eos_rhoav = vec_av
   if (eos_rhoav<1e-3) then
     if (print_node) write(irefwr,*) 'PERROR(init_eos): eos_rhoav is close to 0: ',eos_rhoav
     call instop
   endif
   if (print_node) write(irefwr,*) 'init_eos: eos_rhoav = ',vec_av
else if (ichoice == 2) then
   eos_alphaav = vec_av
else if (ichoice == 3) then
   eos_cpav = vec_av
else if (ichoice == 4) then
   eos_kav = vec_av
else if (ichoice == 6) then
   eos_Ta_av = vec_av
endif

if (eos_type /= 0 .and. ichoice == 5) then
!  special action for conductivity with EoS read from tables. 
!  Compute (rho/rhoav)**pcond
!  write(irefwr,*) 'eos_rhoav: ',eos_rhoav
   call init_eos03(ks(ivec)%sol,pcond,eos_rhoav,npoint)
endif
   
if (ichoice == 5) then
   vec_av = volint(0,2,1,kmesh1,kprob1,ivec,iuser,user,ihelp)
   vec_av = vec_av / volume
   eos_condav = vec_av
endif

if (ichoice == 6 .and. print_node) then
!  verify that adiabat is set up as T1_bar=Tbar-Tbar(0)
   call algebr(6,3,ivec,i1,i1,kmesh1,kprob1,tmin,tmax,p,q,ipoint)
   write(irefwr,'(''init_eos: adiabat tmin/tmax: '',4e15.3)') tmin,tmax,Di
endif
   

end subroutine init_eos

!   INIT_EOS03
!   Compute conductivity = (rho/rho_av)**pcond
subroutine init_eos03(vec,pcond,vecav,npoint)
implicit none
integer,intent(in) :: npoint
real(kind=8),intent(in) :: pcond,vecav
real(kind=8),intent(inout) :: vec(:)
integer :: i

do i=1,npoint
   vec(i) = (vec(i)/vecav)**pcond
enddo

end

!   INIT_EOS01
!   
!   Find values for array vec from lookup tables (stored in pecof900.inc)
!   if eos_type /=0, otherwise use functional description from funccf()
!   Here we read the non-dimensional values prepared in compr_start
!   ichoice = 1   rho
!             2   alpha
!             3   cp
!             4   bulkmodulus Ks
!             5   conductivity
!             6   adiabat
!   For eos_type=1 with ichoice=5 we first store rho, which
!   is converted to conductivity in init_eos03 outside of this routine
subroutine init_eos01(ichoice,coor,npoint,ndim,vec)
use geometry
use convparam
use coeff
implicit none
integer,intent(in) :: npoint,ichoice,ndim
real(kind=8),intent(in) :: coor(ndim,npoint)
real(kind=8),intent(inout) :: vec(:)
integer :: i,find_eos_iz,iz,izmin,izmax
real(kind=8) :: x,y,r,z_d,dz,rl,funccf,adiabat,y0,z,dz_d

izmin=eos_np
izmax=0

do i=1,npoint
   x = coor(1,i)
   y = coor(2,i)
   if (cyl) then
      r = sqrt(x*x+y*y)
!     dimensional depth in km; 
      z_d = (r2-r)*height_dim/1e3
      z   = (r2-r)
   else
      z_d = (1-y)*height_dim/1e3
      z   = (1-y)
   endif
!  if (i==npoint) write(6,*) 'eos_z_d: ',eos_z_d(1),eos_z_d(eos_np),eos_np
   z_d=min(eos_z_d(eos_np),z_d)
   z=min(eos_z(eos_np),z)
!  iz indicates the bottom of the interval in the eos_* lookup tables
!  write(irefwr,*) 'find_eos_iz'
   iz = find_eos_iz(z_d,eos_z_d,eos_np)
   iz = find_eos_iz(z,eos_z,eos_np)
   dz_d = eos_z_d(iz+1)-eos_z_d(iz)
   rl = (z_d-eos_z_d(iz))/dz_d
   dz = eos_z(iz+1)-eos_z_d(iz)
   rl = (z-eos_z(iz))/dz
   if (ichoice == 1) then
!     Compute rho. For ichoice == 5 we'll compute the diffusivity later
      vec(i) = (eos_rho(iz)*(1d0-rl)+eos_rho(iz+1)*rl)
!     if (x<1e-4) write(6,'(''rho: '',i5,4e15.7)') iz,z_d,rl,eos_rho(iz),vec(i)
   else if (ichoice == 2) then
!     alpha
      if (ialphatype > 0) then
         vec(i) = (eos_alpha(iz)*(1d0-rl)+eos_alpha(iz+1)*rl)
      else
         vec(i) = 1d0
      endif
   else if (ichoice == 3) then
!     cp
      vec(i) = (eos_cp(iz)*(1d0-rl)+eos_cp(iz+1)*rl)
   else if (ichoice == 4) then
!     bulk modulus
      vec(i) = (eos_K(iz)*(1d0-rl)+eos_K(iz+1)*rl)
   else if (ichoice == 5) then
!     conductivity
      if (icondtype > 0) then
         vec(i) = (eos_rho(iz)*(1d0-rl)+eos_rho(iz+1)*rl)
      else
         vec(i) = 1d0
      endif
   else if (ichoice == 6) then
!     adiabatic temperature
      vec(i) = (eos_T(iz)*(1d0-rl)+eos_T(iz+1)*rl)
   endif
enddo

end subroutine init_eos01

!   INIT_EOS02
!   Divide components of vec by their average value (vecav)
subroutine init_eos02(coor,npoint,vec,vecav)
use sepmodulecomio
use control
use coeff
implicit none
integer :: npoint
real(kind=8) :: coor(2,npoint),vec(npoint)
integer :: i,find_eos_iz,iz
real(kind=8) :: x,y,r,z,dz,rl,vecav

if (vecav <= 0d0) then
   if (print_node) write(irefwr,*) 'PERROR(init_eos02): vecav <= 0'
   call instop 
endif
 
vec = vec / vecav

end subroutine init_eos02

!   FIND_EOS_IZ
!   Find depth interval in EoS tables for depth z (in km).
integer function find_eos_iz(z,eos_z,eos_np)
use sepmodulecomio
use control
implicit none
integer :: eos_np
real(kind=8) :: z,eos_z(eos_np)
integer :: i
   
if (z <= eos_z(1)) then
   find_eos_iz=1
   return
else if (z >= eos_z(eos_np)) then
   find_eos_iz=eos_np
endif
i=1
do
   i=i+1
   if (z >= eos_z(i-1) .and. z <= eos_z(i) ) then
      find_eos_iz = i-1
      exit
   else if (i >= eos_np) then
      if (print_node) then
         write(irefwr,*) 'PERROR(find_eos_iz): out of range'
         write(irefwr,*) 'i, eos_np = ',i,eos_np
         write(irefwr,*) 'eos_z(1), eos_z(eos_np), z: ',eos_z(1), eos_z(eos_np), z
      endif
      call instop
   endif
enddo
 
end function find_eos_iz

