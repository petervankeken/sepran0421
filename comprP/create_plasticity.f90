subroutine create_plasticity
use sepmodulecomio
use sepmodulevecs
use sepmodulekmesh
use sepran_arrays
use coeff
use control
implicit none
interface
   subroutine create_plasticity2(ndim,npoint,coor,eta,etalin,etaplast,plasticity)
   integer,intent(in) :: ndim,npoint
   real(kind=8),intent(in) :: coor(ndim,npoint),eta(npoint),etalin(npoint),etaplast(npoint)
   real(kind=8),intent(out) :: plasticity(npoint) 
   end subroutine create_plasticity2
end interface

if (iplasticity==0 .or. iviscplast==0 .or. iviscplast==0) then
   if (print_node) then
      write(irefwr,*) 'PERROR(create_plasticity): one of three essential vectors is uncreated'
      write(irefwr,*) 'iplasticity = ',iplasticity
      write(irefwr,*) 'iviscplast  = ',iviscplast
      write(irefwr,*) 'ivisclin    = ',ivisclin
      write(irefwr,*) 'These should be non-zero'
   endif
   call instop
endif

call sepactsolbf1(isol2)
call sepgetmeshinfo(ndim,npoint)
pedebug=.true.
if (pedebug) write(irefwr,*) 'isol2: ',ivisclin,iviscplast,iplasticity
call create_plasticity2(ndim,npoint,coor,ks(ivisc)%sol,ks(ivisclin)%sol, &
     &   ks(iviscplast)%sol,ks(iplasticity)%sol)
pedebug=.false.


end subroutine create_plasticity

subroutine create_plasticity2(ndim,npoint,coor,eta,etalin,etaplast,plasticity)
use control
use sepmodulecomio
implicit none
integer, intent(in) :: ndim,npoint
real(kind=8),intent(in) :: coor(ndim,npoint),eta(npoint),etalin(npoint),etaplast(npoint)
real(kind=8),intent(out) :: plasticity(npoint)

plasticity = 1.0_8/(1.0_8/etalin + 1.0_8/etaplast)
plasticity = plasticity/etaplast
if (print_node) then
   write(irefwr,*) 'eta_lin    min/max: ',minval(etalin),maxval(etalin)
   write(irefwr,*) 'eta_plast  min/max: ',minval(etaplast),maxval(etaplast)
   write(irefwr,*) 'plasticity min/max: ',minval(plasticity),maxval(plasticity)
endif

end subroutine create_plasticity2

