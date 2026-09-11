! FIND_MELT_ZONES
! Figure out which of the nplates zones are melt zones
! by testing whether velocity of the right plate is faster
! than the velocity of the left plate
! PvK Aug 2 2019
subroutine find_melt_zones
use sepmodulecomio
use brandenburg
use control
implicit none
integer :: iplate,ip,ip1,ip2,izone
real(kind=8),dimension(10) :: divV

! A melt zone is defined when plates strongly diverge : v1*v2 < 0
! or weakly diverge : v1,v2<0 and abs(v1)>abs(v2) 
!                     v1,v2>0 and abs(v2)>abs(v1)
! In all cases v2>v1 with velocity going in the clock-wise direction
! or v2<v1 with velocity in the counter-clock-wise direction 
! (which is how JP seems to have defined it)

diverging=.false.
ip=0
nmeltzone=0
ningaszone=0
ingaszone=0
meltzone=0
! zone 1 is north pole at azimuth 0
do izone=1,nplates
   ip1=izone-1
   ip2=izone
   if (izone==1)  then
      ! polar zone has special coordinates
      ip1=nplates
      ip2=1
   endif
   divV(izone)=plate_vel(ip2)-plate_vel(ip1)
   if (plate_vel(ip2)<plate_vel(ip1)) then
      diverging(izone)=.true.
      nmeltzone=nmeltzone+1
      meltzone(nmeltzone)=izone
   else
      ningaszone=ningaszone+1
      ingaszone(ningaszone)=izone
   endif
enddo
if (gable_output_choice>0.and.print_node) then
   write(irefwr,'(''melting now in zones: '',9i8,:)') (meltzone(ip),ip=1,nmeltzone) 
   write(irefwr,'(''diverging           : '',4x,8L8,:)') diverging(1:nplates)
   write(irefwr,'(''V                   : '',9i8,:)') nint(plate_vel(nplates)*0.1),(nint(plate_vel(ip)*0.1),ip=1,nplates)
   write(irefwr,'(''divV                : '',4x,8i8,:)') (nint(divV(ip)*0.1),ip=1,nplates)
endif

end subroutine find_melt_zones
