!   VELOCITY_BC
!
!   Enforce essential boundary conditions for velocity
!
!   PvK 072902
subroutine velocity_bc(ichoice)
use sepmoduleoldrouts
use sepmodulecomio
use sepran_arrays
use geometry
use coeff
use brandenburg
use control
implicit none
integer :: ibp,iall,ichoice,ic,itarget
real(kind=8) :: uval=0.0_8
logical :: exists
character(len=80) :: fname


if (.not.gable_plates) then
   ! deal with default cases first
   if (its_CH94 .and. .not.CH94_fs) then
      ! set essential b.c. at curve 3 using funcbc(10)
      call bvalue(10,2,kmesh1,kprob1,isol1,uval,3,3,1,0)
      return
   endif


   inquire(file='top_velocity',exist=exists)
   if (exists) then
      open(9,file='top_velocity')
      read(9,*) ic,uval
   !  set up top velocity b.c.
      call bvalue(0,1,kmesh1,kprob1,isol1,uval,ic,ic,1,0)
      close(9)
      if (print_node) write(irefwr,'('' set velocity at curve '',i5,'' to '',f12.3)') ic,uval
   endif

   if (ichoice.eq.0) then
   ! Set normal velocity to zero in selected curves
     do ibp=1,iclc
        iall=icloc(ibp)
        call bvalue(0,2,kmesh1,kprob1,isol1,0d0,iall,iall,1,0)
     enddo
   else if (ichoice == 10) then
      ! Set normal velocity to zero in selected curves
     do ibp=1,iclc
        iall=icloc(ibp)
        itarget=isolold1(1)
        call bvalue(0,2,kmesh1,kprob1,itarget,0d0,iall,iall,1,0)
     enddo

   endif

   do ibp=1,nbp
      !  Set velocity to zero in selected points
      iall=iboundpoints(ibp)
      call bvalue(0,1,kmesh1,kprob1,isol1,0d0,iall,iall,1,0)
      call bvalue(0,1,kmesh1,kprob1,isol1,0d0,iall,iall,2,0)
   enddo
   
   return

else  ! gable plates
   ! write(fname,'(''funcbc.'',i1)') gable_stokes_choice
   ! open(109,file=fname)
   if (gable_stokes_choice>0) then
      ! set up boundary conditions in the test solutions
      do ibp=1,iclc
        iall=icloc(ibp)
        call bvalue(0,2,kmesh1,kprob1,plate_sol(gable_stokes_choice),0d0,iall,iall,1,0)
     enddo
     ! now set uniform velocity along this plate
     call bvalue(3,2,kmesh1,kprob1,plate_sol(gable_stokes_choice),0d0,itop,itop,2,0)
     ! and point the points
     do ibp=1,nbp
     !  Set velocity to zero in selected points
        iall=iboundpoints(ibp)
        call bvalue(0,1,kmesh1,kprob1,plate_sol(gable_stokes_choice),0d0,iall,iall,1,0)
        call bvalue(0,1,kmesh1,kprob1,plate_sol(gable_stokes_choice),0d0,iall,iall,2,0)
     enddo
   else  if (gable_stokes_choice==0) then
     ! do the same but now for isol
     do ibp=1,iclc
        iall=icloc(ibp)
        call bvalue(0,2,kmesh1,kprob1,isol1,0d0,iall,iall,1,0)
     enddo
     do ibp=1,nbp
        iall=iboundpoints(ibp)
        call bvalue(0,1,kmesh1,kprob1,isol1,0d0,iall,iall,1,0)
        call bvalue(0,1,kmesh1,kprob1,isol1,0d0,iall,iall,2,0)
     enddo
     call bvalue(0,2,kmesh1,kprob1,isol1,0d0,itop,itop,2,0)
   endif

endif ! gable_plates
      

end subroutine velocity_bc

