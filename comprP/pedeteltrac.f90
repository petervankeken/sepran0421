!   PEDETELTRAC
!   Find element information for velocity interpolation
!     - find element in which this tracer sits
!       first by using the ielfromgrid lookup mechanism
!       if this fails, check neighbors of the indicated elements
!       if that fails too, check all elements (which shouldn't be 
!             necessary; a warning will be printed if too many
!             tracers are not found by the first two look up 
!             mechanisms)
!     - get coordinates
!     - check for and fix curvature if necessary
!     - if necessary, interpolate the solution in (xm,ym)
!
!   The mesh is defined using quadratic triangles that may
!   have curved edges. In order to compute the correct shape functions
!   we straighten the curved elements by recomputing the coordinates
!   of the midpoints and reinterpolate the velocity field using quadratic
!   interpolation. The test whether a tracer is within an element
!   should therefore be for straight elements.
!
!   ic        i  choice parameter to determine 
!                   0 = determine element # and coordinates
!                   1 = correct velocity components, element # and coords
!                   2 = correct first velocity component, element # and coords
!                 >10 = return with ielh=-1 if element couldn't be located
!                       rather than stop
!   ioffset   i  extra offset in user_here
!   xm,ym     i     coordinates of the tracer
!   nodno     o     nodal point numbers of this element
!   nodlin    o     nodal point numbers of linear subelement
!   rl        o     barycentric coordinates
!   xn,yn     o     coordinates of nodal point numbers
!   un,vn     o     velocity components in nodal points (filled if ic>0)
!   ielh      o     element number
!   imissed   o     indicator of number of tracers for which the
!                   tracer lookup mechanism fails
! 
!   PvK 220104
!
!   PvK February 2020
!   Updated for version 1219
subroutine pedeteltrac(ic,ioffset,xm,ym,nodno,nodlin,rl,xn,yn,un,vn,ielh,imissed, & 
      & ifound,phiq,xi,eta,guess_first,calling_routine,itrac)
use sepmodulekmesh
use msper01
use mtime
use geometry
use control
use sepran_interface
use sepran_arrays
implicit none
integer :: ic,ioffset
real(kind=8) :: xm,ym,rl(*),xn(*),yn(*),un(*),vn(*),phiq(*)
real(kind=8) :: xi,eta
logical :: guess_first
integer :: nodno(*),nodlin(*),ielh,isub,imissed
integer :: ifound,ichoice,itrac
character(len=*) calling_routine
integer :: i,iel_now(4)=0,j
logical :: midflag,dotflag,elementcurved,first_go_through
logical :: correct,fail
real(kind=8) :: u,func,xi1,eta1,phiq1(7),r
integer,parameter :: NEIGHBORSMAX=1000
integer,dimension(NEIGHBORSMAX) :: current_neighbors,current_neighbors_neighbors
integer :: iel_first,istep=0,ncurrent,ncurrent2
save iel_first,first_go_through,istep

if (kelmo==0) then
   if (print_node) then
      write(irefwr,*) 'PERROR(pedeteltrac): kelmo=',kelmo
      write(irefwr,*) 'make sure to initialize kmesh part o'
   endif
   call instop
endif
   
call checkbounds(1,xm,ym)
iel_first=ielh
first_go_through=.true.
!if (itrac==2645188 .and. time_now>1.7e-5.and.print_node) then
!    verbose=.true.
!else
!    verbose=.false.
!endif
1000  continue

!if (print_node) verbose=.true.

if (verbose.and.print_node) then
   write(irefwr,'('' guess_first, ielh: '',L5,i10)') guess_first,ielh
   open(109,file='coor.dat')
   write(109,*) xm,ym
   close(109)
   open(110,file='elemfirst.dat')
   open(109,file='elem.dat')
endif
if (guess_first.and.ielh.gt.0) then
!    use the guess provided in ielh
     if (ielh>nelem) then
        if (print_node) write(irefwr,*) 'PERROR(pedeteltrac): ielh>nelem: ',ielh,nelem
        call instop
     endif
     call sper01(kmeshc,coor,nodno,xn,yn,ielh)
     if (curved_elem(ielh)) then
        correct=checkinelem(3,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,ielh)
     else 
        correct=checkinelem(1,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,ielh)
     endif
     if (verbose.and.print_node) then
        write(irefwr,*) 'curved_elem = ',curved_elem(ielh),rl(1:3)
        write(irefwr,*) 'correct: ',correct,xi,eta
        write(irefwr,*) 'xm,ym: ',xm,ym
        write(irefwr,*) 'curved_elem: ',curved_elem(ielh)
        write(irefwr,*) 'correct=',correct
     endif
     if (correct) ifound = ifound+1
endif

if (.not.guess_first.or.ielh.le.0.or..not.correct) then
!  Guess is not correct
!  Look up element number from table prepared by ieltogrid()
   ichoice=1
   if (verbose.and.print_node) write(irefwr,'(''call ielfromgrid: '',2f12.3,4i10)') xm,ym,iel_now(1:4)
   call ielfromgrid(ichoice,xm,ym,kmeshc,coor,nodno,xn,yn,iel_now,ielh,phiq,xi,eta)
   if (verbose.and.print_node) write(irefwr,*) 'after ielfromgrid: ielh = ',ielh
   if (verbose.and.print_node) then
     write(109,'(''>'')') 
     do i=1,6
        write(109,*) xn(i),yn(i)
     enddo
     write(109,*) xn(1),yn(1)
     do i=1,6
        write(110,*) xn(i),yn(i)
     enddo
     write(110,*) xn(1),yn(1)
     close(110)
   endif
   if (ielh.gt.0) ifound = ifound + 1
endif

if (ielh.le.0) then
!  if this fails, check neighboring elements
!  Use kmesh part o to check the neighboring elements
   if (verbose.and.print_node) write(irefwr,*) 'check neighbors'
   r=sqrt(xm*xm+ym*ym)
   if (verbose.and.print_node) write(irefwr,'(''call ielfromgrid: '',5e15.7,2L5)') &
      & xm,ym,r,abs(rtop-r),abs(rbot-r),r.ge.rtop,r.le.rbot
   call check_neighbors(xm,ym,kmeshc,coor,kmesho,kmsho2,nelem, &
     &            nodno,nodlin,xn,yn,iel_now,4,ielh,current_neighbors,ncurrent,rl,phiq,xi,eta,verbose)
   if (verbose.and.print_node) then
     write(irefwr,*) 'current neighbors: ',(current_neighbors(i),i=1,ncurrent)
     write(irefwr,*) 'done with neighbors'
   endif
   needed_neighbors=needed_neighbors+1
endif

if (ielh.le.0) then
   ! in pathetic cases the tracer may be just outside the neighbors of the
   ! element picked by ielfromgrid. So check the neighbor's neighbors
   call check_neighbors(xm,ym,kmeshc,coor,kmsho1,kmsho2,nelem, &
     &            nodno,nodlin,xn,yn,iel_now,ncurrent,ielh,current_neighbors_neighbors,ncurrent2,rl,phiq,xi,eta,verbose)
   if (verbose.and.print_node) then
      write(irefwr,*) 'current neighbors neighbors: ',(current_neighbors_neighbors(i),i=1,ncurrent2)
      write(irefwr,*) 'done with neighbors'
   endif
   needed_neighbors_neighbors=needed_neighbors_neighbors+1

endif

! debugging for when nearest nearest neighbors fails
 if (ielh.le.0) then
!   istep=istep+1
!   if (print_node) verbose=.true.
!   if (istep==1) then
!      goto 1000
!   else
!      write(irefwr,*) 'stopping because check_neighbors failed'
!      return
!   endif


!  Resort to desperate measures: loop over all elements. 
   imissed=imissed+1
   ielh=1
100      continue
     call sper01(kmeshc,coor,nodno,xn,yn,ielh)
     ! in this desperate case assume element is straight because Newton
     ! iteration may fail if tracer is very far outside the test element
     correct=checkinelem(1,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,ielh)
!    if (curved_elem(ielh)) then
!       correct=checkinelem(3,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,ielh)
!    else 
!       correct=checkinelem(1,xn,yn,xm,ym,nodno,nodlin,rl,xi,eta,phiq,isub,ielh)
!    endif
     if (verbose.and.print_node) write(irefwr,'(''ielh, xm, ym = '',i5,2f12.3,2L5)') ielh,xm,ym, curved_elem(ielh),correct
     if (.not.correct) then
        ielh=ielh+1
        if (ielh.le.nelem) then
           goto 100
        else if (first_go_through) then
!          debug: do this all again from the start with more output
           if (print_node) verbose=.true.
           ielh=iel_first
           first_go_through=.false.
           goto 1000
        else
!          bail out
           if (print_node) then
              write(irefwr,*) 'PWARN(pedeteltrac): did not find correct' 
              write(irefwr,*) 'element for tracer : ',ielh,xm,ym,nelem
              write(irefwr,*) 'called from: ',calling_routine
              write(irefwr,*) 'r: ',sqrt(xm*xm+ym*ym)
              open(9,file='failed')
              write(9,*) 'T'
              close(9)
           endif
           call instop()
         endif
     endif
endif

if (ic.eq.1) then
   call sper06(user_here(ioffset+10:ioffset+10+npoint-1),un,nodno)
   call sper06(user_here(ioffset+10+npoint:ioffset+10+2*npoint-1),vn,nodno)
else if (ic.eq.2) then
   call sper06(user_here(ioffset+10:ioffset+10+npoint-1),un,nodno)
   do i=1,6
      vn(i)=0
   enddo
else 
   do i=1,6
      un(i)=0
      vn(i)=0
   enddo
endif

verbose=.false.
   
end subroutine pedeteltrac

