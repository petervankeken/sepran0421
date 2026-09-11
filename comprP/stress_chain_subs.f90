! strchain: A collection of subroutines to compute stress along  
! a marker chain at the base of the lithosphere.       
! J.P. Brandenburg, UofM, 2005/2006                   
!  Initialize the coordinates, tensor, and stress for the stress chain
subroutine strchain_start
use sepmodulecomio
use sepran_arrays
use control
use brandenburg
use geometry
use coeff
implicit none

integer, parameter :: NPMAX=10000
real(kind=8) :: x,y,dx,rdum,ctest
real(kind=8) :: thtest1,thtest2,dth,tangle,test
integer :: i,j,k,idum,icurvs(2),idx
save icurvs

! The methods for setting up the stress chain is somewhat different for
! cartesian and cylindrical coordinates:

if (cyl) then
   ! Setup stress-chain in cylindrical coordinates
   if (print_node) write(irefwr,*)'(in cylindrical) strcurve = ',istress_curve
   ! Use compcr to get nodal points along the stress curve ******
   icurvs(1) = 0
   icurvs(2) = istress_curve
   call compcr(-1,kmesh1,kprob1,isol1,1,icurvs,funcx,funcy)

   ! Number of nodes:
   nsmark = int(funcx(5)/2.0)
   if (print_node) write(irefwr,*) 'nsmark = ',nsmark
   if (2*nsmark > NCHAINMAX) then
      write(irefwr,*) 'PERROR(strchain_start): nsmark too large: ',nsmark
      write(irefwr,*) 'adjust NCHAINMAX = ',NCHAINMAX
      call instop 
   endif

   ! Move the nodal coordinates to the stress chain array:
   idx = 6
   do i=1,2*nsmark
      xstrchain(i) = funcx(idx)
      ystrchain(i) = funcx(idx+1)
      idx = idx + 2
   enddo

   strc_r = DSQRT(xstrchain(1)*xstrchain(1) + ystrchain(1)*ystrchain(1))

   if (print_node) write(irefwr,*)'strc_r : ',strc_r

   ! For now, link the boundaries between plates to the nearest node.
   ! This changes the input plate boundary angles slightly, but 
   ! circumvents problems associated with using intcoor in cylindrical coordinates.
   do i=1,nplates+1
      xbpoint(i) = strc_r*dsin(plate_boundaries(i))
      ybpoint(i) = strc_r*dcos(plate_boundaries(i))
      if (print_node) write(irefwr,*) 'xbpoint : ',xbpoint(i),' ybpoint : ',ybpoint(i)
   enddo

   ! Double-check to make sure this is correct.  Due to variations
   ! in mesh building methods, the curve may be filled right to 
   ! left, or left to right.
   thtest1 = atan(ystrchain(2)/xstrchain(2))
   thtest2 = atan(ystrchain(3)/xstrchain(3))
   if (print_node) write(irefwr,*) 'thtest1 = ',thtest1
   if (print_node) write(irefwr,*) 'thtest2 = ',thtest2

   if (thtest2.ge.thtest1) then
      if (print_node) write(irefwr,*)'chain in wrong order: change splate_sol of strcurve'
      call instop
   endif

   dth = thtest1 - thtest2
   if (print_node) write(irefwr,*)'dtheta stress chain = ',dth

   ! Work out the polar angle (clockwise from 'North') for each node.
   ! Note that the final marker required some special treatment...
   ! Without this, in a full cylinder the algorithm would assplate_sol
   ! an angle of 0 rather than 2pi to the final marker.
   do i=1,nsmark
      tangle = acos(ystrchain(i)/strc_r)
      if (xstrchain(i).lt.0.0_8) then
         tangle = 2.0_8*pi - tangle
      endif
      thetachain(i) = tangle 
   enddo
   if (thetachain(nsmark-1).ge.4.7_8) then
      thetachain(nsmark) = 2.0_8*pi
   endif

   if (print_node) write(irefwr,*)'stress chain starts at theta = ',thetachain(1)
   if (print_node) write(irefwr,*)'stress chain ends at theta = ',thetachain(nsmark)
 
        
   !call strchain_out
else
    ! Setup stress-chain in cartesian coordinates 

    if (print_node) write(irefwr,*)'strcurve = ',istress_curve
    ! Use compcr to get nodal points along the stress curve ******
    icurvs(1) = 0
    icurvs(2) = istress_curve
    call compcr(-1,kmesh1,kprob1,isol1,1,icurvs,funcx,funcy)

    nsmark = int(funcx(5)/2.0)
    if (print_node) write(irefwr,*)' nsmark = ',nsmark
    if (2*nsmark > NCHAINMAX) then
      if (print_node) then
        write(irefwr,*) 'PERROR(strchain_start): nsmark too large: ',nsmark
        write(irefwr,*) 'adjust NCHAINMAX = ',NCHAINMAX
      endif
      call instop 
    endif

    ! Depending on the order that the mesh was constructed, compcr might
    ! go through the points in reverse order.  Double check this --->
    ctest = funcx(8) - funcx(6)
    if (ctest.le.0.0_8) then
       if (print_node) write(irefwr,*)'stress curve in reverse order - change plate_sol of istress_curve'
       call instop
    endif


    idx = 6
    do i=1,2*nsmark
       xstrchain(i) = funcx(idx)
       ystrchain(i) = funcx(idx+1)
       idx = idx + 2
    enddo

    ! call strchain_out
endif

end subroutine strchain_start

subroutine strchain_out
use control
use geometry
use brandenburg
use coeff
implicit none

real(kind=8) :: x,y,alph,dalph,r
integer :: i,j,k

! include 'dpi.inc'
! include 'strchain.inc'
open(50,file='schain.xy')
do i=1,nsmark
   write(50,*) xstrchain(i),'  ',ystrchain(i)
enddo
close(50)

end subroutine strchain_out

! Link the plate boundaries to specific nodes - 
subroutine strchain_analyze()
use sepmodulecomio
use geometry
use control
use brandenburg
use coeff
implicit none

integer :: i,j,k
real(kind=8) :: spctest

! include 'strchain.inc'
! include 'ccc.inc'

if (cyl) then
   if (print_node) write(irefwr,'(''plate boundaries: '',8f10.5)') plate_boundaries(1:nplates)
   ! Cylindrical Linking 
   do i=1,nplates
      ! 1) Loop through tracers and find lower bound
      do j=1,nsmark
         if (thetachain(j).gt.plate_boundaries(i)) then
            intlink(1,i) = j
            exit
         endif
      enddo
      ! 2) Loop through tracers and find upper bound
      do j=1,nsmark
         if (thetachain(j).ge.plate_boundaries(i+1)) then
            intlink(2,i) = j - 1
            exit 
         endif
      enddo
   enddo

   intlink(1,1) = 1
   intlink(2,nplates) = nsmark
   ! write(irefwr,*) 'strchain_analyze: '
   ! do i=1,nplates
        if (print_node) write(irefwr,'(2i5)') intlink(1,i),intlink(2,i)
   ! enddo


   ! Loop through again, and link down to nearest node.  Adjust the plate configuration to suite.
   do i=2,nplates
      intlink(1,i) = intlink(1,i)-1
      plate_boundaries(i) = thetachain(intlink(1,i))
   enddo
   ! do i=1,nplates
   !       write(irefwr,'(2i5,f15.5)') intlink(1,i),intlink(2,i),plate_boundaries(i)
   ! enddo


else 
   ! Cartesian Linking 

   do i=1,nplates
      ! 1 Loop through tracers and find lower bound
      do j=1,nsmark
         if (xstrchain(j).gt.plate_boundaries(i)) then
            intlink(1,i) = j 
            exit
         endif
      enddo
      ! 2) Loop through tracers and find upper bound
      do j=1,nsmark
         if (xstrchain(j).ge.plate_boundaries(i+1)) then
            intlink(2,i) = j - 1
            exit
         endif
      enddo
    enddo

    intlink(1,1) = 1
    intlink(2,nplates) = nsmark 

    ! Just for CH1994 - the slightly offset point between the plates *****
    xbpoint(1) = plate_boundaries(2)
    ybpoint(1) = ystrchain(1)
    if (print_node) write(irefwr,*)'Boundary point at : ',xbpoint(1),'  ',ybpoint(1)
endif


do i=1,nplates
   if (print_node) write(irefwr,*) 'Plate Boundaries ',i,' ~ Markers # ',intlink(1,i),'  ',intlink(2,i)
enddo

end subroutine strchain_analyze


subroutine strchain_compute(user,iuser,iplate,iba)
use sepmodulecomio
use sepran_arrays
use geometry
use brandenburg
use coeff
use control
implicit none


real(kind=8) :: user(*)
real(kind=8) :: xm,ym,xn(6),yn(6),un(6),vn(6),done
real(kind=8) :: tau,theta,btten(3)
real(kind=8) :: tn(6),temp,visc,pefvis,rdum
real(kind=8) :: viscmax,viscmin
integer :: idudy(5),idvdx(5),idudx(5),idvdy(5)
integer :: ihelp,iuser(*),iinmap(1),map(5),iba

! For getting the stress at the boundary point
real(kind=8) :: coorchain(2,10),derchain(1,10)
real(kind=8) :: bdudx,bdudy,bdvdx,bdvdy

integer :: icurvs(2),npts,idx
integer,parameter :: NPMAX=10000
real(kind=8) :: funcy1(NPMAX),funcy2(NPMAX)
real(kind=8) :: fundudx(NPMAX),fundudy(NPMAX)
real(kind=8) :: fundvdx(NPMAX),fundvdy(NPMAX)
real(kind=8) :: dudx(NPMAX),dudy(NPMAX),dvdx(NPMAX),dvdy(NPMAX)

integer :: npoint,nelem,nelgrp,ikelmc,ikelmi,ikelmo
integer :: iniget,inidgt,i,j,k,nodno(6),iplate,idum

save funcy1,funcy2,fundudx,fundudy,fundvdx,fundvdy
save idudy,idvdx,idudx,idvdy,dudx,dudy,dvdx,dvdy,coorchain,derchain,iinmap,map
save xn,yn,un,vn,btten

!     include 'strchain.inc'
!     include 'mysepar.inc'
!     include 'cjp.inc'
!     include 'c1visc.inc'
!     include 'ccc.inc'


if (iplate.eq.0) then
   if (iba.eq.0) then
      open(74,file='tau_before.out')
   else
      open(74,file='tau_after.out')
   endif
endif

if (cyl) then
    if (iplate.eq.0) then
      ! Compute du/dx using DERIVA:  deriva( * , * , xi, ui, ***)
      call deriva(2,1,1,1,1,idudx,kmesh1,kprob1,isol1,isol1,iuser,user,ihelp)
      ! Compute du/dy using DERIVA:  deriva( * , * , xi, ui, ***)
      call deriva(2,1,2,1,1,idudy,kmesh1,kprob1,isol1,isol1,iuser,user,ihelp)
      ! Compute dv/dx using DERIVA
      call deriva(2,1,1,2,1,idvdx,kmesh1,kprob1,isol1,isol1,iuser,user,ihelp)
      ! Compute dv/dx using DERIVA
      call deriva(2,1,2,2,1,idvdy,kmesh1,kprob1,isol1,isol1,iuser,user,ihelp)
    else
       ! For force Balance Test Solutions:
       ! Compute du/dx using DERIVA:  deriva( * , * , xi, ui, ***)
       call deriva(2,1,1,1,1,idudx,kmesh1,kprob1,plate_sol(iplate),plate_sol(iplate),iuser,user,ihelp)
       ! Compute du/dy using DERIVA:  deriva( * , * , xi, ui, ***)
       call deriva(2,1,2,1,1,idudy,kmesh1,kprob1,plate_sol(iplate),plate_sol(iplate),iuser,user,ihelp)
       ! Compute dv/dx using DERIVA
       call deriva(2,1,1,2,1,idvdx,kmesh1,kprob1,plate_sol(iplate),plate_sol(iplate),iuser,user,ihelp)
       ! Compute dv/dx using DERIVA
       call deriva(2,1,2,2,1,idvdy,kmesh1,kprob1,plate_sol(iplate),plate_sol(iplate),iuser,user,ihelp)
    endif

    ! Use compcr to find the value along istress_curve ****************************
    icurvs(1) = 0
    icurvs(2) = istress_curve
    fundudx(1) = NPMAX
    fundudy(1) = NPMAX
    fundvdx(1) = NPMAX
    fundvdy(1) = NPMAX

    call compcr(0,kmesh1,kprob1,idudx,1,icurvs,funcx,fundudx)
    call compcr(0,kmesh1,kprob1,idudy,1,icurvs,funcx,fundudy)
    call compcr(0,kmesh1,kprob1,idvdx,1,icurvs,funcx,fundvdx)
    call compcr(0,kmesh1,kprob1,idvdy,1,icurvs,funcx,fundvdy)

    do j=1,nsmark
       dudx(j) = fundudx(j+5)
       dudy(j) = fundudy(j+5)
       dvdx(j) = fundvdx(j+5)
       dvdy(j) = fundvdy(j+5)
    enddo

    do j=1,nsmark
       !  construct the strain rate tensor ~
       stten(1,j) = 2.0d0*dudx(j)
       stten(2,j) = dudy(j) + dvdx(j)
       stten(3,j) = 2.0d0*dvdy(j)
       tau = 0.5d0*(stten(3,j)-stten(1,j))*sin(2.0d0*thetachain(j))-stten(2,j)*cos(2.0d0*thetachain(j))
       sttau(j) = tau

!      if (iplate.eq.0) then
!          write(74,*)thetachain(j),'  ',sttau(j)
!      endif
    enddo

    ! Get the stress in the boundary points: ********************
!      do j=1,nplates+1
!         iinmap(1) = 0
!         coorchain(1,1) = xbpoint(j)
!         coorchain(2,1) = ybpoint(j)
!
!         call intcoor(kmesh1,kprob1,idudx,derchain,coorchain,1,1,2,
!     v             iinmap,map)
!         bdudx = derchain(1,1)
!
!         call intcoor(kmesh1,kprob1,idudy,derchain,coorchain,1,1,2, iinmap,map)
!         bdudy = derchain(1,1)
!         call intcoor(kmesh1,kprob1,idvdx,derchain,coorchain,1,1,2, iinmap,map)
!         bdvdx = derchain(1,1)
!
!         call intcoor(kmesh1,kprob1,idvdy,derchain,coorchain,1,1,2,iinmap,map)
!         bdvdy = derchain(1,1)
!         btten(1) = 2.0d0*bdudx
!         btten(2) = bdudy + bdvdx
!         btten(3) = 2.0d0*bdvdy
!
!         bttau(j) = 0.5d0*(btten(3)-btten(1)*sin(2.0d0*plate_boundaries(j)) - btten(2)*cos(2.0d0*plate_boundaries(j)) )
!      enddo

!      if (gable_stokes_choice.eq.0) then
!        do j=1,nplates+1
!           write(71,*)plate_boundaries(j),'  ',bttau(j)
!        enddo
!      endif

else 
   ! Cartesian Formulation ************************************************************
   if (iplate.eq.0) then
      ! Compute du/dy using DERIVA:  deriva( * , * , xi, ui, ***)
      call deriva(2,1,2,1,1,idudy,kmesh1,kprob1,isol1,isol1,iuser,user,ihelp)
      ! Compute dv/dx using DERIVA
      call deriva(2,1,1,2,1,idvdx,kmesh1,kprob1,isol1,isol1,iuser,user,ihelp)
   else
      ! Compute du/dy using DERIVA:  deriva( * , * , xi, ui, ***)
      call deriva(2,1,2,1,1,idudy,kmesh1,kprob1,plate_sol(iplate),plate_sol(iplate),iuser,user,ihelp)
      ! Compute dv/dx using DERIVA
      call deriva(2,1,1,2,1,idvdx,kmesh1,kprob1,plate_sol(iplate),plate_sol(iplate),iuser,user,ihelp)
   endif
   ! Use compcr to find the value along curve 13 ****************************
   icurvs(1) = 0
   icurvs(2) = istress_curve
   funcy1(1) = NPMAX
   funcy2(1) = NPMAX

   call compcr(0,kmesh1,kprob1,idudy,1,icurvs,funcx,funcy1)
   call compcr(0,kmesh1,kprob1,idvdx,1,icurvs,funcx,funcy2)
   npts = nint(funcy1(5))
   if (npts.ne.nsmark) then
      if (print_node) write(irefwr,*)'something wrong in strchain'
      call instop
   endif
   do j=1,nsmark
      dudy(j) = funcy1(j+5)
      dvdx(j) = funcy2(j+5)
      ! for Cartesian:
      dudx(j) = 0.0d0
      dvdy(j) = 0.0d0
   enddo

   ! Get the stress in the boundary point: ********************
   iinmap(1) = 0
   coorchain(1,1) = xbpoint(1)
   coorchain(2,1) = ybpoint(1)
   call intcoor(kmesh1,kprob1,idudy,derchain,coorchain,1,1,2,iinmap,map)
   bdudy = derchain(1,1)
   call intcoor(kmesh1,kprob1,idvdx,derchain,coorchain,1,1,2,iinmap,map)
   bdvdx = derchain(1,1)
   sttau(NPMAX/5) = bdudy + bdvdx
   do j=1,nsmark
      ! construct the strain rate tensor ~
      stten(1,j) = 2.0d0*dudx(j)
      stten(2,j) = dudy(j) + dvdx(j)
      stten(3,j) = 2.0d0*dvdy(j)
      ! Now, calculate shear stress in the markers ~
      ! Since this is set up for cartesian in the simplest case,
      ! set theta = 0 for the time being:
       theta = 0.0d0
       tau = 0.5d0*(stten(3,j)-stten(1,j))*sin(2.0d0*theta)+ stten(2,j)*cos(2.0d0*theta)
       sttau(j) = tau
   enddo
   ! Calculate true shear stress instead of strainrate
   if (FBvisc) then
      if (print_node) write(irefwr,*)'FBvisc - not ready now'
      call instop
   endif
endif

close(74)

end subroutine strchain_compute

function strchain_integrate(iplate)
use geometry
use brandenburg
use coeff
implicit none

integer,parameter :: NPMAX=10000
real(kind=8) :: taui,h,strchain_integrate,rampint,r
real(kind=8) :: hc1,hc2,tcor1,tcor2
integer :: iplate,istart,iend,i

taui = 0.0d0

if (cyl) then
   ! Integration in cylindrical coordinates
   r = strc_r

   if (iplate.eq.0) then
      ! Integrate over the whole chain:
      istart = 1
      iend =(nsmark - 1)
      !write(irefwr,*) 'integrate from/to ',istart,iend
      do i=istart,iend
         h = thetachain(i+1) - thetachain(i)
         h = h*r
         taui = taui + h*(0.5d0*sttau(i) + 0.5d0*sttau(i+1))
      enddo
   else
      istart = intlink(1,iplate)
      iend   = intlink(2,iplate)-1
      !write(irefwr,*) 'integrate from/to ',istart,iend
      do i=istart,iend
         h = thetachain(i+1) - thetachain(i)
         h = h*r
         taui = taui + h*(0.5d0*sttau(i) + 0.5d0*sttau(i+1))
      enddo

      ! Apply corrector for when plate boundary is not exactly on a node:
      ! hc1 = thetachain(istart) - plate_boundaries(iplate)
      ! hc2 = plate_boundaries(iplate+1) - thetachain(iend+1)
      ! hc1 = r*hc1
!           hc2 = r*hc2
!           if (iplate.eq.1) then 
!              hc1 = 0.0d0
!           endif
!
!           if (iplate.eq.nplates) then
!              hc2 = 0.0d0
!           endif
 
!           write(irefwr,*)'for plate # ;',iplate
!           write(irefwr,*)'thetachain(istart) : ',thetachain(istart)
!           write(irefwr,*)'thetachain(iend+1)   : ',thetachain(iend+1)
!           write(irefwr,*)'bounds: ',plate_boundaries(iplate),'  ',plate_boundaries(iplate+1)
!           write(irefwr,*)'hc    : ',hc1,'   ',hc2          


! 'lower' end corrector:
!           tcor1 = hc1*(0.5d0*sttau(istart) + 0.5d0*bttau(iplate))
!           tcor1 = hc1*(sttau(istart))
! 'upper' end corrector:
!           tcor2 = hc2*(bttau(iend+1))
 
!           taui = taui + tcor1 + tcor2
!         
     endif

     strchain_integrate = taui

else 

     ! Cartesian 
     if (iplate.eq.0) then
        !  Integrate over the whole chain:
        istart = 1
        iend =(nsmark - 1)
        do i=istart,iend
           h = xstrchain(i+1) - xstrchain(i)
           taui = taui + h*(0.5d0*sttau(i) + 0.5d0*sttau(i+1))
        enddo
     else 

        ! For the moment, the plate boundaries are only precise for two plates.
        ! Future adaptations will have to add a plate boundary correction for 
        ! more than two plates.

         istart = intlink(1,iplate)
         iend   = intlink(2,iplate)-1
         do i=istart,iend
            h = xstrchain(i+1) - xstrchain(i)
            taui = taui + h*(0.5d0*sttau(i) + 0.5d0*sttau(i+1))
         enddo

         ! Apply correction ~ For when the plate boundary is not right on a node
         if (iplate.eq.1) then
            h = xbpoint(1) - xstrchain(iend+1)
            taui = taui + h*(0.5d0*sttau(iend+1) + 0.5d0*sttau(NPMAX/2))
         else if (iplate.eq.2) then
            h = xstrchain(istart) - xbpoint(1)
            taui = taui + h*(0.5d0*sttau(istart) + 0.5d0*sttau(NPMAX/2))
         endif
     endif

     strchain_integrate = taui

endif 

end function strchain_integrate

