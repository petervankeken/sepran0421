      subroutine sepelload2 ( source, kdiag, phi, vec , iseq)
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 13-09-2017
!
!   copyright (c) 2017-2017  "Ingenieursbureau SEPRA"
!   permission to copy or distribute this software or documentation
!   in hard copy or soft copy granted only by written license
!   obtained from "Ingenieursbureau SEPRA".
!   all rights reserved. no part of this publication may be reproduced,
!   stored in a retrieval system ( e.g., in memory, disk, or core)
!   or be transmitted by any means, electronic, mechanical, photocopy
!   recording, or otherwise, without written permission from the
!   publisher.
! **********************************************************************
!
!                       DESCRIPTION
!
!     Fill the element vector with a load
!     Old subroutine el3003
! **********************************************************************
!
!                       KEYWORDS
!
!     element_vector
!     elasticity_equation
! **********************************************************************
!
!                       MODULES USED
!
      use mpetrac
      use tracers
      use geometry
      use coeff
      use sepmodulecactl
      use sepmodulecomio
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer, intent(in) :: kdiag,iseq
      double precision, intent(in) :: source(m), phi(n,m)
      double precision, intent(out) :: vec(n)

!     kdiag          i    if 1 the basis functions are 0 or 1 in the quadrature
!                         points
!                         if  0 they are general
!     phi            i    Contains the values of the basis functions in the
!                         quadrature points
!     source         i    right-hand side
!     vec            o    vector to be filled
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer :: i, k
      double precision :: h

!     h              Help variable to store a constant temporarily
!     i              General loop variable
!     k              Counting variable
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     PRINRL         Print 1d real vector
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
! **********************************************************************
!
!                       PSEUDO CODE
!
!     The part of the element matrix is given by
!
!     / source phi_i  d Omega
!
!     The numerical values are given by
!
!     sum_k w(k) {source phi }(x_k)
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen ( 'sepelload2' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepelload2'
         write(irefwr,1) 'kdiag', kdiag
         call prinrl ( source, m, 'source' )
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      n = n
      if ( debug ) write(irefwr,1) 'n, kstep', n, kstep

      if ( kdiag==1 ) then

!     --- Special case: phi_i(x_j) = delta_ij

         if ( kstep==1 ) then
            vec = vec+source
         else
            do k = kstep, n, kstep
               vec(k) = vec(k)+source(k)
            end do  ! k = kstep, n, kstep
         end if  ! ( kstep==1 )

      else

!     --- General case
 
         !write(6,*) 'going to find rhsd contribution in sepelload2'
         !call instop
         ! PvK: source is rhsd multiplied by density. Make sure to prepare f1/2_tracrhsd with background density too.
         do i = 1, n
            h = sum(source(:)*phi(i,:))
            vec(i) = vec(i)+h
         end do  ! i = 1, n
         if (cyl.and.iseq==4) then
            if ((itracoption==2.and.fieldC) .or. itracoption==2) then
                vec(1:n)=vec(1:n)+f1_tracrhsd(1:n,ielem,1)
            endif
         endif   
         if (iseq==5) then
            if ((itracoption==1.and.tracerC).or.itracoption==2) then
               vec(1:n)=vec(1:n)+f2_tracrhsd(1:n,ielem,1)
            endif
         endif

      end if  ! ( kdiag==1 )

1000  if ( debug ) then

!     --- Debug information

         call prinrl ( vec, n, 'element vector' )
         write(irefwr,*) 'End sepelload2'

      end if  ! ( debug )
      call erclos ( 'sepelload2' )

      end subroutine sepelload2
