      subroutine sepel900rhsd
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
!   or be transmitted by any means, electronic, mechanical, photocopy,
!   recording, or otherwise, without written permission from the
!   publisher.
! **********************************************************************
!
!                       DESCRIPTION
!
!     Compute element vector for Navier-Stokes element due to external
!     forces
!     Old subroutine elm900rhsd
! **********************************************************************
!
!                       KEYWORDS
!
!     element_vector
!     navier_stokes_equation
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulecactl
      use sepmoduleelm
      use sepmodulebasefn
      use sepmodulecomio
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer :: i, iseq, kdiag
      double precision :: rhoc(m), work(m)
      double precision, pointer :: vec(:), phi(:,:)

!     i              Counting variable
!     iseq           Sequence number
!     kdiag          Indicates if integration points are equal to nodal
!                    points (1) or not (2)
!     phi            values of shape functions in quadrature points
!     rhoc           w * rho
!     vec            work array to store a vector
!     work           Work array
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     PRINRL         Print 1d real vector
!     PRINRL1        Print 2d real vector
!     SEPELLOAD2     Fill the element vector with a load
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
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen ( 'sepel900rhsd' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepel900rhsd'
         call prinrl ( E%rho, m, 'rho' )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

!     --- Multiply weights by density rho and store in rhoc

      rhoc = B%w * E%rho

      if ( E%methodupw==0 ) then
         phi => B%phi
         kdiag = B%jdiag
      else
         phi => B%psi
         kdiag = 2
      end if  ! ( E%methodupw==0 )

!     --- each component

      do i = 0, E%ndim-1

         iseq = E%indsource+i

         if ( E%ind(iseq)/=0 ) then

!        --- f1 # 0  store contribution temporary in array work

            vec => E%subvec(i+1:E%ndim*n:E%ndim)
            work = rhoc*E%coeffs(:,iseq)
            call sepelload2 ( work, kdiag, phi, vec, iseq)

            if ( debug ) call prinrl ( vec, n, 'work/1' )

         end if  ! ( E%ind(iseq)/=0 )

      end do  ! i = 0, E%ndim-1

1000  if ( debug ) then

!     --- Debug information

         call prinrl1 ( E%subvec, n, E%ndim, 'elemvc/1' )
         write(irefwr,*) 'End sepel900rhsd'

      end if  ! ( debug )
      call erclos ( 'sepel900rhsd' )

      end subroutine sepel900rhsd
