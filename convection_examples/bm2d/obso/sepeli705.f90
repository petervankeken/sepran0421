      double precision function sepeli705 ( )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 07-06-2017
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
!     Compute element integral for problems with one degree of freedom
!     per point
! **********************************************************************
!
!                       KEYWORDS
!
!     element_integral
!     integral
!     elliptic_equation
!     parabolic_equation
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulebuild
      use sepmodulekprob
      use sepmoduleelm
      use sepmoduleelem
      use sepmodulebasefn
      use sepmoduleoldvecs
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
!

!     sepeli705      o    Computed element integral
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
      double precision :: sepeli705int

!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     PRINRL         Print 1d real vector
!     SEPELBASEFUNC  Compute information about basis functions
!     SEPELDEALLOC   Deallocate coefficients and basis functions for element
!                    subroutines
!     SEPELFILLVAR   Fill space dependent coefficients
!     SEPELI0705     Perform initialization for subroutine sepeli705
!     SEPELI705INT   Compute element integral for problems with one degree
!                    of freedom per point
!     SEPELOLDVECS   Fill old vectors for element
!     SEPELOLDVECSK  Fill old vectors for element (only one vector)
! **********************************************************************
!
!                       I/O
!
!       Depending on icheli the following integrals are computed
!       2-5:  The input is a vector of the structure of the solution
!             i.e. 1 unknown per point
!         6:  The input is a vector of the structure of the derivative
!             vector, i.e. ndim unknowns per point
!          1:   / f(x) dx
!          2:   / f(x) u(x) dx
!          3-5: / f(x) du(x)/dx  dx   i=icheld-2
!                              i
!          6: / f(x) du      (x) dx
!                      jdegfd
!          7: / dx (i.e. compute the volume)
!          8: / f(x) u^2(x) dx
!        11-16:  See 1-6, however, the integration is performed in the
!                positive x-direction only
!        21-26:  See 1-6, however, the integration is performed in the
!                positive y-direction only
!        31-36:  See 1-6, however, the integration is performed in the
!                positive z-direction only
!        Each integral is approximated by a numerical rule of the shape:
!        / int(x) dx = sum  w  int(x )
!        e              k    k      k
!        the factors w  form the weights for the numerical integration
!                     k
!        In the axi-symmetric case they include the factor 2 pi r
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
      call eropen ( 'sepeli705' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepeli705'
         write(irefwr,*) 'first, icheli', E%first, icheli
  1      format ( a, 1x, (10i6) )
  2      format ( a, 1x, (5d12.4) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      sepeli705 = 0d0

!     --- Define element independent quantities
!         only for first element

      if ( E%first==0 ) call sepeli0705

!     --- Fill basis functions

      call sepelbasefunc

      if ( icheli/=7 ) then

!     --- Fill old vectors in quadrature points if necessary

         debug=.true.
         if ( debug ) &
            write(irefwr,1) 'E%compold, iseqin', E%compold, iseqin
         if ( E%compold>0 ) then

!        --- All old vectors are required

            call sepeloldvecs

         else

!        --- Only one vector required

            call sepeloldvecsk ( iseqin )

         end if  ! ( E%compold>0 )
         uprev => oldsols(iseqprevsol)%locvecs(:,1:nunkp)

         if ( debug ) call prinrl ( uprev, B%n, 'uprev' )

      end if  ! ( icheli/=7 )

!     --- Fill all coefficients

      if ( E%isupcfvar>=E%isubcfvar ) call sepelfillvar

!     --- Compute derivatives

      sepeli705 = sepeli705int ( )

!     --- Remove basis functions and coefficients in case of last call

      if ( debug ) write(irefwr,1) 'last', E%last
      if ( E%last==1 .and. icheli/=7 ) call sepeldealloc

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,2) 'sepeli705', sepeli705
         write(irefwr,*) 'End sepeli705'

      end if  ! ( debug )
      call erclos ( 'sepeli705' )

      end function sepeli705
