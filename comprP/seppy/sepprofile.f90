      subroutine sepprofile ( A, usol, rhsd )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 25-03-2016
!
!   copyright (c) 2016-2016  "Ingenieursbureau SEPRA"
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
!     Solve system of real linear equations by Gaussian elimination
!
! **********************************************************************
!
!                       KEYWORDS
!
!     linear_solver
!     real
!     direct_method
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulemat
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      type (rmatrix) :: A
      double precision :: usol(A%neq)
      double precision, intent(in) :: rhsd(A%neq)
      real(kind=8) :: get_resmem

!     A             i/o   Matrix A
!     rhsd           i    right-hand-side vector
!     usol           o    solution vector
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      logical :: decomp

!     decomp         Logical parameter indicating if the decomposition of
!                    the matrix must be computed (true) or not (false)
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     SEPPROFILEDEC  Compute lU decomposition for system of real linear
!                    equations
!     SEPPROFILESOL  Solve system of real linear equations by Gaussian
!                    elimination
!                    The matrix is supposed to be decomposed before
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
      call eropen ( 'sepprofile' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepprofile'

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

!     write(6,*) '     sepprofile: sepprofiledec: ',get_resmem()
      if ( A%decomp ) call sepprofiledec ( A )  ! Compute decomposition

!     --- Next solve system of equations

!     write(6,*) '     sepprofile: sepprofilesol: ',get_resmem()
      call sepprofilesol ( A, usol, rhsd )
!     write(6,*) '     sepprofile: end          : ',get_resmem()

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End sepprofile'

      end if  ! ( debug )
      call erclos ( 'sepprofile' )

      end subroutine sepprofile
