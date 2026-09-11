      subroutine sepsolveinput ( iread, inpsl1, rinsl1 )
! ======================================================================
!
!        programmer    Guus Segal
!        version  2.0  date 23-02-2016 New parameters list
!        version  1.6  date 08-06-2014 Adaption for mumps
!        version  1.5  date 10-04-2014 Extension with mumps
!
!   copyright (c) 2003-2016  "Ingenieursbureau SEPRA"
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
!     Fill arrays inpsol and rinsol for subroutine solvelbf
!
! **********************************************************************
!
!                       KEYWORDS
!
!     linear
!     linear_solver
!     solve
!     input
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulematstr
      use sepmodulekprob
      use sepmodulecons
      use sepmoduleinput
      use sepmodulesolve
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer iread, inpsl1(*)
      double precision rinsl1(*)

!     inpsl1         i    Integer input array for subroutine solvel
!                         See solvelbf for a description
!     iread          i    Defines how the input must be read
!                         Possible values:
!                         -1:  No input is required
!                          0:  All SEPRAN input has been read by
!                              subroutine SEPSTN until
!                              END_OF_SEPRAN_INPUT or
!                              end of file has been found
!                              This input is used
!                          1:  The input is read as described for
!                              EIGENVAL
!     rinsl1         i    Contains information concerning the linear solver
!                         See the Programmers Guide for a description
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERTRACE        Prints trace back of subroutines
!     PRININ         print 1d integer array
!     PRINRL         Print 1d real vector
!     SEPCOPINPSOL   Copy arrays inpsl1 and rinsl1 into inpsol and rinsol
!                    into sepmodulesolve
!     SEPCOPINPSOLIN Copy input for subroutine seplinsol from sepmoduleinput
!     SEPREADINPSOL  Read arrays inpsol and rinsol and store in sepmodulesolve
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
      call eropen ( 'sepsolveinput' )
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepsolveinput'
         call ertrace
         write(irefwr,1) 'iread, ninpsol, nrinsol', iread, ninpsol, nrinsol
         write(irefwr,1) 'nprob', nprob
  1      format ( a, 1x, (10i6) )
         if ( iread==-1 ) then
            call prinin ( inpsl1, inpsl1(1), 'inpsl1' )
            call prinrl ( rinsl1, 20, 'rinsl1' )
         end if

      end if
      if ( ierror/=0 ) go to 1000

      if ( iread==-1 ) then

!     --- iread = -1, fill inpsol from inpsl1

         call sepcopinpsol ( inpsl1, rinsl1 )

      else if ( iread==0 .or. iread==-2 ) then

!     --- iread = 0 or -2 fill inpsol by reading from the standard input file

         call sepreadinpsol ( inpsl1, rinsl1, iread )

      else

!     --- iread>0, fill inpsol by copying from IBUFFR

         call sepcopinpsolin ( iread )

      end if

      debug=.true.
1000  if ( debug ) then
         call prinin ( inpsol, ninpsol, 'inpsol' )
         call prinrl ( rinsol, 20, 'rinsol' )
         write(irefwr,1) 'nprob', nprob
         write(irefwr,*) 'End sepsolveinput'
      end if  ! ( debug )
      debug=.false.
      call erclos ( 'sepsolveinput' )

      end subroutine sepsolveinput
