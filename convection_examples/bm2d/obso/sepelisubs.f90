      double precision function sepelisubs ( )
! ======================================================================
!
!        programmer    Guus Segal
!        version  3.0  date 06-01-2016 Use sepmoduleelem
!        version  2.0  date 06-01-2016 Use sepmoduleelem
!        version  1.2  date 11-11-2014 Extension with type 700
!
!   copyright (c) 2010-2016  "Ingenieursbureau SEPRA"
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
!     Compute value of integral for one element
!
! **********************************************************************
!
!                       KEYWORDS
!
!     element
!     element_group
!     integral
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulebuild
      use sepmoduleelm
      use sepmodulebasefn
      use sepmoduleoldvecs
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!

!     sepelisubs     o    Result of integration in element
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer, pointer :: iuser(:), islold(:)
      double precision, pointer :: user(:)
      double precision :: value

!     islold         User input array in which the user puts information
!                    of all preceding solutions
!     iuser          Integer user array to pass user information from
!                    main program to subroutine. See STANDARD PROBLEMS
!     user           Real user array to pass user information from
!                    main program to subroutine. See STANDARD PROBLEMS
!     value          computed element integral
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
      double precision elintsub, eli100, eli113, eli118, eli400, &
                       eli900, eli903, eli325
      double precision :: sepeli705

!     ELI100         Element subroutine for second order elliptic equation (2D)
!     ELI113         Element subroutine for second order elliptic equation (3D)
!     ELI118         Element subroutine for second order elliptic equation (1D)
!     ELI325         Compute element integral for bearing
!     ELI400         Compute element integral for quadratic triangles
!     ELI900         Compute integral for type numbers 410, 900
!     ELI903         Compute integral for type numbers 903
!     ELINTSUB       Element subroutine for user element
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     SEPELDEFAULT   Set default values for parameters in module sepmodulebasefn
!     SEPELI705      Element subroutine for second order elliptic equation
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!      43   Type number not available
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
      call eropen ( 'sepelisubs' )
      debug = .true.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepelisubs'
         write(irefwr,1) 'itype, ifirst', itype, ifirst
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      iuser => E%iuser
      user => E%user
      islold => E%islold

      if ( ifirst==0 ) call sepeldefault

      select case ( itype )

         case ( 1:99 )

!        --- User element

            value = elintsub ( icheli, jdegfd, iuser, user )

         case ( 300 )

            value = eli100 ( coor, iuser, user, index1, islold )

         case ( 150:151, 153, 155, 158:161, 163, 700:703, 705, 800, 803, &
                807)

            value =  sepeli705 ( )

         case ( 325 )

            if ( icheli<40 ) then

               value =  sepeli705 ( )

            else

              value =  eli325 ( coor, index1 )

            end if  ! ( icheli<40 )

         case ( 400, 404, 420, 421 )

!        --- ielrou = 150  elements of type 400, 404 etc.

            value = eli400 ( icheli, jdegfd, coor, iuser, user, vecold, &
                             index1, index2 )

         case ( 900:902, 410 )

            value = eli900 ( coor, iuser, user, index1, &
                             numold, islold, coeffarr )

         case ( 903 )

            value = eli903 ( coor, iuser, user, islold, coeffarr )

         case default

!        --- Type number not available

            call errint ( itype, 1 )
            call errint ( ielgrp, 2 )
            call errsub ( 43, 2, 0, 0)

      end select  ! case ( itype )
      sepelisubs = value

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End sepelisubs'

      end if  ! ( debug )
      call erclos ( 'sepelisubs' )

      end function sepelisubs
