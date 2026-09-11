      subroutine prinlg ( namear, ivalue, array, iarray )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.8  date 19-12-2014 Print buffer length if ioutp>=2
!        version  1.7  date 04-06-2007 Fortran 90
!        version  1.6  date 06-08-2002 Extension of maximal length to 1000000
!
!   copyright (c) 1990-2014  "Ingenieursbureau SEPRA"
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
!   subroutine prinlg checks the length of by the user declared array
!   lengths and depending on ioutp prints the length if has been extended
! **********************************************************************
!
!                       KEYWORDS
!
!     array_length
!     check
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulecomio
      use sepmoduleerr
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      character * (*) namear
      integer ivalue, iarray(*)
      double precision array(*)

!     array   i   Real array to be checked
!     iarray  i   Integer array to be checked
!     ivalue  i   The sign of ivalue indicates whether the real array "ARRAY"
!                 must be checked (ivalue<0) or the integer array "IARRAY"
!                 |ivalue| gives the computed length of the array
!     namear  i   Name of the array to be checked
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer jvalue, ih

!    ih       Help variable to convert from real to integer
!    jvalue   Absolute value of ivalue
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS  Resets old name of previous subroutine of higher level
!     EROPEN  Produces concatenated name of local subroutine
!     ERRCHR  Put character in error message
!     ERRINT  Put integer in error message
!     ERRSUB  Error messages
!     INSTOP  Stop the program
! **********************************************************************
!
!                       I/O
!
!    none
! **********************************************************************
!
!                       ERROR MESSAGES
!
!       2   Array is too short
!    1826
! **********************************************************************
!
!                       PSEUDO CODE
!
!   if ( ivalue>0 ) then
!      if ( ivalue>iarray(1) ) then
!         error
!      if ( ivalue>iarray(3) ) then
!         iarray(3) := ivalue
!         print new length of array
!   else if ( ivalue<0 ) then
!      if ( -ivalue>array(1) ) then
!         error
!      if ( -ivalue>array(3) ) then
!         array(3) := ivalue
!         print new length of array
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
1     format(' computed array length of array ', a, ' is equal to' , &
               i9/ ' in subroutine ', a )

      jvalue = abs(ivalue)

      if ( ivalue>0 ) then

!     --- ivalue>0, integer array

         if ( iarray(1)<=0 .or. iarray(1)>1d6 ) then

!        --- iarray(1) is incorrect

            call eropen ( 'prinlg' )
            call errint ( iarray(1), 1 )
            call errchr ( namear, 1 )
            call errsub ( 1826, 1, 0, 1 )
            call erclos ( 'prinlg' )

         end if

         if ( iarray(3)>iarray(1) ) iarray(3) = iarray(1)

         if ( jvalue>iarray(3) ) then

!        --- New length longer than preceding one

            if ( iarray(1)<jvalue ) then

!           --- error 2: array length too small

               call eropen ( 'prinlg' )
               call errint ( jvalue,1 )
               call errint ( iarray(1), 2 )
               call errchr ( namear, 1 )
               call errsub ( 2, 2, 0, 1 )
               call erclos ( 'prinlg' )

            end if

            iarray(3) = jvalue

            if ( ioutp>=2 ) write(irefwr,1) namear, jvalue, &
               namsub(levsub)

         end if

      else if ( ivalue<0 ) then

!     --- ivalue<0, real array

         if ( array(1)<=0 .or. array(1)>1d9) then

!        --- array(1) may be incorrect

            call eropen ( 'prinlg' )
            ih = nint(array(1))
            call errint ( ih, 1 )
            call errchr ( namear, 1 )
            call errsub ( 1826, 1, 0, 1 )
            call erclos ( 'prinlg' )

         end if

         if ( array(3)>array(1) ) array(3) = array(1)

         if ( jvalue>array(3) ) then

!        --- New length longer than preceding one

            if ( array(1)<jvalue ) then

!           --- error 2: array length too small

               call eropen ( 'prinlg' )
               call errint ( jvalue,1 )
               ih = nint(array(1))
               call errint ( ih, 2 )
               call errchr ( namear, 1 )
               call errsub ( 2, 2, 0, 1 )
               call erclos ( 'prinlg' )

            end if

            array(3) = jvalue

            if ( ioutp>=2 ) write(irefwr,1) namear, jvalue, &
               namsub(levsub)

         end if

      end if

      if ( ierror/=0 ) call instop

      end
