      subroutine sepcopyvc ( iold, inew )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 28-10-2015
!
!   copyright (c) 2015-2015  "Ingenieursbureau SEPRA"
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
!     copies the contents of array array iold in inew
!
! **********************************************************************
!
!                       KEYWORDS
!
!     copy
!     vector
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulemain
      use sepmodulekprob
      use sepmodulevecs
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer, intent(in) :: iold
      integer, intent(inout) :: inew

!     inew      i/o   new array (gets the same structure as iold)
!     iold       i    old array (must be a SEPRAN array)
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer :: itypvc, fstblk, i, leng, iprobsav, iseqkpsav
      character (len=16) :: namarr1

!     fstblk         First blank found in name
!     i              Counting variable
!     itypvc         Indication of the type of solution vector.
!                    110: solution vector
!                    115, 116: vectors of special structure
!     leng           size of array
!     namarr1        name of array
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     SEPCRSOLINFO   Create a new vector in ks
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!      13   array iold has not been filled.
! **********************************************************************
!
!                       PSEUDO CODE
!
!   trivial
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen ( 'sepcopyvc' )
      itypvc = ks(iold)%typevec
      if ( itypvc/=110 .and. itypvc/=115 .and. itypvc/=116 .and. &
           itypvc/=126 .and. itypvc/=127 .and. itypvc/=129 .and. &
           itypvc/=801 ) then
         call errint ( itypvc, 1 )
         call errsub ( 13, 1, 0, 0 )
      end if
      if ( ierror/=0 ) go to 1000

      iprobsav = iprob
      iseqkpsav = iseqkp

      namarr1 = namesol(iold)
      fstblk = index(namarr1,' ')
      if ( fstblk>1 .and. fstblk<15 ) then

!     --- Still space left at end of name

         namarr1(fstblk:16) = '_cp'

      else if ( fstblk==1 ) then

!     --- Array has not yet a name

         namarr1 = 'inew'

      end if

      if ( itypvc==129 ) then

!     --- Vector of type 129, copy both integer as real part

         ivectin => ks(iold)%isol
         leng = size(ivectin)
         call sepcrsolinfo ( inew, 0, leng )
         ks(inew)%isol(1:leng) = ivectin(1:leng)
         usoltemp => ks(iold)%sol
         leng = size(usoltemp)
         allocate ( ks(inew)%sol(leng) )
         ks(inew)%sol(1:leng) = usoltemp(1:leng)
         ks(inew)%typevec = itypvc
         ks(inew)%nusol = ks(iold)%nusol
         ks(inew)%ivec = ks(iold)%ivec

      else

!     --- Standard case

         usolt => ks(iold)%sol
         iseqkp = ks(iold)%iseqkp
         usoltemp => ks(iold)%sol
         leng = size(usoltemp)
         call sepcrsolinfo ( inew, leng, 0 )
         !ks(inew)%sol(1:leng) = usoltemp(1:leng)
         do i=1,leng
            ks(inew)%sol(i) = usoltemp(i)
         enddo
 
         ks(inew)%typevec = itypvc
         ks(inew)%nusol = ks(iold)%nusol
         ks(inew)%ivec = ks(iold)%ivec
         usoltemp => ks(inew)%sol

      end if  ! ( itypvc==129 )
      ks(inew)%n = ks(iold)%n
      ks(inew)%jcart = ks(iold)%jcart
      ks(inew)%nphys = ks(iold)%nphys
      ks(inew)%iprob = ks(iold)%iprob
      ks(inew)%compl = ks(iold)%compl
      ks(inew)%iseqkp = iseqkpsav

      iprob = iprobsav
      iseqkp = iseqkpsav

1000  call erclos ( 'sepcopyvc' )
      end subroutine sepcopyvc
