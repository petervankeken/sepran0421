      subroutine torenum
! ======================================================================
!
!        programmer    Guus Segal
!        version  2.3  date 06-01-2016 Long integers
!        version  2.2  date 27-10-2010 New call to toreordr
!        version  2.1  date 04-02-2005 seppointer --> seppointer
!
!   copyright (c) 1998-2016  "Ingenieursbureau SEPRA"
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
!    Renumber right-hand side from numbering as used by matrix to numbering
!    used in standard solution arrays, i.e. without renumbering
!    This routine is part of subroutine BUILD
! **********************************************************************
!
!                       KEYWORDS
!
!     renumbering
!     right_hand_side
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulevecs
      use sepmodulebuild
      use sepmodulekprob
      use sepmoduleelem
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer :: i,j
      real(kind=8),allocatable :: rhsdbuf(:)

!     i              Counting variable
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS          Resets old name of previous subroutine of higher level
!     EROPEN          Produces concatenated name of local subroutine
!     ERTRACE        Prints trace back of subroutines
!     TOREORDC        Map array of solution type from matrix sequence to vector
!                     sequence. Complex version
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
      call eropen ( 'torenum' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from torenum'
         call ertrace
         write(irefwr,1) 'numrhs, nusol', numrhs, nusol
         write(irefwr,*) 'compl ', compl
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

!     --- Renumber rhsd

      allocate(rhsdbuf(nusol))
      do i = 1, numrhs

         rhsd => rhsds(i)%rh
         if ( compl ) then

!        --- complex case

            call toreordc ( kprobr, rhsd, nusol )

         else

!        --- real case

            do j=1,nusol
               rhsdbuf(j) = rhsd(kprobr(j))
            enddo
            do j=1,nusol
               rhsd(j)=rhsdbuf(j)
            enddo

         end if

      end do
      if (allocated(rhsdbuf)) deallocate(rhsdbuf)

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End torenum'

      end if  ! ( debug )
      call erclos ( 'torenum' )

      end subroutine torenum
