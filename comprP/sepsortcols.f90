      subroutine sepsortcols
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 22-03-2020
!
!   copyright (c) 2020-2020  "Ingenieursbureau SEPRA"
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
!     Sort column numbers per row of internal matrix
!
! **********************************************************************
!
!                       KEYWORDS
!
!     sort
!     matrix_structure
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulekprob
      use sepmodulematstr
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
      ! PvK: make icheck and istore allocatable
      integer,allocatable :: icheck(:),istore(:)
      integer i, ifirst, isum, ncols, ilast, j !, icheck(nusol), istore(nusol)

!     i              Loop variable
!     icheck         work array for sorting
!     ifirst         Lower bound of row loop
!     ilast          upper bound of row loop
!     istore         work array for sorting
!     isum           Integer variable to compute a sum
!     j              Loop variable
!     ncols          Number of columns in a row
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     CHSORT         Sort integer array for increasing sequence
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
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
      allocate(icheck(nusol),istore(3*nusol))
      icheck=0
      istore=0

      call eropen ( 'sepsortcols' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepsortcols'
         write(irefwr,*) 'associated(irowsii), size(irowsii) ', &
                          associated(irowsii), size(irowsii)

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      isum = irowsii(1)
      ifirst = 2
      ilast = nrusol+1

      do i = ifirst, ilast
         ncols = irowsii(i) - isum
         do j = 1, ncols
            istore(j) = icolsii(isum+j)
            icheck(j)  = j
         end do
         if ( ncols>0 ) call chsort ( istore, icheck, ncols )
         do j = 1, ncols
            icolsii(isum+j) = istore(icheck(j))
         end do
         isum = irowsii(i)
      end do
      if (allocated(icheck)) deallocate(icheck)
      if (allocated(istore)) deallocate(istore)


1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End sepsortcols'

      end if  ! ( debug )
      call erclos ( 'sepsortcols' )

      end subroutine sepsortcols
