      subroutine sepcommatfill
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 31-01-2020
!
!   copyright (c) 2020-2020  "Ingenieursbureau SEPRA"
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
!     Fill information concerning the storage of the large matrix
!     Old subroutine incommatfill
! **********************************************************************
!
!                       KEYWORDS
!
!     matrix_structure
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulekprob
      use sepmoduleintmat
      use sepmodulesimple
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
      type (intmatrix), pointer :: IM
      integer :: i

!     IM             Refers to intmt(intmat)
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     SEPCOMMATREV   Fill arrays irowsnew and icolsnew with upper triangular
!                    part if icolsii and irowsii contain the lower triangular
!                    part and vice versa
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
      call eropen ( 'sepcommatfill' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepcommatfill'

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

!     --- Compute length of array copy parts of matrices to
!         temporary arrays

      IM => intmt(intmat)

      IM%ninterface = ninterface
      IM%jmethod = jmethod
      IM%signmethod = 1

      if ( associated(irowsii) ) then

!     --- S_ii, rows

         if ( IM%symmetric .and. .not. (IM%direct .or. IM%simple) ) then

!        --- Symmetric case, use column-type storage instead or row-type

            length = size(irowsii)
            allocate ( IM%rowsii(size(irowsii)), IM%colsii(size(icolsii)) )
            call sepcommatrev ( IM%rowsii, IM%colsii, irowsii, icolsii, nrusol )
            if ( size(irowsii)>1 ) &
               IM%diag = IM%rowsii(1:size(irowsii)-1)+1

         else

!        --- Standard case

            length = size(irowsii)
            allocate ( IM%rowsii(length) )
            IM%rowsii = irowsii(1:length)
            length = size(icolsii)
            allocate ( IM%colsii(length) )
 
            !this causes segfault IM%colsii = icolsii(1:length)
            do i=1,length
               IM%colsii(i) = icolsii(i)
            enddo

         end if  ! ( IM%symmetric .and. .not. (IM%direct .or. IM%simple) )

      end if  ! ( associated(irowsii) )

      IM%typecoupl = typecoupl

      nu = IM%nu
      np = IM%np

      if ( debug ) write(irefwr,*) 'before deallocate'

      if ( associated(irowsii) ) deallocate ( irowsii, stat = error )
      if ( associated(icolsii) ) deallocate ( icolsii, stat = error )
      if ( allocated(kprobhglobinv) ) deallocate ( kprobhglobinv, stat = error )

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End sepcommatfill'

      end if  ! ( debug )
      call erclos ( 'sepcommatfill' )

      end subroutine sepcommatfill
