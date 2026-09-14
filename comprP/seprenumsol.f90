      subroutine seprenumsol ( usol, itrans )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 17-03-2016
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
!     Renumber solution vector such that prescribed boundary conditions are
!     stored at the end.
!     Local transformations are carried out if necessary
!     Real case
! **********************************************************************
!
!                       KEYWORDS
!
!     real
!     ordering
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulekprob
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer, intent(in) :: itrans
      double precision, intent(inout) :: usol(nusol)
      integer :: i
      real(kind=8),allocatable :: usolbuf(:)

!     itrans         i    If > 0 local transformations are carried out if
!                         npltra>0
!     usol          i/o   Solution vector to be transformed
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
!     SEPLOCTRNRL    Compute (back-)transformation in case of local
!                    transformations ( real case )
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
      call eropen ( 'seprenumsol' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from seprenumsol'
         write(irefwr,1) 'npltra, itrans, indprh', npltra, itrans, indprh
         write(irefwr,1) 'nusol', nusol
         call ertrace
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

!     --- Transform essential boundary conditions in case of local transform

      if ( npltra>0 .and. itrans>0 ) call seploctrnrl ( 1, usol )
!      if ( indprh>0 ) usol = usol(kprobhinv)  ! reorder usol
      allocate(usolbuf(nusol))
      do i=1,nusol
         usolbuf(i)=usol(kprobhinv(i))
      enddo
      do i=1,nusol
         usol(i)=usolbuf(i)
      enddo
      if (allocated(usolbuf)) deallocate(usolbuf)

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End seprenumsol'

      end if  ! ( debug )
      call erclos ( 'seprenumsol' )

      end subroutine seprenumsol
