      subroutine sepcopytoprec9 ( A )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 20-02-2020
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
!     Copy non-symmetric matrix into preconditioning matrix
!     Row-compact storage
! **********************************************************************
!
!                       KEYWORDS
!
!     preconditioning
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulesolve
      use sepmodulemat
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      type (rmatrix) :: A
      integer :: i

!     A             i/o   Structure of type matrix containing the necessary
!                         information of the matrix
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
!     INSTOP         Stop the program
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
      call eropen ( 'sepcopytoprec9' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepcopytoprec9'

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      if ( iseqpre/=iseqintmat ) then

!     --- the structure of the preconditioning
!         matrix and the standard matrix is different

    !     call preccopy ( A%neq, A%rows, A%Mrows, A%L, A%U, &
    !                     A%ML, A%MU, A%cols, A%Mcols, A%MD, A%D, 6 )
      print *, 'option preccopy not yet available'
      call ertrace
      call instop


      else

!     --- both matrices have the same structure
!         Copy matrix into premat

         if ( size(A%MD)==A%neq ) then
            A%MD = A%D(A%diag)
         else
            do i=1,size(A%MD)
               A%MD(i) = A%D(i)
            enddo
         end if  ! ( size(A%MD)==A%neq )
         A%Mrows => A%rows
         A%Mcols => A%cols
         A%Mdiag => A%diag

      end if  ! ( inpsol(40)>0 )

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End sepcopytoprec9'

      end if  ! ( debug )
      call erclos ( 'sepcopytoprec9' )

      end subroutine sepcopytoprec9
