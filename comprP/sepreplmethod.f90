      subroutine sepreplmethod
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 25-02-2020
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
!     Translate jmethod to logicals in type intmatrix
!
! **********************************************************************
!
!                       KEYWORDS
!
!     matrix_structure
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulematstr
      use sepmoduleintmat
      use sepmodulekprob
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

!     IM             Refers to intmt(intmat)
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERTRACE        Prints trace back of subroutines
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
      call eropen ( 'sepreplmethod' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepreplmethod'
         call ertrace
         write(irefwr,1) 'jmethod, iseqim', jmethod, iseqim
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      IM => intmt(iseqim)
      ! restore to 1020 version. Current version sets nested to true for profile methods
      IM%nested = jmethod==60 ! nested
      IM%direct = jmethod>=1 .and. jmethod<=4
      IM%full = jmethod>=20 .and. jmethod<=23
      IM%lapshift = jmethod>=24 .and. jmethod<=27
      IM%mumps = jmethod>=52 .and. jmethod<=54
      IM%posdev = jmethod==54
      IM%nosepmat = jmethod==56
      IM%petsc = jmethod==55 .or. jmethod==56
      IM%simple = jmethod>=13 .and. jmethod<=16
      IM%symmetric = jmethod==1 .or. jmethod==3 .or. jmethod==5 .or. &
                     jmethod==7 .or. jmethod==20 .or. jmethod==22 .or. &
                     jmethod==53 .or. jmethod==54
      IM%Qsymmetric = jmethod==13 .or. jmethod==15
      IM%Gsymmetric = jmethod==13 .or. jmethod==14
      if ( IM%simple ) jmethod = 5
      if ( IM%lapshift ) jmethod = jmethod-19
      IM%complexm = jmethod==3 .or. jmethod==4 .or. jmethod==7  .or. &
                    jmethod==8 .or. jmethod==22 .or. jmethod==23

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,1) 'jmethod, iseqim', jmethod, iseqim
         write(irefwr,*) 'IM%simple ', IM%simple
         write(irefwr,*) 'End sepreplmethod'

      end if  ! ( debug )
      call erclos ( 'sepreplmethod' )

      end subroutine sepreplmethod
