      subroutine sepfillmatstr1 ( iincommt )
! ======================================================================
!
!        programmer    Guus Segal
!        version  2.0  date 06-01-2016 Use sepmodulematstr
!        version  1.0  date 27-12-2015
!
!   copyright (c) 2015-2016  "Ingenieursbureau SEPRA"
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
!     Fill parameters from iincommt in module sepmodulematstr
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
      use sepmoduleintmat
      use sepmodulemain
      use sepmodulecpack
      use sepmodulematstr
      use sepmodulekprob
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer, target, intent(in) :: iincommt(*)

!     iincommt       i    Integer input array containing information concerning
!                         the structure of the large matrix
!                         See subroutine read00
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer :: istart

!     istart         First position used at a part
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     PRININ         print 1d integer array
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!      94   jmetod has incorrect value
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
      call eropen ( 'sepfillmatstr1' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepfillmatstr1'
         call prinin ( iincommt, 14, 'iincommt' )
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      iprob  = iincommt(5)
      iseqkp = kprob
      kprobpart => kp(iseqkp)%kprobsub(iprob)
      kprobloc => kprobpart%kprobloc

      jmethod = iincommt(2)
      if ( jmethod==9 .or. jmethod==51 ) jmethod = 6
      if ( jmethod==12 ) jmethod = 8
      iprint  = iincommt(3)
      typecoupl = 0
      extrafillin = iincommt(6)
      iphys1 = iincommt(7)
      iphys2 = min(iincommt(8),nphys)
      numdec = iincommt(10)   ! number of decoupled unknowns
      scheme = iincommt(15)   ! storage/renumbering scheme
      nested = .false.
      write(6,*) 'sepfillmatstr1: ',scheme,scheme,jmethod,kprobpart%iinputloc(14)
      if ( scheme==0 .or. scheme==20 ) then
         nested = jmethod==1 .or. jmethod==2 .and. kprobpart%iinputloc(14)==0
      else if ( jmethod<=2 ) then
         nested = scheme==2 .or. scheme==3
      end if  ! ( scheme==0 )
      if ( numdec>0 ) then
         istart = iincommt(11)   ! starting address information
         decoupled => iincommt(istart:istart+numdec-1)
      end if  ! ( numdec>0 )
      write(6,*) 'sepfillmatstr1: ',nested
      call instop

      fillbouncond = iincommt(12)==0 .and. nbound>0 .or. npboun>0

      if ( parallel .and. iincommt(14)>0 ) then

!     --- We need an adaptation in case of dirichlet-dirichlet coupling

         typecoupl = iincommt(14)

      end if  ! ( parallel .and. iincommt(14)>0 )
      if ( parallel .and. jmethod>=55 .and. jmethod<=56 ) typecoupl = 1 ! petsc

      if ( debug ) write(irefwr,1) 'typecoupl', typecoupl

      if ( jmethod<1 .or. &
           jmethod>=57 .and. jmethod<=59 .or. &
           jmethod>=17 .and. jmethod<=19 .or. &
           jmethod>=28 .and. jmethod<=49 ) then

!     --- jmethod has wrong value

         call errint ( jmethod, 1 )
         call errsub ( 94, 1, 0, 0 )
         go to 1000

      end if  ! ( jmethod<1 .or. jmethod>56 .or. ... )

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End sepfillmatstr1'

      end if  ! ( debug )
      call erclos ( 'sepfillmatstr1' )

      end subroutine sepfillmatstr1
