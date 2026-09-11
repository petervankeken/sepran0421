      subroutine sepaddmat ( matr1, matr2, &
                             alpha1, alpha2, beta1, beta2 )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 18-01-2018
!
!   copyright (c) 1999-2018  "Ingenieursbureau SEPRA"
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
!    add alpha times matrix amat1 to beta times matrix amat2 and
!    store in amat1
!    amat1 must be a large matrix
!    amat2 may be either a large matrix or a diagonal matrix
! **********************************************************************
!
!                       KEYWORDS
!
!     matrix
!     multiplication
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulematr
      use sepmodulemain
      use sepmoduleintmat
      use sepmodulematstr
      use sepmodulekprob
      use sepmodulecactl
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer :: matr1, matr2
      double precision :: alpha1, alpha2, beta1, beta2

!     alpha1  i    multiplication factor:  alpha = alpha1 + i alpha2
!     alpha2  i    multiplication factor:  alpha = alpha1 + i alpha2
!     beta1   i    multiplication factor:  beta = beta1 + i beta2
!     beta2   i    multiplication factor:  beta = beta1 + i beta2
!     matr1   i    integer array containing information concerning the
!                  large matrix amat1
!     matr2   i    integer array containing information concerning the
!                  matrix amat2
!                  amat2 may be a large matrix or a diagonal matrix
!                  stored as right-hand-side vector
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer :: k, minnum, nloc, lengmt, maxnum, &
                 iseqmat1, iseqmat2
      double precision :: timebefore(2)
      logical :: diag, compl, bouncond
      double precision, pointer :: matrix2(:)
      type (intmatrix), pointer :: IM
      integer :: i

!     IM             Refers to intmt(intmat)
!     bouncond       If true boundary conditions have been read
!     compl          indication whether the arithmetic is real (false)
!                    or complex (true)
!     diag           Diagonal of the large matrix
!     iseqmat1       Sequence number of matr1 in realmt
!     iseqmat2       Sequence number of matr2 in realmt
!     k              Counting variable
!     lengmt         Length of matrix
!     matrix2        Second matrix to be added
!     maxnum         Maximum row number for part of decomposition
!     minnum         Minimum row number for part of decomposition
!     nloc           number of free degrees of freedom
!     timebefore     Time at entrance of subroutine
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ER0019C        Add diagonal matrix to large matrix (complex case)
!     ER0020C        Add two large matrices (complex case)
!     ERCLMN         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERSETTIME      Get the CPU time at start of subroutine
!     PRININ         print 1d integer array
!     SEPFILLKPROB   Fill parameters in module sepmodulekprob
!     SEPFILLMATSTR  Fill parameters in module sepmodulematstr
!     SEPPRINTMAT    Print matrix
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
      call eropen ( 'sepaddmat' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepaddmat'
         write(irefwr,2) 'alpha1, alpha2, beta1, beta2', &
                          alpha1, alpha2, beta1, beta2
         call prinin ( matr1, 1, 'matr1' )
         call prinin ( matr2, 1, 'matr2' )
  1      format ( a, 1x, (10i6) )
  2      format ( a, 1x, (5d12.4) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

!     --- tests

      iseqmat1 = matr1
      intmat = realmt(iseqmat1)%intmatref

      IM => intmt(intmat)

      call sepfillmatstr
      iprob = realmt(iseqmat1)%iprob
      iseqkp = kprob
      kprobpart => kp(iseqkp)%kprobsub(iprob)
      kprobloc => kprobpart%kprobloc
      kprobf => kprobpart%kprobf
      indprf = kprobloc(19)
      call sepfillkprob

      if ( itime>0 ) call ersettime ( timebefore )
      nusolcac=nusol
      nloc=nusolcac-nbound-npboun

      if ( debug ) then
         write(irefwr,1) 'nusol, nbound, n', nusolcac, nbound, nloc
         write(irefwr,1) 'indprh', indprh
      end if  ! ( debug )

      S_ii => realmt(iseqmat1)%S_ii
      iseqmat2 = matr2
      matrix2 => realmt(iseqmat2)%S_ii

      diag = realmt(iseqmat2)%diagonal
      compl = realmt(iseqmat1)%complexmat
      bouncond = associated(realmt(iseqmat1)%S_ip)
      if ( bouncond ) S_ip => realmt(iseqmat1)%S_ip

      if ( debug ) write(irefwr,*) 'diag', diag

      minnum=1

      lengmt = size ( S_ii )
      maxnum=nloc

      if ( ierror/=0 ) go to 1000

      if ( diag ) then

         if ( compl ) then

!        --- Complex first matrix

            call er0019c ( matrix2, S_ii, &
                           IM%rowsii, lengmt, alpha1, &
                           beta1, jmethod, minnum, maxnum )
         else

!        --- Both matrices real

            S_ii = alpha1*S_ii
            S_ii(IM%diag) = S_ii(IM%diag)+beta1*matrix2

         end if

      else

         if ( compl ) then

!        --- Complex first matrix

            call er0020c ( S_ii, matrix2, lengmt, &
                           alpha1, beta1 )
         else

!        --- Both matrices real

            ! PvK write out do loop to avoid segfault
            do i=1,lengmt
               S_ii(i) = alpha1*S_ii(i)+beta1*matrix2(i)
            enddo

         end if
      end if
      minnum=maxnum+1

      realmt(iseqmat1)%decomp = .true.

1000  if ( debug ) then

!     --- Debug information

         call sepprintmat ( matr1 )
         write(irefwr,*) 'End sepaddmat'

      end if  ! ( debug )
      call erclmn ( 'sepaddmat', timebefore )

      end subroutine sepaddmat
