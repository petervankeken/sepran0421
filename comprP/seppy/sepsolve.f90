      subroutine sepsolve ( A, usol, rhsd )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.1  date 23-02-2016 New calls to several subroutines
!        version  1.0  date 02-04-2015
!
!   copyright (c) 2015-2016  "Ingenieursbureau SEPRA"
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
!     Solve linear system of equations
!
! **********************************************************************
!
!                       KEYWORDS
!
!     linear_solver
!     solve
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulesolve
      use sepmodulematr
      use sepmodulemat
      use sepmodulekprob
      implicit none

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      type (rmatrix) :: A
      double precision :: usol(A%neq), rhsd(A%neq)

!     A             i/o   Structure of type matrix containing the necessary
!                         information of the matrix
!     rhsd           i    right-hand-side vector
!     usol           o    solution vector
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer method, iseqsolsav

!     iseqsolsav     Help parameter to save the value of iseqsol
!     method         defines type of linear solver
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     SEPAL          Solve system of linear equations of special form using the
!                    Augmented Lagrange appoach
!     SEPBLOCKTRIA   Solve a system of linear equations of special form using
!                    the Block Triangular preconditioner
!     SEPKRYLOV      Solve system of linear equations by a krylov method
!     SEPMUMPS       Solve a linear system of equation by MUMPS
!     SEPOVERRL      Solve system of linear equations by the overrelaxation
!                    method
!     SEPPROFILE     Solve system of real linear equations by
!                    Gaussian elimination
!     SEPSCHUR       Solve system of linear equations of special form using
!                    schur type methods
!     SEPSIMPLE      Solve system of linear equations of special form using the
!                    SIMPLE method
!     SEPSOLVENEST   Solve system of real linear equations by
!                    Gaussian elimination
!                    The nested dissection approach is applied
!     SEPSOLVEPETSC  Solve system of equations by a multigrid method
!                    using the petsc library as interface
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!    2987   Iterative solver not available for real matrices
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
      call eropen ( 'sepsolve' )
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepsolve'
         write(irefwr,1) 'method', A%typesolv
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      method = A%typesolv     ! type of solution method
      write(6,*) 'sepsolvel method: ',method

      select case  (method)

      case (0)

!     --- method =0, profile method

         call sepprofile ( A, usol, rhsd )

      case (1:4, 7, 9, 10, 17, 31 )

!     --- method = 1-4, 7, 9, 10, 17, standard krylov method

         call sepkrylov ( A, usol, rhsd )
         inpsol(1,1) = A%numiter

      case (5)

!     --- method =5, overrelaxation

         call sepoverrl ( A, usol, rhsd )

      case (6, 8, 12, 13, 14, 15)

!     --- method = 6, 8, 12-15, block type method for (Navier)-Stokes

         call sepsimple ( A, usol, rhsd )

      case (11)

!     --- method =11, AL type method for (Navier)-Stokes

         call sepal ( A, usol, rhsd )

      case (16)

!     --- method = 16, schur type method for (Navier)-Stokes

         call sepschur ( A, usol, rhsd )

      case (18)

!     --- method = 18, BLOCK_TRIANGULAR preconditioner for (Navier)-Stokes

         call sepblocktria ( A, usol, rhsd )

      case (19)

!     --- method = 19, MUMPS solver

         call sepmumps ( A, usol, rhsd )

      case (20)

!     --- method = 20, Petsc solver

         call sepsolvepetsc ( A, usol, rhsd )
         inpsol(1,1) = niter

      case (21)

!     --- method = 21, Direct solver using nested dissection

         call sepsolvenest ( A, usol, rhsd )

      case default

!     --- Solver not available

         call errint ( method, 1 )
         call errsub ( 2987, 1, 0, 0 )

      end select  ! case  (method)

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End sepsolve'

      end if  ! ( debug )
      call erclos ( 'sepsolve' )

      end
