      subroutine seplinsolreal ( A, u, f )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.1  date 13-12-2017 New storage of intmat
!        version  1.0  date 05-03-2016
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
!     Solve system of linear equations by either a direct or an iterative
!     method
!     Real matrix
! **********************************************************************
!
!                       KEYWORDS
!
!     linear
!     linear_solver
!     solve
!     real
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulevecs
      use sepmodulematr
      use sepmoduleintmat
      use sepmoduleproj
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      type (rmatrix) :: A
      double precision :: u(A%neq)
      double precision, intent(in) :: f(A%neq)

!     A              i    Structure of type matrix containing the necessary
!                         information of the matrix
!     f         i    right-hand-side vector
!     u         o    solution vector
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
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     PRINRL         Print 1d real vector
!     SEPPARGLOB     Send information to neighboring processors in order to give
!                    all local arrays the global value in the common points
!     SEPSOLVE       Solve system of linear equations
!     SEPSOLVEPAR    Solve system of linear equations in a parallel environment
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!    2140   ISEQ_RES_EXACT has been given, but not ISEQ_EXACT
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
      call eropen ( 'seplinsolreal' )
      debug = .false. .and. ioutp>=0
      debug=.true.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from seplinsolreal'
         write(irefwr,1) 'typesolv, type_prec',  A%typesolv, A%type_prec
         write(irefwr,1) 'typecoupl', typecoupl
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      debug=.false.
      if ( ierror/=0 ) go to 1000

      if ( parallel .and. typecoupl==0 ) then

!     --- parallel case, make right-hand side global

         length = size(A%sendinfo)
         call sepparglob ( A, f, A%neq ) ! make global rhsd

      end if  ! ( parallel )

!     --- Store matrix in type matrix

      if ( A%typesolv==0 .or. A%typesolv==21 ) &
         A%decomp = realmt(iseqmt)%decomp
      copyprec = .true.   ! Default

      if ( debug ) write(irefwr,*) 'A%decomp ', A%decomp

      if ( iseqres_exact>0 ) then

!     --- iseqres_exact>0, compute residual of exact solution

         if ( iseq_exact<=0 ) then

!        --- iseq_res_exact given, but not iseq_exact

            call errint ( iseqres_exact, 1 )
            call errsub ( 2140, 1, 0, 0 )

         else

!        --- Compute residual

            resusolex = residual ( A, usolex, f )

         end if  ! ( iseqexact<=0 )

      end if  ! ( iseqres_exact>0 )

      if ( iseq_start_res>0 ) startres = residual ( A, u, f )

      if ( parallel .and. typeparsolv>1 ) then

!     --- Parallel case, new type of krylov solver

         call sepsolvepar ( A, u, f )

      else

!     --- Serial case, new type of krylov solver

         if ( A%typesolv==19 ) A%newmumps = intmt(iseqim)%signmethod
         call sepsolve ( A, u, f )

         if ( A%typesolv==19 ) intmt(iseqim)%signmethod = A%newmumps
         if ( A%typesolv==0 .or. A%typesolv==21 ) &
              realmt(iseqmt)%decomp = A%decomp

      end if  ! ( new )

      nullify(A%projrows)
      nullify(A%projcols)
      nullify(A%projU)
      nullify(A%projE)
      nullify(A%projAu)

1000  if ( debug ) then

!     --- Debug information

         call prinrl ( u, A%neq, 'u' )
         write(irefwr,*) 'End seplinsolreal'

      end if  ! ( debug )
      call erclos ( 'seplinsolreal' )

      end subroutine seplinsolreal
