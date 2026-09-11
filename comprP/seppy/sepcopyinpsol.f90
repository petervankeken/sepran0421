      subroutine sepcopyinpsol ( A, iseqinp )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 04-07-2015
!
!   copyright (c) 2015-2015  "Ingenieursbureau SEPRA"
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
!     Put information of inpsol and rinsol in matrix structure
!
! **********************************************************************
!
!                       KEYWORDS
!
!     copy
!     matrix
! **********************************************************************
!
!                       MODULES USED
!
      use sepmoduleintmat
      use sepmodulematr
      use sepmodulemat
      use sepmodulekprob
      use sepmodulematstr
      use sepmodulesolve
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer, intent(in) :: iseqinp
      type (rmatrix), intent(inout) :: A

!     A             i/o   Matrix A
!     iseqinp        i    Sequence number of second entry of inpsol/rinsol
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer :: iaccur, iprec
      type (intmatrix), pointer :: IM

!     IM             Refers to intmt(intmat)
!     iaccur         defines how accuracy is stored
!     iprec          Type of preconditioning.
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
      double precision :: getscal

!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERTRACE        Prints trace back of subroutines
!     GETSCAL        Return with value of scalar
!     PRININ         print 1d integer array
!     PRINRL         Print 1d real vector
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
      call eropen ( 'sepcopyinpsol' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepcopyinpsol'
         call ertrace
         write(irefwr,1) 'iseqinp', iseqinp
         call prinin ( inpsol(1,iseqinp), 40, 'inpsol' )
         call prinrl ( rinsol(1,iseqinp), 5, 'rinsol ')
  1      format ( a, 1x, (10i11) )
  2      format ( a, 1x, (5d12.4) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      IM => intmt(iseqim)

      iaccur = inpsol(11,iseqinp)
      nmax = max(inpsol(5,iseqinp),inpsol(43,iseqinp))
      typesolv = inpsol(3,iseqinp)
      write(6,*) 'sepcopyinpsol b: ',typesolv
      if ( IM%nested ) typesolv = 21  !nested dissection
      write(6,*) 'sepcopyinpsol a: ',typesolv
      call instop
      if ( nmax==0 ) then
         if ( parallel ) then
            nmax = 10*nusolglob
         else
            if ( typesolv/=0 ) nmax = 2*A%neq
         end if  ! ( parallel )
      end if  ! ( nmax==0 )
      ipreco = inpsol(4,iseqinp)
      if ( debug ) write(irefwr,1) 'ipreco, typesolv', ipreco, typesolv

      if ( typesolv==0 .or. typesolv==19 ) ipreco = 0

!     --- Select a standard (=default) preconditioner:

      iprec = mod(ipreco,100)  ! Type of preconditioning

      if ( iprec==6 .and. typesolv/=20 ) iprec = 3

!     --- Put new value of iprec in ipreco
!         The old value is multiplied by 100 and added to ipreco

      ipreco = iprec+100*(ipreco/100)

!     --- Set accuracy eps

      A%eps = rinsol(1,iseqinp)
      A%subsolve = rinsol(2,iseqinp)

      if ( iaccur>1000 ) then

!     --- ireltot>1000, hence eps is stored in a scalar

         A%eps = getscal(mod(iaccur,10))

      end if  ! ( iaccur>1000 )

      if ( mod(iaccur,1000)/=2 ) then
         A%eps1=0d0
         A%eps2=abs(A%eps)
      else
         A%eps1=A%eps
         A%eps2=0d0
      end if  ! ( ireler/=2 )

      if ( mod(iaccur,1000)>=100 ) A%eps1 = max(A%eps1,rinsol(7,iseqinp))

      if ( debug ) then
         write(irefwr,2) 'eps, eps1, eps2', A%eps, A%eps1, A%eps2
      end if  ! ( debug )
      A%alammin = 0d0
      if ( inpsol(15,iseqinp)==1 ) A%alammin = rinsol(4,iseqinp)

      A%type_prec = ipreco ! Type of preconditioner

      A%typesolv = typesolv            ! type of solver
      A%ldim = inpsol(14,iseqinp)      ! parameter l in bicgstab(l)
      A%nmax = nmax                    ! maximum number of iterations
      A%iprint = inpsol(6,iseqinp)     ! print level
      A%istartsol= inpsol(9,iseqinp)
      A%imax = max(1,inpsol(8,iseqinp)) ! Maximal dimension of Krylov space
      A%startvec = inpsol(9,iseqinp)   ! defines start vector
      A%iaccuracy = mod(iaccur,10)     ! defines how accuracy is stored
      A%iabserr = mod(iaccur,100)      ! check on absolute error
      A%irelerr = mod(iaccur,1000)     ! relative/abs error
      A%ntreigv = max(20,inpsol(13,iseqinp)) ! number of iterations used
                                       ! to estimate eigenvalue
      A%ntrunk = inpsol(13,iseqinp)    ! restart value for gmresr
                                       ! or maximum number of gcr iterations
      A%ninner = inpsol(14,iseqinp)    ! maximum number of inner iterations
                                       ! or sparam for idr
      A%aterror = inpsol(22,iseqinp)   ! defines what to do in case of error
      A%typeparsolv = typeparsolv      ! type of parallel solver

      A%eistat = A%type_prec==2 .or. A%type_prec==4 .or. A%type_prec==7
      A%precon = A%type_prec==3 .or. A%type_prec==5 .or. A%type_prec==8 &
                 .or. A%type_prec==9 .or. A%type_prec==10
      A%globshadow = A%ninner==1 .and. parallel .and. typesolv==10
      A%factor = rinsol(3,iseqinp)
      A%adapt = rinsol(2,iseqinp)
      A%iseqmt = iseqmt
      A%simple = IM%simple
      A%iseqim = iseqim

!     --- parameters for overrelaxation

      A%iconstrain = inpsol(17,iseqinp)
      A%iminimum = inpsol(18,iseqinp)
      A%imaximum = inpsol(19,iseqinp)
      A%idegfd = inpsol(20,iseqinp)
      A%minlev = rinsol(7,iseqinp)
      A%maxlev = rinsol(8,iseqinp)

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End sepcopyinpsol'

      end if  ! ( debug )
      call erclos ( 'sepcopyinpsol' )

      end subroutine sepcopyinpsol
