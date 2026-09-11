      subroutine sepmumpsbody ( numnonzeros, mumps_par, new , usol, rhsd )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 14-04-2016
!
!   copyright (c) 2016-2016  "Ingenieursbureau SEPRA"
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
!     Actual body for the call to the mumps solver
!     Old subroutine prmumpsbody
! **********************************************************************
!
!                       KEYWORDS
!
!     linear_solver
!     solve
!     mumps
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulematr
      use sepmodulekprob
      use sepmodulematstr
      use sepmodulesolve
      use sepmodulecpack
      use control
      implicit none
      include 'mpif.h'
      include 'dmumps_struc.h'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer, intent(in) :: numnonzeros
      type (dmumps_struc) :: mumps_par
      logical :: new
      double precision :: usol(nrusol)
      double precision, intent(in) :: rhsd(nrusol)

!     mumps_par      i    standard mumps structure
!     new            i    if true a new structure of the matrix has been made
!     numnonzeros    i    number of nonzeros in matrix
!     rhsd           i    right-hand-side vector
!     usol           o    solution vector
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer :: ierr, iref
      logical :: extradebug, firstcall
      save iref

!                    they are not
!     extradebug     If true more debug statements are carried out otherwise
!                    they are not
!     firstcall      Indicates if this is the first call of prmupsbody (true)
!                    or not (false)
!     ierr           Error indicator of mpi
!     iref           File reference number for mumps information
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
      integer inisetnewref

!     DMUMPS         MUMPS solver
!     ERALLOC        Produce error message in case allocate went wrong
!     ERCLOS         Resets old name of previous subroutine of higher level
!     ERDEALLOC      Produce error message in case deallocate went wrong
!     EROPEN         Produces concatenated name of local subroutine
!     INISETNEWREF   Create new reference number
!     MPI_INIT       Start MPI
!     PRININ         print 1d integer array
!     PRINRL         Print 1d real vector
!     PRMUMPSINFO    Write error message or warning in case of a mumps error
!     PRMUMPSROWIND  Fill array rowindex for mumps using irowsii
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
      data firstcall /.true./
! ======================================================================
!
      call eropen ( 'sepmumpsbody' )
      debug = .false. .and. ioutp>=0
      extradebug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepmumpsbody'
         write(irefwr,1) 'numnonzeros', numnonzeros
         write(irefwr,*) 'firstmumps, new', firstmumps, new
         write(irefwr,*) 'mpistart', mpistart
         write(irefwr,1) 'mumps_par%JOB', mumps_par%JOB
         write(irefwr,1) 'keep', keep
         if ( extradebug .and. .not. skipmumps ) then
            call prinrl ( rhsd, nrusol, 'rhsd' )
            call prinrl ( S_ii, numnonzeros, 'S_ii' )
            call prinin ( icolsii, numnonzeros, 'icolsii' )
         end if  ! ( .not. skipmumps )
  1      format ( a, 1x, (10i6) )
      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      usemumps = .true.

      if ( firstmumps .and. .not. mpistart ) then

!     --- First call of mumps, open mpi

         call mpi_init(ierr)
         firstmumps = .false.
         mpistart = .true.

      end if  ! ( firstmumps )

      if ( firstcall ) then

!     --- first call to sepmumpsbody, open file mumps.info

         iref = inisetnewref ()
         open ( unit = iref, file  = 'mumps.info' )
         firstcall = .false.

      end if

      if ( new .and. mumps_par%JOB/=0 ) then

!     --- New structure, but the mumps structure for this equation already exists
!         Destroy old structure

         mumps_par%JOB = -2
         call dmumps(mumps_par)

         if ( mumps_par%info(1)/=0 ) call prmumpsinfo ( mumps_par%info(1) )

      end if  ! ( new )

      if ( new .and. keep<2 ) then

!     --- Set some parameters
!         new matrix
!         Define a communicator for the package.

         mumps_par%COMM = MPI_COMM_WORLD
         mumps_par%PAR = 1    ! host is involved in factorization/solve phases

         if ( debug ) write(irefwr,1) 'mumps_par%SYM/1', mumps_par%SYM

!        --- Set type of matrix

         select case ( jmethod )
         case(52)
            mumps_par%SYM = 0  ! unsymmetric
         case(53)
            mumps_par%SYM = 2  ! symmetric
         case(54)
            mumps_par%SYM = 1  ! symmetric positive definite
         end select  ! case ( jmethod )
         if ( debug ) then
            write(irefwr,1) 'mumps_par%SYM/2', mumps_par%SYM
            write(irefwr,1) 'mumps_par%JOB', mumps_par%JOB
         end if  ! ( debug )

!        --- Initialize mumps

         mumps_par%JOB = -1   ! initialize mumps
         call dmumps(mumps_par)

         mumps_par%ICNTL(3) = iref  ! file number
         mumps_par%ICNTL(7) = mumps_ren  ! renumbering scheme
         mumps_par%ICNTL(14) = 400  ! percentage increase in the estimated
                                    ! work space
         if (no_mumps_info) then
            !write(6,*) 'mumps_par%ICNTL: ',mumps_par%ICNTL(1:7)
            mumps_par%ICNTL(2) = 0 ! set diagnostic output to none...
            mumps_par%ICNTL(3) = 0 
         endif

      end if  ! ( keep<2 ) then

!     --- Define problem on the host (processor 0)

      if ( mumps_par%MYID==0 ) then

!     --- Host only

         mumps_par%N = nrusol
         mumps_par%NZ = numnonzeros

         allocate ( mumps_par%IRN(numnonzeros), mumps_par%RHS(nrusol), &
                    stat = error )
         if ( error/=0 ) &
            call eralloc ( error, numnonzeros+nrusol,'rowindex' )
         if ( ierror/=0 ) go to 1000

         mumps_par%JCN => icolsii   ! jcn is identified as icolsii
         mumps_par%A => S_ii            ! A is identified as S_ii
         mumps_par%RHS(1:nrusol) = rhsd(1:nrusol)  ! RHS is identified as rhsd

!        --- Create array mumps_par%IRN with the row numbers of the matrix
!            entries

         call prmumpsrowind ( nrusol, numnonzeros, irowsii, &
                              mumps_par%IRN )

      end if  ! ( mumps_par%MYID==0 )

!     --- Perform analysis

      if ( new .and. keep<2 ) then

!     --- new matrix

         mumps_par%JOB = 1
         call dmumps(mumps_par)
         if ( mumps_par%info(1)/=0 ) call prmumpsinfo ( mumps_par%info(1) )
         if ( ierror/=0 ) go to 1000

      end if  ! ( new .and. keep<2 )

      if ( keep<2 ) then

!     --- Compute the LU factorization

         mumps_par%JOB = 2
         call dmumps(mumps_par)

         if ( mumps_par%info(1)/=0 ) call prmumpsinfo ( mumps_par%info(1) )
         if ( ierror/=0 ) go to 1000

      end if  ! ( keep<2 )

!     --- Solve system of equations

      mumps_par%JOB = 3
      call dmumps(mumps_par)

      if ( mumps_par%info(1)/=0 ) call prmumpsinfo ( mumps_par%info(1) )
      if ( ierror/=0 ) go to 1000

      if ( mumps_par%MYID==0 ) then

!     --- Host only
!         Copy solution into usol

         usol(1:nrusol) = mumps_par%RHS(1:nrusol)

         deallocate ( mumps_par%IRN, mumps_par%RHS, stat = error )
         if ( error/=0 ) call erdealloc ( error, 'rowindex' )

      end if  ! ( mumps_par%MYID==0 )

1000  if ( debug ) then

!     --- Debug information

         if ( extradebug .and. .not. skipmumps ) &
            call prinrl ( usol, nrusol, 'usol' )
         write(irefwr,*) 'End sepmumpsbody'

      end if  ! ( debug )
      call erclos ( 'sepmumpsbody' )

      end
