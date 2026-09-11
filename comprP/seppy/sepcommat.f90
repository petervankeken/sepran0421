      subroutine sepcommat ( iincommt )
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
!     compute structure of the large matrix
!     Actual body of commat and matstruc
!     Old subroutine incommat
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
      use sepmodulekmesh
      use sepmodulecpack
      use sepmodulematstr
      use sepmodulekprob
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer iincommt(*)

!     iincommt   i    Integer input array containing information concerning
!                     the structure of the large matrix
!                     See subroutine read00 for a description
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer isum, maxmat, isumnp, isumbc, lenim2
      type (intmatrix), pointer :: IM
      double precision :: timebefore(2)

!     IM             Refers to intmt(intmat)
!     isum           length of array imat part 2a
!     isumbc         length of array imat part 2b
!     isumnp         length of array imat part 2a minus the periodical boundary
!                    conditions part
!     lenim2         length of the array imat part 2
!     maxmat         Help parameter to store the number of entries in the large
!                    matrix
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLMN         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERSETTIME      Get the CPU time at start of subroutine
!     INCOMMATEXTEND Update the integer description of S_ii in the parallel
!                    case if necessary
!     INCOMMATPETSC  Fill array o_nnz with number of nonzeros in off-diagonal
!                    blocks for petsc
!     PRININ         print 1d integer array
!     PRININTMAT     Print array INTMAT
!     SEPCOMMATBOUN  Fill parts of array intmat with respect to
!                    boundary conditions
!     SEPCOMMATCOMP  Fill structure of matrix in case of a compact matrix
!     SEPCOMMATFILL  Fill information concerning the storage of the large matrix
!                    in array intmatinf
!     SEPCOMMATINTER Fill parts of array intmat with respect to interface
!                    unknowns
!                    Is only used in case of parallel computations
!     SEPCOMMATINV   Create structure of LU-decomposition matrix in case of
!                    a direct solver
!     SEPCOMMATREAC  Fill parts of array intmat with respect to boundary
!                    conditions.  Part reaction forces
!     SEPCOMMATRNES  If necessary renumber nodes for nested dissection
!                    Adapt kprobh
!     SEPCRINTMT     create a new array in intmt
!     SEPFILLCOMMAT  Fill modules for subroutine incommat
!     SEPFILLMATSTR1 Fill parameters from iincommt in module sepmodulematstr
!     SEPPROBLOOP2   information depending on node numberin with respect to
!                    the problem is computed and filled in array kprobg
!     SEPREPLMETHOD  Translate jmethod to logicals in type intmatrix
!     SEPSIMPLEPAR   Create two arrays sendinfo_u and sendinfo_p from sendinfo
!                    in case of parallel computing
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
!     --- Initialize some parameters

      call eropen ( 'sepcommat' )
      if ( itime>0 ) call ersettime ( timebefore )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepcommat'
         call prinin ( iincommt, 14, 'iincommt' )
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000
      nocallprobdf = .true.

!     --- Fill parameters from iincommt in module sepmodulematstr

      call sepfillmatstr1 ( iincommt )
      if ( ierror/=0 ) go to 1000

      debug=.true.
      if ( debug ) then
         write(irefwr,1) 'iprob, jmethod, iseqkm', iprob, jmethod, iseqkm
         write(irefwr,*) 'associated(km(iseqkm)%kmeshj), nested ', &
                          associated(km(iseqkm)%kmeshj), nested
      end if  ! ( debug )
      debug=.false.

!      --- Renumber nodes if necessary

      kmeshj => km(iseqkm)%kmeshj
      if ( .not. associated(km(iseqkm)%kmeshj) ) kelmj = 0
      if ( nested ) call sepcommatrnes

!     --- Fill kprob part depending on node numbering

      call sepprobloop2
      if ( ierror/=0 ) go to 1000

!     --- Fill modules

      call sepfillcommat ( iincommt )
      if ( ierror/=0 ) go to 1000

!     --- Create intmat structure

      call sepcrintmt ( intmat )
      call sepreplmethod          ! translate jmethod to logicals
      if ( ierror/=0 ) go to 1000

      IM => intmt(intmat)

      fillbouncond = iincommt(12)==0 .and. nbound>0 .or. npboun>0

      isumnp = 0

      isum = 0
      isumbc = 0
      isumnp = 0

      if ( debug ) then
         write(irefwr,1) 'nbound, npboun, nusol, nrusol, npoint', &
                          nbound, npboun, nusol, nrusol, npoint
         write(irefwr,1) 'iprob', iprob

      end if

      if ( fillbouncond ) then

!     --- Fill information with respect to boundary conditions

         call sepcommatboun ( isum, isumnp )

!        --- Fill information with respect to reaction forces

         call sepcommatreac ( isumbc )
         if ( ierror/=0 ) go to 1000

      end if  ! ( fillbouncond )

      if ( parallel .and. typecoupl>=1 ) then

!     --- Parallel case: new data structure
!         Fill part referring to interface unknowns

         call sepcommatinter
         if ( ierror/=0 ) go to 1000

      end if  ! ( typecoupl==1 )

      if ( debug ) &
         write(irefwr,1) 'isum, isumbc, isumnp', isum, isumbc, isumnp

!     --- Fill data structure for internal unknowns (row-compact storage)

      call sepcommatcomp ( lenim2 )
      if ( debug ) write(irefwr,*) 'after sepcommatcomp'
      if ( ierror/=0 ) go to 1000

!     --- Parallel case: dirichlet-dirichlet coupling
!         Check if irowsii and icolsii must be updated

      if ( typecoupl>=1 ) call incommatextend

!     --- In case of parallel petsc, fill array onnz
!         and store in kprob (par)

      if ( IM%petsc .and. parallel ) call incommatpetsc
      if ( ierror/=0 ) go to 1000

!     --- Fill information of storage in IM and remove irowsii/icolsii

      call sepcommatfill
      if ( debug ) write(irefwr,*) 'after sepcommatfill'
      if ( ierror/=0 ) go to 1000

!     --- Parallel case, neumann-neumann coupling (simple)

      if ( parallel .and. IM%simple ) call sepsimplepar
      kprobpart => kp(iseqkp)%kprobsub(iprob)
      kprobloc => kprobpart%kprobloc

!     --- Fill storage for LU-decomposition in irowsiiprec

      if ( IM%direct .and. nrusol>0 ) call sepcommatinv
      if ( debug ) write(irefwr,*) 'after sepcommatinv'
      if ( ierror/=0 ) go to 1000

      if ( ioutp>=1 .or. mod(iprint,10)>=1 ) then

!     --- ioutp>1   print size of large matrix

         if ( IM%complexm ) then

!        --- complex matrix

            write(irefwr, 610) size(IM%colsii)
610         format(' the large matrix contains', i10,' complex entries')

         else

!        --- real matrix

            if ( IM%nested ) then
               write(irefwr, 620) size(IM%colsii), IM%rowsiiprec(nrusol+1)
            else if ( IM%direct ) then
               write(irefwr, 620) size(IM%colsii), IM%rowsiiprec(nrusol)
            else
               write(irefwr, 630) size(IM%colsii)
            end if  ! ( IM%direct )
620         format(' the large matrix in compact form contains', i10, &
                   ' real entries'/ &
                   ' the large matrix in decomposed form contains', i10, &
                   ' real entries')
630         format(' the large matrix contains', i10, ' real entries')

         end if

      end if

      if ( debug .or. mod(iprint,10)>=2 ) then

!     --- print first part of array intmat

         call prinintmat ( 0 )

      end if ! ( debug .or. mod(iprint,10)>=2 )

      if ( allocated(kprobhglobinv) ) deallocate(kprobhglobinv)

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End sepcommat'

      end if  ! ( debug )

      call erclmn ( 'sepcommat', timebefore )
      end subroutine sepcommat
