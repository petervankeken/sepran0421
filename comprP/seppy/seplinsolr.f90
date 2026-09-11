      subroutine seplinsolr ( matr, isol, irhsd )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.1  date 11-12-2017 New storage of intmat
!        version  1.0  date 23-02-2016
!
!   copyright (c) 2016-2017  "Ingenieursbureau SEPRA"
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
!     Solve system of linear equations by either a direct or an iterative
!     method (Real matrix)
! **********************************************************************
!
!                       KEYWORDS
!
!     linear
!     linear_solver
!     solve
! **********************************************************************
!
!                       MODULES USED
!
      use sepmoduleintmat
      use sepmodulesolve
      use sepmodulematr
      use sepmodulevecs
      use sepmodulekmesh
      use sepmoduleproj
      use sepmodulekprob
      use sepmodulecactl
      implicit none

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer :: matr(*), isol(*), irhsd(*)

!     irhsd        i    Standard SEPRAN array, containing information of the
!                         large vector.
!                         If defect correction is applied, irhsd must contain
!                         two vectors and hence must be declared as irhsd(5,2)
!                         irhsd(.,1) corresponds to matr(.,1) including the
!                                   effect of boundary conditions
!                         irhsd(.,2) corresponds to matr(.,2) including the
!                                   effect of boundary conditions
!     isol           i    Standard SEPRAN array containing information of the
!                         solution array
!     matr           i    Standard SEPRAN array, containing information of the
!                         large matrix
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer :: iinvec(6), idummy(2)
      double precision :: rinvec(2)
      type (globmat), pointer :: A
      type (rmatrix), pointer :: Ar
      type (intmatrix), pointer :: IM
      double precision :: timebefore(2)

!     A              Structure of type globmat containing the necessary
!                    information of the matrix
!     Ar             Structure of type rmatrix containing the necessary
!                    information of the real matrix
!     IM             Refers to intmt(intmat)
!     idummy         Dummy parameter
!     iinvec         Array containing the integer input for subroutine MANVEC
!     rinvec         Array containing the real input for subroutine MANVEC
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLMN         Resets old name of previous subroutine of higher level
!     ERCLOS         Resets old name of previous subroutine of higher level
!     ERDEALLOC      Produce error message in case deallocate went wrong
!     EROPEN         Produces concatenated name of local subroutine
!     ERSETTIME      Get the CPU time at start of subroutine
!     PRESIDUAL      Compute residual of system of linear equations
!     SEPCOMPREAC1   Actual computing of reaction force
!     SEPCOPPROJ     Copy information of projection vectors in type matrix
!     SEPCOPYINPSOL  Put information of inpsol and rinsol in matrix structure
!     SEPLINSOLBCK   Renumber solution vectors back for subroutine linsol
!     SEPLINSOLREAL  Solve system of linear equations (real matrix)
!     SEPPETSCMAT    Transform sepran matrix into petsc matrix
!     SEPPJPOST      contour plots are made of the projection vectors
!     SEPPRINTV      Print solution vector or vector of special structure
!     SEPPROJMK      Creates projection vectors for the deflated iccg algorithm
!     SEPRENUMSOLBCK Renumber solution vector back to original form
!                    Real case
!     SEPRMMATSTR2   nullify arrays in module sepmodulematstr
!     SEPSOLDIR      Actual body of subroutine SOLVE
!     SEPSOLVEDEFCOR apply defect correction
!     SEPSOLVELITREF iterative refinement and/or
!                    Compute and print the residual
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
!   If ( typesolv=0 ) then
!      Solve by profile method
!   else if ( typesolv=5 ) then
!      if ( inpsol(17)=0 ) then
!         Solve by overrelaxation (OVERRL)
!      else
!         Solve by constraint overrelaxation (OVERCS)
!   else if ( typesolv>0 ) then
!      Solve by CG-type iterative method
!      if ( inpsol(17)=1 ) then
!         Create projection vectors
!   else ( typesolv=-1 ) then
!      Solve by Y12M
!   if ( direct solver is applied and (nmax>0 and profile method
!        is applied or irescomp>0 ) then
!      niter := 0
!      maxitr := max (1, nmax)
!      While ( niter<maxitr ) do
!         niter := niter + 1
!         Compute residual
!         if ( profile method and nmax>0 ) then
!            Perform one stpe iterative improvement
!         if ( irescomp>0 ) then
!            Print residual
!   if ( defect correction must be applied ) then
!      niter := 0
!      maxitr := max (1, numdef_corr)
!      While ( niter<maxitr ) do
!         niter := niter + 1
!         Compute residual (matr(2) usol - rhsd(2))
!         Solve: matr(1) ucorr = -residual
!         Add ucorr to usol
!   if ( iseq_react>0 ) then
!      Compute reaction force
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen ( 'seplinsolr' )
      if ( itime>0 ) call ersettime ( timebefore )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from seplinsolr'
         write(irefwr,1) 'ireadinpsol,iseqintmat', ireadinpsol,iseqintmat
  1      format ( a, 1x, (10i6) )

      end if
      if ( ierror/=0 ) go to 1000

      iseqmt = matr(1)
      A => realmt(iseqmt)
      Ar => A%S_i
      iseqim = A%intmatref
      copyprec = .true.   ! Default

      IM => intmt(iseqim)

      if ( A%diagonal ) then

!     --- Diagonal matrix

         call sepsoldir ( matr, isol, irhsd )
         go to 1000

      end if  ! ( A%diagonal )

!     --- Clear arrays

      projdata = 0

      if ( imetproj==1 ) then

!     --- Create projection vectors

         call sepprojmk
         if ( plot_lv>=1 .and. ndim==2 .and. .not. parallel ) &
            call seppjpost

      end if  ! ( imetproj==1 )

!     --- Make global rhsd in case of parallel computing

      Ar%decomp = .false.
      if ( nproj>0 ) call sepcopproj ( Ar )

!     --- Put information of inpsol and rinsol in matrix structure

      write(6,*) 'Ar b  : ',Ar%typesolv
      call sepcopyinpsol ( Ar, 1 )
      write(6,*) 'Ar a : ',Ar%typesolv

!     --- Store sepran matrix into petsc matrix
!         in case of jmethod = 56, this is done in build

      if ( IM%petsc .and. .not. IM%nosepmat ) call seppetscmat
      if ( ierror/=0 ) go to 1000

      call seplinsolreal ( Ar, usol, rhsd )

      if ( typesolv<=0 .and. (nmax>0 .or. irescomp>0) )  then

!     --- nmax>0: iterative refinement and/or
!         irescomp>0 Compute and print the residual

         call sepsolvelitref ( Ar, usol, rhsd )

      end if

      if ( numdef_corr>0 ) call sepsolvedefcor  !defect correction

!     --- Copy solution array back into usolbuf

      call seprenumsolbck ( usol, 1, xglob )
      if ( indprh>0 ) then
         deallocate ( usol, rhsd, stat = error )
         if ( error/=0 ) call erdealloc ( error, 'usol' )
      else
         nullify(usol)
         nullify(rhsd)
      end if  ! ( indprh>0 )

      if ( parallel ) then

!     --- Deallocate arrays

         deallocate ( xglob, sendwork, stat = error )
         if ( error/=0 ) call erdealloc ( error, 'sendwork' )

         if ( associated(sendunknowns) ) nullify (sendunknowns)
         if ( associated(recunknowns)) nullify (recunknowns)

      end if  ! ( parallel )

      if ( iseqreac/=0 ) then

!     --- Compute reaction force

         Smat => A
         call sepcompreac1 ( usolbuf, reac, rhsdbuf )

      end if

      if ( ilimit>0 .and. valmax>valmin ) &
         usolbuf = max(min(usolbuf,valmax),valmin) ! limit solution

      if ( iseq_end_res>0 ) then

!     --- iseq_end_res>0, compute residual of end solution

         call presidual ( matr, isol(iseqsol), &
                          irhsd, isol(iseq_end_res), 1 )

      end if

      numiters(iprob) = inpsol(1,1)

      if ( (keep==0 .or. lenpremat<=0) .and. &
           (keepmat==0 .or. keepmat==2) ) then

!     --- preconditioning matrix must be removed

         if ( associated(realmt(iseqmt)%matinv) ) &
            deallocate(realmt(iseqmt)%matinv)
         nullify (realmt(iseqmt)%S_i%MD )

      end if  ! ( keep==0 .or. lenpremat<=0 )

!     --- Renumber solutions back to standard sequence

      if ( indprh>0 ) call seplinsolbck ( .false. )
      call seprmmatstr2  ! remove associations with arrays in sepmodulematstr

1000  if ( debug ) then

!     --- Debug information

         call sepprintv ( isol(iseqsol), 'solution' )
         write(irefwr,*) 'End seplinsolr'

      end if
      if ( ierror==0 ) then
         call erclmn ( 'seplinsolr', timebefore )
      else
         call erclos ( 'seplinsolr' )
      end if  ! ( ierror==0 )

      end subroutine seplinsolr
