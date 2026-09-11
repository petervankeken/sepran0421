      subroutine dtm_p_mass( pp, w, npres, viscos, psi)
! ======================================================================
!
!        programmer    Guus Segal
!                           10-11-2017 copied from elm900mat to
!                           produce only scaled pressure mass matrix
!
!   copyright (c) 2007-2016  "Ingenieursbureau SEPRA"
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
!     Store extra element matrices which are necessary for
!     schur complement type methods
! **********************************************************************
!
!                       KEYWORDS
!
!     element_matrix
!     navier_stokes_equation
! **********************************************************************
!
!                       MODULES USED
!
      use sepmoduleelem
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer npres
      double precision w(2*m), viscos(m), psi(npres,m), pp(npres,npres)

!     elemmt        i/o   Element matrix to be filled
!                         In this subroutine we fill the part corresponding
!                         to scaled pressure mass matrix, diagonal of
!                         pressure mass matrix and velocity mass matrix
!     iseqvel        i    Array containing the positions of the velocities in
!                         the element vector in the sequence required by el2005
!     npres          i    Number of pressure parameters per element
!     psi            i    Array of length npsi x m containing the values of the
!                         basis functions for the pressure in the integration
!                         points
!     viscos         i    Array of size m containing the viscosity in the
!                         integration points
!     w              i    array of length m containing the weights for
!                         integration
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, ipweipmm, ipdiagpmm, ipdiagvmm, k, nvel, nsave
      integer :: j,jdiagsav
      double precision kappa(m), sum, f(npres), g(m), work(m)

!     i              Counting variable
!     ipdiagpmm      Starting address of diagonal of pressure mass matrix
!                    in elemmt
!     ipdiagvmm      Starting address of diagonal of velocity mass matrix
!                    in elemmt
!     ipweipmm       Starting address of scaled pressure mass matrix
!                    in elemmt
!     k              Counting variable
!     kappa          Cis used to store the inverse of the viscosity
!     nsave          Help variable to store old value of n
!     nvel           Number of velocity points
!     sum            Help parameter to compute a sum
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ELFILLPRESMMT  Fill pressure mass matrix weighted by coefficient kappa
!     ELMASSDIAG     Copy components-wise defined diagonal element mass matrix
!                    into diagonal element mass matrix
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
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
      call eropen ( 'dtm_p_mass' )
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from dtm_p_mass'

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

!     --- Store weighted pressure mass matrix, diagonal of pressure mass
!         matrix and diagonal of velocity mass matrix into elemmt
!         from position icount^2+1
!         icount^2+1 ... icount^2+1+npres^2 weighted pressure mass matrix
!         icount^2+1+npres^2+1 ... icount^2+1+npres^2+npres diag press mm
!         icount^2+1+npres^2+npres+1 ... icount^2+1+npres^2+npres+icount vel mm


!     --- Fill weighted pressure mass matrix

      kappa = 1d0/viscos
      nsave = n
      n = npres
      jdiagsav = jdiag
      jdiag = 2

!     --- Fill the pressure matrix according to
!         PP(i,j) = / kappa * psi(i) * psi(j) dx
!                   x

      !print_matrix write(irefwr,*) 'pressure mass matrix'
      do j = 1, npres

!     --- Re-initialise array WORK2

         f = 0d0

         g = kappa*psi(j,1:m)

         call el3003( g, f, w, psi, m, npres, work )

         pp(1:npres,j) = f
         !print_matrix  write(irefwr,'(3e15.7)') pp(1:npres,j)

      end do

      n = nsave
      jdiag = jdiagsav



1000  call erclos ( 'dtm_p_mass' )
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End dtm_p_mass'

      end if  ! ( debug )

      end
