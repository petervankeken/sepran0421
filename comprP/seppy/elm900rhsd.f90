      subroutine elm900rhsd ( elemvc, rho, phi, f, &
                              stress, rhospecial, ugauss, dudx, &
                              goertler, kappa, detdsc, w, dphidx, cn, &
                              mconv, jtime, mcoefconv, mcont, modelv, &
                              modelrhsd, secinv, xgauss, psiupw, x )
! ======================================================================
!
!        programmer    Guus Segal
!        version  3.0  date 06-01-2016 Use sepmoduleelem
!        version  2.0  date 27-02-2010 Change of parameter list
!        version  1.0  date 17-04-2007
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
!     Compute element vector for Navier-Stokes element
!
! **********************************************************************
!
!                       KEYWORDS
!
!     element_vector
!     navier_stokes_equation
! **********************************************************************
!
!                       MODULES USED
!
      use coeff
      use sepmoduleelem
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision elemvc(ndim,n), rho(m), phi(*), f(m,*), &
                       rhospecial(*), stress(*), &
                       ugauss(*), dudx(*), goertler(*), kappa(*), &
                       detdsc(*), w(m), cn, dphidx(*), secinv(*), &
                       xgauss(*), psiupw(*), x(n,ndim)
      integer mconv, jtime, mcoefconv, mcont, modelv, modelrhsd

!     cn             i    Power in power law model
!     detdsc         i    Derivative of second invariant per integration point
!     dphidx         i    Array of length n x m x ndim  containing the
!                         derivatives
!                         of the basis functions in the sequence:
!                         d phi / dx (xg ) = dphidx (i,k,1);
!                              i        k
!                         d phi / dy (xg ) = dphidx (i,k,2);
!                              i        k
!                         d phi / dz (xg ) = dphidx (i,k,3);
!                              i        k
!                         If the element is a linear simplex, the
!                         derivatives are constant and only k = 1 is filled.
!     dudx           i    Real m x ndim  array in which the values of the
!                         derivatives of the solution vector in the integration
!                         points are stored in the sequence:
!                         du / dx  (x ) = dudx ( k, j )
!                           j     k
!                         If jdercn = 0 in common block celint, then only
!                         dudx (1,*) is filled because dudx is then a constant
!                         per element
!                         dudx is only filled if jderiv>0
!     elemvc         o    Element vector
!     f              i    array of length m*ndim containing the right-hand side
!                         vector components
!     goertler       i    Goertler number
!     jtime          i    Integer parameter indicating if the mass matrix for
!                         the time-dependent equation must be computed (>0)
!                         or not ( = 0)
!     kappa          i    Parameter for the Goertler equations
!     mcoefconv      i    Defines how the convective term is treated:
!                         0: standard case
!                         1: special case: the convective term is defined as
!                               c_conv_1 u du/dx + c_conv_2 v du/dy
!                               c_conv_3 u dv/dx + c_conv_4 v dv/dy
!     mcont          i    Defines the type of continuity equation
!                         Possible values:
!                         0: incompressible flow
!                         1: compressible flow: div (rho u ) = 0
!                            Not applicable with penalty function method
!     mconv          i    see USER INPUT in subroutine ELM600
!     modelrhsd      i    Defines how the stress tensor is treated in the rhs
!                         Possible values:
!                         0:    if there is a given stress tensor the reference
!                               to this tensor is stored in position 19
!                         1:    It is supposed that there are 6 given stress
!                               tensor components stored in positions 19-24
!     modelv         i    Defines the type of viscosity model
!                         Possible values:
!                           1:  Newtonian liquid
!                           2:  power law liquid
!                           3:  Carreau liquid
!                           4:  plastico viscous liquid
!                         100:  user provided model depending on grad u
!                         101:  user provided model depending on the second
!                               invariant
!                         102:  user provided model depending on the second
!                               invariant, the co-ordinates and the velocity
!     phi            i    array of length n * m, values of shape
!                         functions in integration points in the sequence:
!                         phi_i(xg^k) = phi (i,k)
!                         array phi is only filled when gauss integration is
!                         applied
!     psiupw         i    Array of size n*m containing the adapted basis
!                         functions for upwind
!     rho            i    density rho
!     rhoc                Work array to store the adapted rho
!     rhospecial     i    Array with coefficients c_conv
!                         c_conv_1 = rho(.,1,1)
!                         c_conv_2 = rho(.,2,1)
!                         c_conv_3 = rho(.,1,2)
!                         c_conv_4 = rho(.,2,2)
!     secinv         i    Array containing the second invariant in the
!                         integration points
!     stress         i    Array of size 6 x nnodes containing the stress tensor
!     ugauss         i    Preceding solution in Gauss points
!     w              i    array of length m containing the weights for
!                         integration
!     x              i    array of length n x ndim containing the co-ordinates
!                         of the nodal points
!     xgauss         i    array of length m x ndim containing the co-ordinates
!                         of the gauss points.  array xgauss is only filled when
!                         gauss integration is applied.
!                         xg_i = xgauss (i,1);  yg_i = xgauss (i,2);
!                         zg_i = xgauss (i,3)
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer jdiagsav
      double precision adams, rhoc(m), work(n), work1(m*3)

!     adams          Multiplication factor for the Adams-Bashfort method
!     jdiagsav       Help parameter to save the value of jdiag
!     work           work array
!     work1          work array
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     EL2011         compute convective terms of Navier-stokes
!                    for newton iteration   (k-d)  (right-hand side)
!     EL2011G        compute convective terms of Goertler equations
!                    for newton iteration   (k-d)  (right-hand side)
!     EL2011SP       compute convective terms of Goertler equations
!                    for newton iteration   (k-d)  (right-hand side)
!                    (special form)
!     EL3003         Compute integrals of the shape / f(x) phi (x) dx
!                                                   e         i
!     EL4029         Fill element vector for newton linearization of power law
!                    viscosity terms
!     ELSTRRHS       Compute contribution of initial stress to right-hand side
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERRCHR         Put character in error message
!     ERRSUB         Error messages
!     PRINRL         Print 1d real vector
!     PRINRL1        Print 2d real vector
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!    2873   upwind not yet implemented
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
      call eropen ( 'elm900rhsd' )
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from elm900rhsd'
         write(irefwr,1) 'metupw, jdiag, ndim, n, m', &
                          metupw, jdiag, ndim, n, m
         write(irefwr,1) 'mconv, mcoefconv, jtime, modelv, mcont', &
                          mconv, mcoefconv, jtime, modelv, mcont
         write(irefwr,2) 'cn', cn
         call prinrl1 ( f, ndim, m, 'f' )
  1      format ( a, 1x, (10i6) )
  2      format ( a, 1x, (5d12.4) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      if ( metupw>0 ) then

!     --- upwind = true, basis functions, not diagonal
!         set jdiag = 2

         jdiagsav = jdiag
         jdiag = 2

      end if  ! ( upwind )

!     --- First set element vector equal to zero

      elemvc = 0d0

!     --- Next compute contribution due to external forces
!         Multiply weights by density rho and store in rhoc

      if ( debug ) call prinrl ( rho, m, 'rho' )
      if ( ind(3)/=0 .or. ind(4)/=0 .or. ind(5)/=0 ) then

!     --- At least one body force vector available

         rhoc = w * rho

      end if

!     --- First component

      if ( ind(3)/=0 ) then

!     --- f1 # 0  store contribution temporary in array work

         jadd = 0
         if ( metupw==0 ) then

!        --- No upwind required

             call el3003 ( f(1,1), work, rhoc, phi, m, n )
             if (ibuoy_trac>=1) then
                jadd=1
                call pel3003(f(1,1),work,rhoc,phi,m,n,work1,x,1,xgauss)
             endif

         else

!        --- Upwind

            call el3003 ( f(1,1), work, rhoc, psiupw, m, n )

         end if  ! ( metupw==0 )

         elemvc(1,1:n) = work
         if ( debug ) call prinrl ( work, n, 'work/1' )

      end if

!     --- Second component

      if ( ind(4)/=0 ) then

!     --- f2 # 0  store contribution temporary in array work

         jadd = 0
         if ( metupw==0 ) then

!        --- No upwind required

            call el3003 ( f(1,2), work, rhoc, phi, m, n )
             if (ibuoy_trac>=1) then
                jadd=1
                call pel3003(f(1,2),work,rhoc,phi,m,n,work1,x,2,xgauss)
             endif

         else

!        --- Upwind

            call el3003 ( f(1,2), work, rhoc, psiupw, m, n )

         end if  ! ( metupw==0 )

         elemvc(2,1:n) = work

      end if

!     --- Third component

      if ( ind(5)/=0 ) then

!     --- f3 # 0  store contribution temporary in array work

         jadd = 0
         if ( metupw==0 ) then

!        --- No upwind required

            call el3003 ( f(1,3), work, rhoc, phi, m, n )

         else

!        --- Upwind

            call el3003 ( f(1,3), work, rhoc, psiupw, m, n )

         end if  ! ( metupw==0 )

         elemvc(3,1:n) = work

      end if
      if ( debug ) call prinrl1 ( elemvc, n, ndim, 'elemvc/1' )

      if ( ind(10)/=0 ) then

!     --- Stress term corresponding to right-hand side

         if ( metupw==0 ) then

!        --- No upwind required

            call elstrrhs ( stress, dphidx, elemvc, &
                            w, m, n, ndimlc, phi, xgauss, &
                            modelrhsd  )

         else

!        --- Upwind not yet implemented

            call errchr ( 'Stress right-hand side', 1 )
            call errsub ( 2873, 0, 0, 1 )

         end if  ! ( metupw==0 )

      end if
      if ( mconv==2 .or. mconv ==4 ) then

!     --- compute contribution of convective terms

         if ( mconv==4 ) then
            if ( jtime>0 ) then
               adams = -23d0/12d0
            else
               adams = -1d0
            end if
         else
            adams = 1d0
         end if

         if ( mcoefconv>0 ) then

!        --- compute convective terms where the coefficients are stored
!            in the special way

            call el2011sp ( adams, elemvc, phi, rhospecial, &
                            w, ugauss, dudx, work )

         else if ( mcont==2 ) then

!        --- compute convective terms in case of Goertler equations

            call el2011g ( adams, elemvc, phi, rho, w, &
                           ugauss, dudx, xgauss, work, &
                           goertler, kappa )

         else if ( metupw==0 ) then

!        --- compute convective terms in standard case

            call el2011 ( adams, elemvc, phi, rho, w, &
                          ugauss, dudx, xgauss, work )

         else

!        --- compute convective terms in upwind case

            call el2011 ( adams, elemvc, psiupw, rho, w, &
                          ugauss, dudx, xgauss, work )

         end if

      end if

      if ( modelv==2 .and. cn>1d0 ) then

!     --- modelv = 2, and cn>1:  Newton linearization

         call el4029 ( elemvc, xgauss, ugauss, &
                       dudx, secinv, detdsc, &
                       work, dphidx, w )
      end if

      if ( metupw>0 ) then

!     --- reset jdiag

         jdiag = jdiagsav

      end if  ! ( upwind )

1000  if ( debug ) then

!     --- Debug information

         call prinrl1 ( elemvc, n, ndim, 'elemvc/end' )
         write(irefwr,*) 'End elm900rhsd'

      end if  ! ( debug )
      call erclos ( 'elm900rhsd' )

      end
