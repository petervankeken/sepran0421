      subroutine elfillcoeff ( number, user, coeffnodes, x, &
                               xgauss, phi, ichois, index, iuser, &
                               unodes, vecloc, coeffintpnts )
! ======================================================================
!
!        programmer    Guus Segal
!        version  2.1  date 12-12-2012 Long integers
!        version  2.0  date 17-07-2011 Remove xgauss, linear subinterpolation
!        version  1.0  date 19-02-2010
!
!   copyright (c) 2010-2012  "Ingenieursbureau SEPRA"
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
!     Fills function values in array coeffnodes for the nodal points
!     and in coeffintpnts for the integration points
!     Number corresponds to the number th coefficient in common block celiar
! **********************************************************************
!
!                       KEYWORDS
!
!     coefficient
!     element
!     fill
!     vector
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulebuild
      use sepmoduleoldvecs
      use sepmodulecactl
      use sepmodulecomio
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer number, ichois, index(n), iuser(*)
      double precision user(*), coeffnodes(*), x(n,*), xgauss(m,*), &
                       phi(n,*), unodes(*), &
                       vecloc(numold,maxunk,*), coeffintpnts(*)

!     coeffintpnts  i/o   Real array in which the values of the
!                         coefficients in the integration points are stored
!     coeffnodes    i/o   Real array in which the values of the
!                         coefficients in the nodal points are stored
!     ichois         i    Choice parameter to be used when
!                         1000<= ind(number)<2000
!                         Possible values:
!                         1:   In each node there are exactly ncomp unknowns
!                         2:   The number of unkowns per point may vary
!                              Array index contains the number of nodes per
!                              point accumulated. Hence
!                              number of unkowns in node 1 is index(1)
!                              number of unkowns in node i is
!                             index(i)-index(i-1)
!     index          i    Index array to be used if 1000<=ind(number)<2000 and
!                         ichois = 2
!     iuser          i    Integer user array to pass user information from
!                         main program to subroutine. See STANDARD PROBLEMS
!     number         i    The number-th unknown in each nodal point of the array
!                         is considered. If the array is complex number=2*i-1
!                         corresponds to the real part of array i and
!                         number=2*i to the imaginary part of array i
!     phi            i    array of length n * m, values of shape
!                         functions in integration points in the sequence:
!                         phi_i(xg^k) = phi (i,k)
!                         array phi is only filled when gauss integration is
!                         applied
!     unodes         i    Preceding solution in nodal points
!     user           i    Real user array to pass user information from
!                         main program to subroutine. See STANDARD PROBLEMS
!     vecloc         i    Work array in which all old solution vectors for the
!                         integration points are stored
!     x              i    array of length n x ndim containing the co-ordinates
!                         of the nodes
!     xgauss         i    array of length m x ndim containing the co-ordinates
!                         of the gauss points.  array xgauss is only filled when
!                         gauss integration is applied.
!                         xg_i = xgauss (i,1);  yg_i = xgauss (i,2);
!                         zg_i = xgauss (i,3)
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, j, istart, ifunc, ih, ndegfd, ivec, interpol
      double precision help, work(max(ncomp,(n+1)*m))
      logical fillgauss
      double precision, pointer :: vecoldlc(:)

!     fillgauss      If true Gauss points must be filled
!     help           help variable to compute a product
!     i              Counting variable
!     ifunc          Output parameter of elflr8 to be used if ind(number)==2004
!                    It indicates if the parameter which with the function
!                    is multiplied is a function (>0) or a constant (<0)
!     ih             Help variable to store a constant temporarily
!     interpol       If 0 the standard interpolation is applied, &
!                    if 1 the interpolation is restricted to linear
!                    subelements
!     istart         Starting address in array VECOLD
!     ivec           Sequence number of vector of special structure
!     j              Counting variable
!     ndegfd         Number of degrees of freedom per point
!     work           work array
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
      integer inivec
      double precision funccf, funcc1, funcc3

!     ELFILLCOEFF1   Fills array coeff with a coefficient, where the coefficient
!                    is equal to a solution or derivative of a preceding
!                    computation.
!     ELFILLCOEFFINT Fill coefficients computed in nodes in integration points
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     FUNCC1         User function subroutine to compute the value of a
!                    coefficient in a point depending on a preceding solution
!     FUNCC3         User function subroutine to compute the value of a
!                    coefficient in a point depending on a series of preceding
!                    solutions
!     FUNCCF         User function subroutine to compute the value of a
!                    coefficient in a point
!     INIVEC         Return with vector number corresponding to ivecin
!     PRININLONG     print 1d integer array
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
!    2048   Number of integration points in old solution not equal to that in
!           element
!    2050   Vector of type 129 multiplied by function not yet implemented
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
      call eropen ( 'elfillcoeff' )
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from elfillcoeff'
         write(irefwr,1) 'ifirst, number', ifirst, number
         write(irefwr,1) 'ich(number), ind(number)', &
                          ich(number), ind(number)
         write(irefwr,1) 'ist(number)', ist(number)
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      fillgauss = .true.  ! Fill also Gauss points
      interpol = irule/100

      if ( ifirst==0 ) then

!     --- ifirst = 0, first call for this element group

         if ( ich(number)>2000 ) then

!        --- Vector of type 129, check m

            ipos  = ivecnums(number)
            ivec = inivec ( iuser(ipos) )
            iseqtype = oldsols(ivec)%iseqtype
            index4n => indoldsols(iseqtype)%index4

            if ( index4n(2)/=m ) then
               call errint ( ivec, 1 )
               call errint ( index4n(2), 2 )
               call errint ( itype, 3 )
               call errint ( m, 4 )
               call errint ( ielgrp, 5 )
               call errsub ( 2048, 5, 0, 0 )
               go to 1000
            end if  ! ( index4n(2)/=m )

         end if  ! ( ich(number)>2000 )

      end if  ! ( ifirst==0 )

!     --- First coeffnodes is filled with values in the nodal points
!         if possible

      if ( ind(number)>=10000 ) then

!     --- ind(number)>10000, coefficients are computed by funcc3

         if ( interpol==1 ) then

!        --- interpol==1, possibly linear subinterpolation

            do i = 1, n
               coeffnodes(i) = funcc3 ( ind(number)-10000, x(i,1), &
                                        x(i,2), x(i,3), numold, &
                                        maxunk, vecloc(1,1,i) )
            end do  ! i = 1, n

         else

!        --- interpol==0, no linear subinterpolation

    !        do i = 1, n
    !           coeffnodes(i) = funcc3 ( ind(number)-10000, x(i,1), &
    !                                    x(i,2), x(i,3), numold, &
    !                                    maxunk, vecloc(1,1,i) )
     !       end do  ! i = 1, n

            do i = 1, n
               coeffintpnts(i) = funcc3 ( ind(number)-10000, &
                                          xgauss(i,1), xgauss(i,2), &
                                          xgauss(i,3), numold, maxunk, &
                                          vecloc(1,1,i) )
            end do  ! i = 1, n
            fillgauss = .false.  ! Gauss points have been filled

         end if  ! ( interpol==1 )

      else if ( ind(number)>=2003 ) then

!     --- ind(number)>2002, special treatment

         call elfillcoeff1 ( number, coeffnodes, iuser, ifunc )
         if ( ierror/=0 ) go to 1000

         if ( ich(number)>2000 ) then

!        --- Special situation of vector of type 129 with ich(number)-2000
!            degrees of freedom per element

           iseqtype = oldsols(ivec)%iseqtype
           index3n => indoldsols(iseqtype)%index3
            vecoldlc => oldsols(ivec)%oldvecs
            do i = 1, (ich(number)-2000)*m
               coeffintpnts(i) = vecoldlc ( index3n(i) )
            end do
            fillgauss = .false.  ! Gauss points have been filled

         end if  ! ( ich(number)>2000 )

         if ( ind(number)==2004 ) then

!        --- ind(number)=2004, coefficient must be multiplied

            if ( ich(number)<1000 ) then

!           --- Standard situation, one coefficient at a time

               if ( ifunc<0 ) then

!              --- ifunc<0, multiplication factor is constant

                  help = user(-ifunc)
                  do i = 1, n
                     coeffnodes(i) = coeffnodes(i) * help
                  end do

               else

!              --- ifunc>0, multiplication factor is a function

                  do i = 1, n
                     coeffnodes(i) = coeffnodes(i) * &
                        funccf ( ifunc, x(i,1), x(i,2), x(i,3) )
                  end do

               end if

            else if ( ich(number)<2000 ) then

!           --- Vectors of type 126, more coefficients at a time

               ndegfd = ich(number)-1000
               do j = 1, ndegfd

                  ih = (j-1)*n

                  if ( ifunc<0 ) then

!                 --- ifunc<0, multiplication factor is constant

                     help = user(-ifunc)
                     do i = 1, n
                        coeffnodes(i+ih) = coeffnodes(i+ih) * help
                     end do

                  else

!                 --- ifunc>0, multiplication factor is a function

                     do i = 1, n
                        coeffnodes(i+ih) = coeffnodes(i+ih) * &
                           funccf ( ifunc, x(i,1), x(i,2), x(i,3) )
                     end do

                  end if

               end do

            else

!           --- The option multiply vector of type 129 by constant or function
!               has not yet been implemented

               call errsub ( 2050, 0, 0, 0 )
               go to 1000

            end if

         end if

      else if ( ist(number)>0 ) then

!     --- ist(i)>0, copy coefficient from preceding one

         do i = 1, n
            coeffnodes(i) = coeffnodes(i+(ist(number)-number)*n)
         end do

      else if ( ind(number)<=1000 ) then

!     --- ind(number)<1000,  coefficient is a function of space alone

         if ( debug ) then
            write(irefwr,1) 'jdiag, m, kstep', jdiag, m, kstep
            write(irefwr,1) 'ind(number)', ind(number)
            call prinrl1 ( x, 2, n, 'x' )
         end if  ! ( debug )

         do i = 1, n
            coeffnodes(i) = funccf (ind(number), x(i,1), x(i,2), x(i,3))
         end do

         if ( debug ) call prinrl1 ( xgauss, 2, m, 'xgauss' )
         do i = 1, m
            coeffintpnts(i) = funccf (ind(number), xgauss(i,1), &
                                      xgauss(i,2), xgauss(i,3))
         end do
         fillgauss = .false.

         if ( debug ) call prinrl ( coeffnodes, n, 'coeffnodes' )

      else if ( ind(number)<=2000 ) then

!     --- ind(number)<2000,  coefficient is a function of space and a
!                            preceding solution

         if ( ichois==1 ) then

!        --- ichois = 1,  Constant number of degrees of freedom per point
!                         Store unknowns per point in work

            do i = 1, n
               do j = 1, ncomp
                  work(j) = unodes(i+(j-1)*n)
               end do
               coeffnodes(i) = funcc1 ( ind(number)-1000, x(i,1), &
                                        x(i,2), x(i,3), work )
            end do

         else

!        --- ichois = 2,  Variable number of degrees of freedom per point
!                         Store unknowns per point in work

            istart = 1
            do i = 1, n
               do j = istart, index(i)
                  work(j-istart+1) = unodes(j)
               end do
               istart = istart + index(i)+1
               coeffnodes(i) = funcc1 ( ind(number)-1000, x(i,1), &
                                        x(i,2), x(i,3), work )
            end do

         end if

      else if ( ind(number)==2001 ) then

!     --- ind(number) = 2001,  coefficient is a function of space alone
!                              it is given in all nodes

         ! index1 comes from sepmodulebuild
         !write(6,*) 'index1: ',index1(1:n)
         do i = 1, n
            coeffnodes(i) = user ( mst(number) + index1(i) )
         end do
         if ( inpelm/=n ) then

!        --- n # inpelm, i.e. 6-point triangle with unkowns in centroid

            coeffnodes(7) = ( -coeffnodes(1)-coeffnodes(3)-coeffnodes(5) &
               +4d0 * (coeffnodes(2)+coeffnodes(4)+coeffnodes(6)) )/9d0
         end if

      else

!     --- ind(number) = 2002,  coefficient is a function of space alone
!                              it is constant per element

         do i = 1, n
            coeffnodes(i) = user ( mst(number) + irelem )
         end do

      end if

!     --- Finally fill coefficients in integration points

      if ( .not. fillgauss ) go to 1000  ! jump, allready filled

      if ( jdiag==1 ) then

!     --- jdiag = 1, integration points and nodes coincide

         if ( ich(number)>1001 .and. ich(number)<2000 ) then

!        --- ndegfd degrees of freedom

            ndegfd = ich(number)-1000

         else

!        --- 1 degree of freedom

            ndegfd = 1

         end if  ! ( ich(number)>1001 .and. ich(number)<2000 )

         coeffintpnts(1:ndegfd*n) = coeffnodes(1:ndegfd*n)

      else

!     --- jdiag = 2, integration points and nodes are different

         call elfillcoeffint ( coeffintpnts, coeffnodes, phi, &
                               number )

      end if  ! ( jdiag==1 )

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End elfillcoeff'

      end if  ! ( debug )
      call erclos ( 'elfillcoeff' )

      end
