      subroutine sepelfillvar1 ( number, vecloc )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 12-08-2015
!
!   copyright (c) 2015-2015  "Ingenieursbureau SEPRA"
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
!     Fill one space dependent coefficient
!
! **********************************************************************
!
!                       KEYWORDS
!
!     element
!     fill
!     coefficient
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulecactl
      use sepmoduleelm
      use sepmodulebasefn
      use sepmodulecomio
      use sepmoduleoldvecs
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer, intent(in) :: number
      double precision, intent(in) :: vecloc(1:numoldvecs,1:E%maxunk,1:m)

!     number         i    Parameter indicating which coefficient
!                         must be considered
!     vecloc         i    Work array in which all old solution vectors for the
!                         quadrature points are stored
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer :: ivec, ifunc, i, j, ndegfd, i2, i3
      double precision, pointer :: coeff(:,:), coeffnodes(:,:)
      integer, pointer :: index1(:)
      double precision :: work(n), work1(E%maxunk)

!     coeff          Array containing the coefficients in the quadrature points
!     coeffnodes     Real array in which the values of the
!                    coefficients in the nodal points are stored
!     i              Loop variable
!     i2             min (2,ndim)
!     i3             min (3,ndim)
!     ifunc          Indicates if the parameter which the function
!                    is multiplied by is a function (>0) or a constant (<0)
!     index1         Contains the node numbers of the element
!     ivec           Sequence number of vector in set of old vectors
!     j              Loop variable
!     ndegfd         Number of degrees of freedom per point
!     work           Work array
!     work1          Work array to store components of vector in 1 point
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
      integer :: inivec
      double precision :: funccf, funcc1, funcc3

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
!     PRINRL         Print 1d real vector
!     SEPELFILLVAR2  Fills array coeff with a coefficient, where the
!                    coefficient is equal to a solution or derivative of
!                    a preceding computation
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
      call eropen ( 'sepelfillvar1' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from sepelfillvar1'
         write(irefwr,1) 'ifirst, number', E%first, number
         write(irefwr,1) 'E%ich(number), E%ind(number)', &
                          E%ich(number), E%ind(number)
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      i2 = min(2,E%ndim)
      i3 = min(3,E%ndim)

      if ( E%first==0 ) then

!     --- ifirst = 0, first call for this element group

         if ( E%ich(number)>2000 ) then

!        --- Vector of type 129, check m

            ivec = inivec ( E%ivecnums(number) )
            iseqtype = oldsols(ivec)%iseqtype

            if ( indoldsols(iseqtype)%index4(2)/=m ) then
               call errint ( ivec, 1 )
               call errint ( indoldsols(iseqtype)%index4(2), 2 )
               call errint ( E%itype, 3 )
               call errint ( m, 4 )
               call errint ( E%ielgrp, 5 )
               call errsub ( 2048, 5, 0, 0 )
               go to 1000
            end if  ! ( index4(ivec,2)/=m )

         end if  ! ( E%ich(number)>2000 )

      end if  ! ( ifirst==0 )
      coeff => E%coeffs
      coeffnodes => E%coeffsnode
      index1 => E%nodes(1:E%npelm)

      if ( E%ind(number)>=10000 ) then

!     --- E%ind(number)>10000, coefficients are computed by funcc3

         if ( B%jdiag==1 ) then

!         --- jdiag = 1,  Nodal points

            do i = kstep, m, kstep
               coeff(i,number) = funcc3 ( E%ind(number)-10000, B%xyz(i,1), &
                                          B%xyz(i,i2), B%xyz(i,i3), &
                                          numoldvecs, E%maxunk, vecloc(1,1,i) )
            end do

         else

!        --- jdiag = 2,  Integration points

            do i = 1, m
               coeff(i,number) = funcc3 ( E%ind(number)-10000, &
                                          B%xyzg(i,1), B%xyzg(i,i2), &
                                          B%xyzg(i,i3), numoldvecs, &
                                          E%maxunk, vecloc(1,1,i) )
            end do

         end if  ! ( B%jdiag==1 )


      else if ( E%ind(number)==2003 .or. E%ind(number)==2004 ) then

!     --- E%ind(number)>2003, special treatment

         call sepelfillvar2 ( number )
         if ( ierror/=0 ) go to 1000
         if ( E%ind(number)==2004 ) then

!        --- E%ind(number)=2004, coefficient must be multiplied

            ipos = E%ivecnums(number)
            ifunc = E%iuser(ipos+2)

            if ( debug ) write(irefwr,1) 'ipos, ifunc', E%ipos, ifunc

            if ( E%ich(number)<1000 ) then

!           --- Standard situation, one coefficient at a time

               if ( ifunc<0 ) then

!              --- ifunc<0, multiplication factor is constant

                  do i = B%kstep, m, B%kstep
                     coeff(i,number) = coeff(i,number) * E%user(-ifunc)
                  end do

               else

!              --- ifunc>0, multiplication factor is a function

                  if ( B%jdiag==1 ) then

!                 --- jdiag = 1,  Nodal points

                     do i = B%kstep, m, B%kstep
                        coeff(i,number) = coeff(i,number) * &
                           funccf ( ifunc, B%xyz(i,1), B%xyz(i,i2),&
                                    B%xyz(i,i3) )
                     end do
                  else

!                 --- jdiag = 2,  Integration points

                     do i = 1, m
                        coeff(i,number) = coeff(i,number) * &
                             funccf ( ifunc, B%xyzg(i,1), B%xyzg(i,i2), &
                             B%xyzg(i,i3) )
                     end do
                  end if
               end if
            else if ( E%ich(number)<2000 ) then

!           --- Vectors of type 126, more coefficients at a time

               ndegfd = E%ich(number)-1000
               do j = number-1, number+ndegfd-1
                  if ( ifunc<0 ) then

!                 --- ifunc<0, multiplication factor is constant

                     do i = B%kstep, m, B%kstep
                        coeff(i,j) = coeff(i,j) * E%user(-ifunc)
                     end do
                  else

!                 --- ifunc>0, multiplication factor is a function

                     if ( B%jdiag==1 ) then

!                    --- jdiag = 1,  Nodal points

                        do i = B%kstep, m, B%kstep
                           coeff(i,j) = coeff(i,j) * funccf ( ifunc, &
                                         B%xyz(i,1), B%xyz(i,i2), B%xyz(i,i3) )
                        end do
                     else

!                    --- jdiag = 2,  Integration points

                        do i = 1, m
                           coeff(i,j) = coeff(i,j) * funccf ( ifunc, &
                                         B%xyzg(i,1), B%xyzg(i,i2), &
                                         B%xyzg(i,i3) )
                        end do
                     end if
                  end if
               end do
            else

!           --- The option multiply vector of type 129 by constant or function
!               has not yet been implemented

               call errsub ( 2050, 0, 0, 0 )
               go to 1000
            end if

         end if

      else if ( E%ist(number)>0 ) then

!     --- E%ist(i)>0, copy coefficient from preceding one

         do i = B%kstep, m, B%kstep
            coeff(i,number) = coeff(i,E%ist(number))
         end do

      else if ( E%ind(number)<=1000 ) then

!     --- E%ind(number)<1000,  coefficient is a function of space alone

         if ( debug ) write(irefwr,1) 'jdiag, m, kstep', B%jdiag, m, B%kstep
         if ( B%jdiag==1 ) then

!        --- jdiag = 1,  Nodal points

            do i = B%kstep, m, B%kstep
               coeff(i,number) = &
                  funccf ( E%ind(number), B%xyz(i,1), B%xyz(i,i2), B%xyz(i,i3) )
            end do
            coeffnodes(:,number) = coeff(:,number)

         else

!        --- jdiag = 2,  Integration points

            if ( itype>=900 .and. itype<=902 ) then
               if ( debug ) call prinrl1 ( B%xyz, 2,n, 'x' )
               do i = 1, n
                  coeffnodes(i,number) = funccf ( E%ind(number), B%xyz(i,1), &
                                               B%xyz(i,i2), B%xyz(i,i3) )
               end do
            end if  ! ( itype>=900 .and. itype<=902 )

            if ( debug ) call prinrl1 ( B%xyzg(1,1), 2, m, 'xgauss' )
            do i = 1, m
               coeff(i,number) = funccf ( E%ind(number), B%xyzg(i,1), &
                                          B%xyzg(i,i2), B%xyzg(i,i3) )
            end do

         end if

      else if ( E%ind(number)<=2000 ) then

!     --- E%ind(number)<2000,  coefficient is a function of space and a
!                              preceding solution

         do i = 1, m

            work1(1:E%maxunk) = E%ugauss(i,1:E%maxunk)

            coeff(i,number) = funcc1 ( E%ind(number)-1000, B%xyzg(i,1), &
                                       B%xyzg(i,i2), B%xyzg(i,i3), work1 )

         end do  !  i = 1, m


      else if ( E%ind(number)==2001 ) then

!     --- E%ind(number) = 2001,  coefficient is a function of space alone
!                              it is given in all nodes

         if ( B%jdiag==1 ) then

!        --- jdiag = 1

            do i = B%kstep, m, B%kstep
               coeff(i,number) = E%user ( E%mst(number) + index1(i) )
            end do
            if ( E%npelm/=n ) then

!           --- n # inpelm, i.e. 6-point triangle with unkowns in centroid

               coeff(7,number) = &
                  ( -coeff(1,number)-coeff(3,number)-coeff(5,number)+ &
                 4d0*(coeff(2,number)+coeff(4,number)+coeff(6,number))) /9d0

            end if
            coeffnodes(1:7,number)=coeff(1:7,number)

         else

!        --- jdiag = 2

            do i = 1, n
               work(i) = E%user ( E%mst(number) + index1(i) )
            end do
            if ( E%npelm/=n ) then

!           --- n # inpelm, i.e. 6-point triangle with unkowns in centroid

               work(7) = ( -work(1)-work(3)-work(5)+ &
                            4d0 * (work(2)+work(4)+work(6)) )/9d0
            end if
            coeffnodes(1:7,number)=work(1:7)
            do i = 1, m
               coeff(i,number) = dot_product ( work(1:n), B%phi(1:n,i) )
            end do

         end if

      else

!     --- E%ind(number) = 2002,  coefficient is a function of space alone
!                                it is constant per element

         do i = B%kstep, m, B%kstep
            coeff(i,number) = E%user ( E%mst(number) + E%irelem )
         end do

      end if

1000  if ( debug ) then

!     --- Debug information

         call prinrl ( coeffnodes(1,number), m, 'coeffnodes' )
         call prinrl ( coeff(1,number), m, 'coeffintpnts' )
         write(irefwr,*) 'End sepelfillvar1'

      end if  ! ( debug )
      call erclos ( 'sepelfillvar1' )

      end subroutine sepelfillvar1
