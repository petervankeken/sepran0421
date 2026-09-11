      subroutine chsloa ( mstart, level, nstats, istmet, ifirst, ilast )
! ======================================================================
!
!        programmer    Guus Segal
!        version  4.4  date 11-11-2014 Extra debug
!        version  4.3  date 13-02-2008 Extra debug
!        version  4.2  date 17-04-2006 Fortran 90
!
!   copyright (c) 1987-2014  "Ingenieursbureau SEPRA"
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
!    Renumbering algorithm according to Sloan
!    See: Int Journal of Num. Meth. in Engng. Vol. 23,  pp. 239-251  (1986)
! **********************************************************************
!
!                       KEYWORDS
!
!     renumbering
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulework
      use sepmodulekmesh
      use sepmodulecomio
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer, target :: level(*)
      integer :: mstart,  nstats, istmet, ifirst, ilast

!     ifirst  i    The nodes from ifirst must be renumbered
!     ilast   i    The nodes until ilast must be renumbered
!     istmet  i    Indicates the method for computing the start for renumbering
!                  0:  default
!                  1:  user point
!                  2:  curve
!                  3:  surface
!     level   i    Integer work array of length npoint to store the levels
!     mstart  i    Number of points in set kmeshj at input
!     nstats  o    Number of points in kmeshj at output
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, mset, iw1, iw2, kpoint, maxdeg, ideg, &
              inum, ndist, ipoint, nprio, k1, l, jdist, nbr, nbr1, j, &
              k2, k, k3, ncount, ilow, ihigh, ineigh, inodp, jnodp, &
              in, ih, ihigh1, jstart, maxwid, iend, nlev, nlpoint ! , &
!             itemp(npoint), iprio(npoint), idist(npoint)
      integer, allocatable :: itemp(:),iprio(:),idist(:)

!     i              Counting variable
!     ideg
!     idist          integer work array of length npoint in which the distance
!                    of all points to the end point is stored
!     iend           End point of rooted level structure
!     ih             Help variable to store a constant temporarily
!     ihigh
!     ihigh1         ?
!     ilow
!     in
!     ineigh         set of new neighbours
!     inodp          Nodal point number
!     inum
!     ipoint
!     iprio          integer work array of length npoint positions to store the
!                    priorities of each node
!     itemp          Work array
!     iw1            ?
!     iw2            ?
!     j              Counting variable
!     jdist          ?
!     jnodp          ?
!     jstart         Point in end level
!     k              Counting variable
!     k1
!     k2
!     k3             ?
!     kpoint         Largest node number that has been renumbered
!     l              General loop variable
!     maxdeg         ?
!     maxwid         Maximum number of elements in end level if iend is start
!                    point
!     mset
!     nbr            number of neighbours of nodal point
!     nbr1           ?
!     ncount         Counter
!     ndist          ?
!     nlev           Number of elements in end level
!     nlpoint        Largest node number that may been renumbered in this
!                    connected set
!     nprio          ?
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     CHROLS         Generate root level structure corresponding to root istart
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     PRININ         print 1d integer array
!     SEPNEXTLEVEL   Find neighbors corresponding to a set of nodes
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
      allocate(idist(npoint),iprio(npoint),itemp(npoint))
      call eropen ( 'chsloa' )
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from chsloa'

         call prinin ( kmeshj, npoint, 'kmeshj' )

  1      format ( a, 1x, (10i8) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

!     --- Sort starting set and store in kmeshj

      do i = 1, ifirst-1
         marknodes(i) = -10
      end do
      do i = ifirst, ilast
         marknodes(i) = -3
      end do

      if ( debug ) then

!     --- debug statements

         write(irefwr,*) 'in chsloa'
         write(irefwr,1) 'ifirst, ilast', ifirst, ilast
         write(irefwr,1) 'istmet,mstart', istmet,mstart

      end if

!     --- Set weights;  values are chosen according to Sloan

      iw1 = 2
      iw2 = 1
      if ( istmet==0 .and. mstart>1 ) then

!     --- The starting set is reduced to 1 point if istmet = 0

         mstart = (mstart+2)/2+ifirst-1
         do i = ifirst, mstart
            itemp(i) = kmeshj(i)
         end do
         maxwid = npoint

         if ( debug ) then

!        --- debug statements

            write(irefwr,1) 'ifirst, mstart, npoint, maxwid', &
                             ifirst, mstart, npoint, maxwid
            call prinin ( itemp(ifirst), mstart-ifirst+1, 'itemp' )

         end if

         do i = ifirst, mstart

            jstart = itemp(i)
            call chrols ( jstart, level, nlev, &
                          npoint, kpoint, iw2, idist )

            if ( debug ) then

!           --- debug statements

               write(irefwr,1) 'i, jstart, nlev, maxwid', &
                                i, jstart, nlev, maxwid

               write(irefwr,1) 'itemp(i)', itemp(i)
            end if

            if ( nlev<maxwid ) then
               iend = itemp(i)
               maxwid = nlev
            end if

         end do
         kmeshj(ifirst) = iend
         mstart = ifirst

      else if ( istmet==0 ) then

         mstart = ifirst

      end if
      if ( debug ) &
         write(irefwr,1) 'istmet, mstart, ifirst', &
                          istmet, mstart, ifirst

      do i = ifirst, mstart
         marknodes(kmeshj(i)) = 0
         level(i-ifirst+1) = kmeshj(i)
      end do

      if ( debug ) then

!     --- debug statements

         call prinin ( marknodes, npoint, 'marknodes' )
         call prinin ( level, mstart-ifirst+1, 'level' )
         call prinin ( kmeshj(ifirst), mstart-ifirst+1, 'start set' )

      end if

!     --- Find end level in order to set the distances

      mset = mstart-ifirst+1
      call chrols ( -1, level, mset, &
                    ilast, kpoint, 0, idist )
      kpoint = kpoint+ifirst-1

      if ( debug ) then

!     --- debug statements

         write(irefwr,1) 'mset, mstart, ifirst, kpoint', &
                          mset, mstart, ifirst, kpoint

      end if

!     --- iend is part of the set level.
!         The point of this set that ends with the smallest starting level
!         is used as iend

      nlpoint = kpoint
      if ( mset>1 ) then

!     --- mset>1

         mset = (mset+2)/2
         do i = 1, mset
            itemp(i) = level(i)
         end do
         maxwid = nlpoint
         do i = 1, mset
            jstart = itemp(i)
            if ( debug ) write(irefwr,1) 'i, jstart', i, jstart
            call chrols ( jstart, level, nlev, &
                          ilast, kpoint, iw2, idist )
            if ( nlev<maxwid ) then
               iend = itemp(i)
               maxwid = nlev
            end if
         end do

         level(1) = iend
         mset = 1
         if ( debug ) write(irefwr,1) 'mset, iend', mset, iend

      end if

!     --- Fill priority array iprio
!         First the rooted level structure corresponding to iend is build and
!         the part due to the distances is set into iprio

      call chrols ( -1, level, mset, ilast, kpoint, iw2, idist )
      if ( debug ) then
         call prinin ( level, nlevel, 'level' )
         call prinin ( idist, npoint, 'idist' )
         write(irefwr,1) 'ifirst, nlpoint, mset, ilast', &
                          ifirst, nlpoint, mset, ilast
      end if  ! ( debug )

!     --- Next iw1 * degree(i) is subtracted and maximum degree is computed

      maxdeg = 0
      do i = ifirst, nlpoint
         ideg = kmshe1(i+1)-kmshe1(i)
         iprio(i) = idist(i)-iw1*ideg
         maxdeg = max(maxdeg, ideg)
      end do
      if ( debug ) &
         call prinin ( iprio(ifirst), nlpoint-ifirst+1, 'iprio/1' )

!     --- Finally maxdeg * iw1 is added to iprio

      do i = ifirst, nlpoint
         iprio(i) = iprio(i)+iw1*maxdeg
      end do
      if ( debug ) &
         call prinin ( iprio(ifirst), nlpoint-ifirst+1, 'iprio/2' )

!     --- Set status of all nodes inactive

      do i = ifirst, nlpoint
         if ( marknodes(i)>-10 ) marknodes(i) = 0
      end do
      do i = ifirst, mstart
         marknodes(kmeshj(i)) = -10
      end do
      if ( debug ) then
         write(irefwr,1) 'ifirst, nlpoint, mstart', &
                          ifirst, nlpoint, mstart
         call prinin ( marknodes, npoint, 'marknodes' )
         call prinin ( kmeshj, npoint, 'kmeshj' )
      end if  ! ( debug )

!     --- nstats is equal to the number of inactive points
!         Initialize set of eligible nodes E  (array level)

      nlevel = mstart
      nstats = nlpoint - mstart

!     --- Start node is preactive

      lastlevel => kmeshj(1:mstart)
      nextlevel => level(1:npoint)
      call sepnextlevel ( -1, nlevel )
      do i = ifirst, nlpoint
         if ( marknodes(i)>-10 ) marknodes(i) = 0
      end do
      do i = 1, nlevel
         marknodes(level(i)) = 1
      end do
      nstats = nstats - nlevel
      if ( debug ) then
         write(irefwr,1) 'mstart, nlevel', mstart, nlevel
         call prinin ( marknodes, npoint, 'marknodes' )
         call prinin ( level, nlevel, 'level' )
      end if  ! ( debug )

!     --- kmeshj is empty  (last position used:  inum)

      inum = mstart

!     --- While loop as long as set E not empty

      ndist = nlpoint
      do ipoint = ifirst, nlpoint
         if ( nlevel==0 ) go to 900

!        --- Find node i in E with highest priority

         nprio = -1
         kpoint = -1
         do k1 = 1, nlevel
            l = level(k1)
            jdist = idist(l)
            if ( iprio(l)>nprio ) then
               nprio = iprio(l)
               i = l
               nbr = kmshe1(i+1)-kmshe1(i)
               kpoint = k1
               ndist = min(ndist, jdist)
            else if (ndist>jdist .and. iprio(l)==nprio) then
               nbr1 = kmshe1(l+1)-kmshe1(l)
               if ( nbr1<nbr ) then
                  nbr = nbr1
                  i = l
                  kpoint = k1
                  ndist = min(ndist, jdist)
               end if
            end if
         end do

!        --- Check if there are any points with positive priority
!            If not all points are preactive, active or post active
!            In that case the last points are ordered according to a special
!            algorithm

         if ( idist(i)<=0 .or. kpoint<0 ) goto 400

!        --- Remove node i from set E

         do k1 = kpoint+1, nlevel
            level(k1-1) = level(k1)
         end do
         ndist = idist(i)
         nlevel = nlevel-1
         if ( marknodes(i)==1 ) then

!        --- i is preactive, find all neighbours j of i

            do k1 = kmshe1(i)+1, kmshe1(i+1)
               j = kmshe2(k1)
               iprio(j) = iprio(j)+iw1
               if ( marknodes(j)==0 ) then

!              --- j is inactive, make it preactive and put in set E

                  marknodes(j) = 1
                  nstats = nstats-1
                  if ( nstats==0 ) then

!                 --- nstats = 0   There are no active nodes anymore
!                                  All nodes which are preactive get the
!                                  negative priority -10*npoint

                     do k2 = 1, nlpoint
                        if ( marknodes(k2)==1 ) then
                           iprio(k2) = -10*npoint
                        end if
                     end do
                  end if
                  nlevel = nlevel+1
                  level(nlevel) = j
               end if
            end do
         end if

!        --- Put node i in renumbered set and make i postactive

         inum = inum+1
         kmeshj(inum) = i
         marknodes(i) = -10

         if ( debug ) then

!        --- debug statements

            write(irefwr,1) 'i, inum', i, inum

         end if

!        --- Loop over all neighbours of i

         do k1 = kmshe1(i)+1, kmshe1(i+1)
            j = kmshe2(k1)
            if ( marknodes(j)==1 ) then

!           --- j is preactive,  raise priority and make j active

               iprio(j) = iprio(j)+iw1
               marknodes(j) = 2

!              --- Consider all neighbours of j

               do k2 = kmshe1(j)+1, kmshe1(j+1)
                  k = kmshe2(k2)
                  if ( marknodes(k)==1 .or. marknodes(k)==2 ) then

!                 --- k is active or preactive, raise priority

                     iprio(k) = iprio(k)+iw1

                  else if ( marknodes(k)==0 ) then

!                 --- k is inactive, raise priority, make k preactive and
!                     put k in set E

                     iprio(k) = iprio(k)+iw1
                     nlevel = nlevel+1
                     level(nlevel) = k
                     marknodes(k) = 1
                     nstats = nstats-1
                     if ( nstats==0 ) then

!                    --- nstats = 0   There are no active nodes anymore
!                                     All nodes which are preactive get the
!                                     negative priority -10*npoint

                        do k3 = 1, ilast
                           if ( marknodes(k3)==1 ) then
                              iprio(k3) = -10*npoint
                           end if
                        end do
                     end if
                  end if
               end do
            end if
         end do
400   end do

!     --- Only points with negative priority are left
!         Use special algortihm to sort remaining points
!         These points are sorted such that points with smallest new neighbour
!         number get lowest nodal point number
!         First store new numbers in array istat
!         Clear marknodes

      do i = ifirst, nlpoint
         marknodes(i) = 0
      end do

!     --- Fill marknodes for numbered points

      do i = 1, inum
         marknodes(kmeshj(i)) = i
      end do

!     --- Put set of non-numbered points in iprio and set corresponding position
!         in marknodes equal to npoint+1

      if ( debug ) then
         call prinin ( kmeshj, inum, 'kmeshj' )
         call prinin ( marknodes, npoint, 'marknodes' )
      end if  ! ( debug )

      ncount = 0
      do i = ifirst, nlpoint
         if ( marknodes(i)==0 ) then
            ncount = ncount+1
            iprio(ncount) = i
            marknodes(i) = npoint+1
         end if
      end do

      if ( debug ) then
         write(irefwr,1) 'ncount, ifirst, nlpoint', &
                          ncount, ifirst, nlpoint
         call prinin ( iprio, ncount, 'iprio' )
         call prinin ( marknodes, npoint, 'marknodes' )
      end if  ! ( debug )

!     --- Sort array iprio

      do i = 1, ncount-1
         ilow = nlpoint+1
         ihigh = ilow
         ineigh = 0
         do j = i, ncount
            inodp = iprio(j)

!           --- Check all neighbours of inodp
!               Find if there is one with smaller number than ilow

            do k1 = kmshe1(inodp)+1, kmshe1(inodp+1)
               jnodp = kmshe2(k1)
               if ( marknodes(jnodp)<ilow ) then

!              --- marknodes(jnodp)<ilow, so new number found

                  ilow = marknodes(jnodp)
                  ineigh = j
                  ihigh = 0
               else if ( marknodes(jnodp)==ilow ) then

!              --- marknodes(jnodp) = ilow
!                  Point "j" and point "ineigh" are connected to the same
!                  "low" number.
!                  Check which of the two is related to the highest number
!                  and reject that one

                  if ( ihigh==0 .and. ineigh>0 ) then

!                 --- ihigh = 0, so highest number corresponding to "ineigh"
!                     has not yet been found

                     in = iprio(ineigh)
                     do k2 = kmshe1(in)+1, kmshe1(in+1)
                        ih = marknodes(kmshe2(k2))
                        if ( ih<=npoint ) ihigh = max ( ihigh, ih )
                     end do
                  end if

!                 --- Find highest number ihigh1 corresponding to point "j"

                  ihigh1 = 0
                  do k2 = kmshe1(inodp)+1, kmshe1(inodp+1)
                     ih = marknodes(kmshe2(k2))
                     if ( ih<=npoint ) ihigh1 = max ( ihigh1, ih )
                  end do
                  if ( ihigh>ihigh1 ) then

!                 --- ihigh>ihigh1:  interchange both points

                     ineigh = j
                     ihigh = ihigh1
                  end if
               end if
            end do
         end do

!        --- ineigh gives position in iprio with smallest neighbour number
!            Replace ineigh with position 1 in iprio and store iprio(ineigh)
!            in array kmeshj

         inum = inum + 1
         kmeshj(inum) = iprio(ineigh)
         marknodes(kmeshj(inum)) = inum
         iprio(ineigh) = iprio(i)

         if ( debug ) then

!        --- debug statements

            write(irefwr,1) 'inum, kmeshj(inum)', inum, kmeshj(inum)

         end if

      end do

!     --- Put last position in kmeshj

      kmeshj(inum+1) = iprio(ncount)
      inum = inum+1
 900  nstats = inum

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,1) 'ncount, inum, nstats', ncount, inum, nstats
         write(irefwr,*) 'End chsloa'

      end if  ! ( debug )
      call erclos ( 'chsloa' )
      if (allocated(idist)) deallocate(idist)
      if (allocated(iprio)) deallocate(iprio)
      if (allocated(itemp)) deallocate(itemp)
      end
