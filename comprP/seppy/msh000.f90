      subroutine msh000 ( iinput, krenum, ichois )
! ======================================================================
!
!        programmer    Guus Segal
!        version  9.5  date 30-04-2015 Correction error message
!        version  9.4  date 03-03-2015 Extra debug
!        version  9.3  date 30-11-2014 Change isurf 23 in 99
!
!   copyright (c) 1982-2014  "Ingenieursbureau SEPRA"
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
!     fill some parameters from array iinput
!     If the length of arrays in this subroutine is changed, adapt all
!     calling subroutines, including pipelayers !!!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh
! **********************************************************************
!
!                       MODULES USED
!
      use sepmoduleplot
      use sepmodulekmesh
      use sepmodulemesh
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer iinput(*), krenum(*), ichois

!     ichois         i    Choise parameter JCHOIS in subroutine MESH
!                         Possible values
!                          0 standard method,  input from standard input file
!                            the mesh is created for the first time
!                          1 input is read from the arrays iinput and rinput
!                            the mesh is created for the first time
!                          2 input is read from the standard input file and
!                            stored in the arrays iinput and rinput
!                            the mesh is created for the first time
!                          3 input is read from the arrays iinput and rinput
!                            the co-ordinates of the existing mesh are changed
!                            the topological description remains unchanged
!                          4 input is read from the arrays iinput and rinput.
!                            the coordinates and the topological description of
!                            the mesh are changed, however the input as stored
!                            in iinput and rinput is the same as in a preceding
!                            call of mesh.
!                            The number of points and curves remains the same, &
!                            as well as the number of elements along the curves.
!                            The coordinates of the user points and the curves
!                            are changed by a call of subroutine chanbn and as a
!                            consequence the number elements in the surfaces may
!                            be changed.
!                            This possibility must be preceded by calls of
!                            subroutine mesh and subroutine chanbn
!                          5 See jchois = 4, however, the topological
!                            description remains unchanged
!                          6 Only the arrays iinput and rinput are filled
!                            The subroutine is returned without actually
!                            creating a mesh
!                          7 The input must be stored in the arrays IINPUT and
!                            RINPUT
!                            The subroutine creates the curves and stores the
!                            information of curves and user points via KMESH
!                            After this step the subroutine returns
!                            This possibility is meant for program SEPGEOM
!                         10 The boundary of the mesh has been changed by
!                            subroutine ADAPBN. Information with respect to this
!                            change has been stored in the arrays ICURVS, &
!                            CURVES,IINPUT and RINPUT, which are all part of
!                            array KMESH
!                            The input/output parameters IINPUT and RINPUT are
!                            not used.
!                            The topological description of the mesh may be
!                            changed
!                         11 See ICHOIS=10, however, the topological description
!                            may not be changed
!     iinput         i    Standard integer input array for subroutine MESH
!     krenum         o    Array of length 10 with information about the
!                         renumbering.
!                         Array must be filled as follows:
!                          1: if>=0, only position 1 is filled.
!                             This position contains the parameter method
!                             according to the old SEPRAN definition, i.e.
!                             jrenum = 9: Use default depending on mesh
!                             else:
!                             jrenum = ispecial*100 + method + 3 * ipref
!                             ispecial = 0: standard
!                                        1: ivertx=1, imidp=1, icentr=3
!                             If<0, -jrenum(1) contains the sequence number
!                             of the last entry filled by the user
!                             All other values are filled with defaults
!                          2: method Indication of the renumbering method used
!                                    The algorithms are compared to the original
!                                    numbering:
!                                    0  No renumbering
!                                    1  Cuthill-Mckee algorithm
!                                    2  Sloan algorithm
!                                    3  Both algorithms,  the best is used
!                                    9  Choose depending on mesh
!                             Default: 9
!                          3: ipref  Indication of the criterion to be used
!                                    to find the best algorithm:
!                                    0  Choose the smallest profile
!                                    1  Choose the smallest band width
!                                    2  Use only the renumbered algorithm
!                                       (method=1 or 2)
!                             Default: 0
!                          4: iprint Indicates if information about the
!                                    renumbering must always be printed
!                                    (1) or only depending on ioutp (0)
!                             Default: 0
!                          5: ivertx Priority for vertex nodes (1-3)
!                             Default: 3
!                          6: imidp  Priority for mid point nodes (1-3)
!                             Default: 3
!                          7: icentr Priority for centroids (1-3)
!                             Default: 3
!                          8: ilevel Indicates if renumbering of specia nodes
!                                    must be performed per level (1) or not (0)
!                             Default: 0
!                          9: istmet Indicates how the start of the renumbering
!                                    algorithm must be performed
!                                    0  Automatic
!                                    1  Start user point given
!                                    2  Start curve given
!                                    3  Start surface given
!                             Default: 0
!                         10: istval Value of start point, curve f surface
!                             Default: 1
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer jmaxin, i, icurve, ndefsf, iprec, kinp, k, is, &
              iinp, ifirst, ilast, istart, ndefvl, iinpsav, &
              nchancoor
      integer, allocatable :: ipointsurf(:)
      logical part4

!     i              Counting variable
!     icurve         Type of curve
!     ifirst         First curve (surface) number read
!     iinp           last position used in array iinput
!     iinpsav        Help variable to save the value of iinpx
!     ilast          Last curve (surface) number read
!     ipointsurf     Array containing the relative starting addresses
!                    of all surfaces in the part of iinput referring
!                    to surfaces
!     iprec          Surface number of surface to be copied etc.
!     is             ?
!     istart         Starting value - 1 of plot information in IINPUT
!     jmaxin         if>0 its value indicates where the maximum numbers of
!                    points, curves, surfaces and volumes are stored in iinput
!     k              Counting variable
!     kinp           ?
!     nchancoor      Number of parts with change coordinates
!     ndefsf         Number of curves in surface
!     ndefvl         Number of surfaces in volume
!     part4          If true rinput part 4 is used
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERALLOC        Produce error message in case allocate went wrong
!     ERCLOS         Resets old name of previous subroutine of higher level
!     ERDEALLOC      Produce error message in case deallocate went wrong
!     EROPEN         Produces concatenated name of local subroutine
!     ERRCHR         Put character in error message
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     MSH046         Raise iinp as function of icurve
!     MSH049         Compute starting address of specific surface in array
!                    IINPUT
!     MSH050         Raise iinp as function of ivolm
!     MSHFILLSHAPE1  Fill shape numbers of a surface that has not been
!                    filled yet
!     PRININ         print 1d integer array
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!      15   jmark has in correct value
!     119   The number of curves is incorrect
!     120   The number of nodal points or the nodal point number read is
!           incorrect
!     121   number of line elements<=0 whereas ndim = 1
!     150   The number of surfaces is incorrect
!     154   The number of volumes is incorrect
!     298   code for generation of volumes incorrect
!     379   Incorrect type of surface generator found
!    1001   Parameter notopol out of range
!    1240   The number of times the mesh should be refined is out of range
!    1328   The type number read after transform has incorrect value
!    1947   TRANSF or refine not permitted for changing mesh
!    2147   Parameter out of range
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
      call eropen ( 'msh000' )
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from msh000'
         call prinin ( iinput, 21, 'iinput' )
  1      format ( a, 1x, (10i6) )

      end if
      if ( ierror/=0 ) go to 1000

      part4 = .false.
      ndim = iinput(1)
      nuspnt = iinput(2)
      ncurvs = iinput(3)
      nsurfs = iinput(4)
      nvolms = iinput(5)
      iplot  = iinput(9)
      jrenum = iinput(10)
      if ( debug ) then
         write(irefwr,1) 'ndim, nuspnt, ncurvs, nsurfs, nvolms', &
                          ndim, nuspnt, ncurvs, nsurfs, nvolms
         write(irefwr,1) 'iplot, jrenum', iplot, jrenum
      end if  ! ( debug )

      allocate ( ipointsurf(max(1,nsurfs)), stat = error )
      if ( error/=0 ) &
         call eralloc (  error, max(1,nsurfs), 'ipointsurf' )
      if ( ierror/=0 ) go to 1000

      if ( jrenum>=0 ) then
         krenum(1) = jrenum
      else
         krenum(1) = -10
         do i = 1, 9
            krenum(1+i) = iinput(-iinput(10)-6+i)
         end do
      end if
      jmaxin = iinput(13)-5
      if ( iinput(14)==0 ) then
          iinp9 = 0
          nspec = 0
      else
          iinp9 = iinput(14)-5
          nspec = 1
      end if
      iinp8 = iinp9
      if ( jmaxin>0 ) then
         maxpnt = iinput(jmaxin)
         maxcrv = iinput(jmaxin+1)
         maxsrf = iinput(jmaxin+2)
         maxvlm = iinput(jmaxin+3)
         iinp8 = max(iinp8,jmaxin+3)
      else
         maxpnt = 10000
         maxcrv = 1000
         maxsrf = 1000
         maxvlm = 500
      end if

      if ( debug ) then
         write(irefwr,1) 'maxpnt, maxcrv, maxsrf, maxvlm', &
                          maxpnt, maxcrv, maxsrf, maxvlm
         write(irefwr,1) 'iinp8, iinp9',iinp8, iinp9
      end if  ! ( debug )

      iinp10 = 0
      iinp11 = 0
      irefin = 0
      icheck = 0
      itrans = 0
      notopol = 0
      maxlencur = 9
      maxlensur =10
      maxlenvol = 0
      iinp12 = 0
      iinp13 = 0
      iinp14 = 0
      iinp15 = 0
      irnp5 = 0
      nintprop = 0
      nrealprop = 0
      iwritfaces = 0
      itype_scale = 0
      iscale = 0
      iparallel = 0
      nelgrpdum = 0
      if ( iinput(15)>0 ) then

!     --- iinput(15)>0, hence extra information

         iinp10 = iinput(15)-5
         irefin = iinput(iinp10)
         itrans = iinput(iinp10+1)
         icheck = iinput(iinp10+2)
         notopol = iinput(iinp10+3)
         iinp11 = max(0,iinput(iinp10+4)-5)
         iinp12 = max(0,iinput(iinp10+5)-3)
         iinp13 = max(0,iinput(iinp10+6)-5)
         iinp14 = max(0,iinput(iinp10+7)-5)
         iwritfaces = iinput(iinp10+8)
         iscale = iinput(iinp10+9)
         itype_scale = iinput(iinp10+10)
         iparallel = iinput(iinp10+11)
         iinp15 = max(0,iinput(iinp10+12)-5)

         if ( debug ) then

!        --- Debug information

            write(irefwr,1) 'iinp10, irefin, itrans, icheck, notopol', &
                             iinp10, irefin, itrans, icheck, notopol
            write(irefwr,1) 'iinp11, iinp12, iinp13, iinp14', &
                             iinp11, iinp12, iinp13, iinp14
            write(irefwr,1) 'iwritfaces, iscale, itype_scale', &
                             iwritfaces, iscale, itype_scale
            write(irefwr,1) 'iparallel, iinp15', &
                             iparallel, iinp15

         end if

         if ( iscale<0 .or. iscale>1 ) then
            call errchr ( 'iscale', 1 )
            call errint ( iscale, 1 )
            call errint ( 0, 2 )
            call errint ( 1, 3 )
            call errsub ( 2147, 3, 0, 1 )
         end if
         if ( itype_scale<0 .or. itype_scale>0 ) then
            call errchr ( 'itype_scale', 1 )
            call errint ( itype_scale, 1 )
            call errint ( 0, 2 )
            call errint ( 1, 3 )
            call errsub ( 2147, 3, 0, 1 )
         end if
         if ( iparallel<0 .or. mod(iparallel,100)>7 .or. &
              iparallel/100>1024 ) then
            call errchr ( 'iparallel', 1 )
            call errint ( iparallel, 1 )
            call errint ( 0, 2 )
            call errint ( 7, 3 )
            call errsub ( 2147, 3, 0, 1 )
         end if
         iinp8 = max(iinp8,iinp10+19)
         if ( irefin<0 .or. irefin>=100 ) then
            call errint ( irefin, 1 )
            call errsub ( 1240, 1, 0, 0 )
         end if
         if ( itrans/=0 .and. (itrans<4 .or. itrans>7) ) then
            call errint ( itrans, 1 )
            call errsub ( 1328, 1, 0, 0 )
         end if
         if ( notopol<0 .or. notopol>1 ) then
            call errint ( notopol, 1 )
            call errsub ( 1001, 1, 0, 0 )
         end if
         if ( debug ) write(irefwr,1) 'iinp11, iinp8', iinp11, iinp8
         if ( iinp11>0 ) then

!        --- interface elements found

            iinp8 = iinp11+4*(iinput(iinp11+1)+iinput(iinp11+2))
            if ( debug ) then
               write(irefwr,1) 'iinp8', iinp8
               write(irefwr,1) 'interface_curves, interfaces_surfs', &
                                iinput(iinp11+1), iinput(iinp11+2)
            end if  ! ( debug )

         end if
         if ( iinp12>0 ) then

!        --- properties found

            nintprop = iinput(iinp12-2)
            nrealprop = iinput(iinp12-1)
            if ( debug ) &
               write(irefwr,1) 'nintprop, nrealprop', &
                                nintprop, nrealprop

         end if
         if ( iinp15>0 ) then

!        --- dummy elements found

            nelgrpdum = iinput(iinp15)
            iinp8 = iinp15+2*nelgrpdum

         end if

      end if

      if ( iplot>0 ) then

!     --- iplot>0, old situation

         numsub = iplot/10
         jmark = iplot-10*numsub-1
         iorien = 1

      else if ( iplot<0 ) then

!     --- iplot<0, new situation

         istart = -iinput(9)-6
         jmark  = iinput(istart+1)
         numsub = iinput(istart+2)
         jcurvp = iinput(istart+3)
         jcurvn = iinput(istart+4)
         juserp = iinput(istart+5)
         dcolor = iinput(istart+6)
         i3d = iinput(istart+7)
         irenpl = iinput(istart+8)
         nomesh = iinput(istart+9)
         iorien = iinput(istart+10)
         iinp8 = max(iinp8,istart+10)
         if ( debug ) then

!        --- Debug information

            write(irefwr,1) 'istart, jmark, numsub, dcolor', &
                             istart, jmark, numsub, dcolor
            write(irefwr,1) 'iorien, iinp8', iorien, iinp8
            call prinin ( iinput(istart), 11, 'iinput_plot' )

         end if  ! ( debug )

      else

!     --- iplot = 0

         numsub = 0
         jmark  = 0
         iorien = 0

      end if
      if ( ichois>=3 .and. ichois<=5 .or. ichois>=10 ) then

!     --- mesh must be changed

         irenpl = 0
         iplot = 0
         if ( itrans/=0 .or. irefin/=0 ) then

!        --- Refine not implemented

            call errsub ( 1947, 0, 0, 0 )

         end if

      end if
      if ( jmark==-1) jmark = 0
      if ( jmark>5 ) then
         call errint ( jmark, 1 )
         call errint ( 5, 2 )
         call errsub ( 15, 2, 0, 0 )
      end if
      if ( ncurvs<=0 .or. ncurvs>maxcrv ) then
         call errint ( ncurvs, 1 )
         call errint ( maxcrv, 2 )
         call errsub ( 119, 2, 0, 0 )
      end if
      if ( nuspnt<=0 .or. nuspnt>maxpnt ) then
         call errint ( nuspnt, 1 )
         call errint ( maxpnt, 2 )
         call errsub ( 120, 2, 0, 0 )
      end if
      if ( nsurfs<0 .or. nsurfs>maxsrf ) then
         call errint ( nsurfs, 1 )
         call errint ( maxsrf, 2 )
         call errsub ( 150, 2, 0, 0 )
      end if
      if ( nvolms<0 .or. nvolms>maxvlm ) then
         call errint ( nvolms, 1 )
         call errint ( maxvlm, 2 )
         call errsub ( 154, 2, 0, 0 )
      end if
      if ( ierror/=0 ) go to 1000

      nlines = iinput(6)
      nsurel = iinput(7)
      nvolel = iinput(8)
      irnp1 = 3
      irnp2 = (ndim+jcoars)*nuspnt+2*jcoars+irnp1
      irnp3 = irnp2+3*ncurvs-1+4*nsurfs
      irnp4 = irnp3+1
      if ( ndim==3 .and. i3d/=0 ) irnp3 = irnp3+3
      if ( debug ) then
         write(irefwr,1) 'ndim, nuspnt, ncurvs, nsurfs, nvolms', &
                          ndim, nuspnt, ncurvs, nsurfs, nvolms
         write(irefwr,1) 'nlines, nsurel, nvolel, jcoars', &
                          nlines, nsurel, nvolel, jcoars
      end if  ! ( debug )

!     --- When these starting positions are changed, also consider chanbn!!!!

      iinp1 = 16
      iinp3 = iinp1
      do i = 1,ncurvs
         icurve = iinput(iinp3)
         iinpsav = iinp3
         call msh046(iinput,iinp3,icurve)

         if ( debug ) then
            write(irefwr,1) 'i, icurve', i, icurve
            call prinin ( iinput(iinpsav), iinp3-iinpsav, 'iinput' )
         end if  ! ( debug )
         maxlencur = max ( maxlencur, iinp3-iinpsav )
      end do
      if ( debug ) write(irefwr,1) 'iinp1, iinp3', iinp1, iinp3

      iinp5 = iinp3

!     --- surface elements

      do i = 1, nsurfs
         ipointsurf(i) = iinp5
         isurf = iinput(iinp5)
         iinpsav = iinp5
         call msh049 ( iinput, iinp5, isurf )
         maxlensur = max ( maxlensur, iinp5-iinpsav )
      end do
      iinpsav = iinp5

      if ( debug ) then
         write(irefwr,1) 'iinp5sav', iinpsav
         call prinin ( ipointsurf, nsurfs, 'ipointsurf' )
      end if  ! ( debug )

      do i = 1, nsurfs
         iinp5 = ipointsurf(i)
         isurf = iinput(iinp5)

         if ( isurf<0 .or. isurf==6 .or. isurf>24 .and. isurf/=99 ) then

!        --- Error 379:  Type of surface incorrect

            call errint ( isurf, 1 )
            call errint ( i, 2 )
            call errsub ( 379, 2, 0, 0)
            go to 1000

         end if

         if ( isurf==0 ) then
            go to 400
         else if ( isurf==19 ) then
            part4 = .true.
         end if
         ishape = iinput(iinp5+1)
         if ( debug ) &
            write(irefwr,1) 'iinp5, isurf, ishape', iinp5, isurf, ishape

         if ( ishape==0 ) then

!        --- ishape = 0   either error or ishape must be copied

            call  mshfillshape1 ( iinput, ipointsurf, isurf, iinp5, i )

         end if

400   end do
      iinp5 = iinpsav

      iinp7 = iinp5
      if ( debug ) write(irefwr,1) 'iinp5', iinp5

!     --- volume elements

      do i = 1, nvolms

         ivol = iinput(iinp7)
         if ( ivol==0 ) then
            iinp7 = iinp7 + 1
            go to 500
         end if
         if ( debug ) write(irefwr,1) 'i, ivol, iinp7', i, ivol, iinp7

         if ( ivol<0 .or. ivol>11 ) then
            call errint ( ivol, 1 )
            call errsub ( 298, 1, 0, 0 )
            go to 1000
         end if

         if ( ivol>5 .and. ivol<9 ) then

!        --- Place shape number by copying from volume to be
!            translated, rotated or reflected

            ndefvl = iinput(iinp7+2)
            iprec  = iinput(iinp7+3+ndefvl)
            kinp   = iinp5
            do k = 1, iprec-1
                is = iinput(kinp)
                call msh050(iinput, kinp, is)
            end do
            ishape = iinput( kinp+1 )
            iinput(iinp7+1) = ishape

         end if

         iinpsav = iinp7
         call msh050 ( iinput, iinp7, ivol )
         maxlenvol = max ( maxlenvol, iinp7-iinpsav )

500   end do

      iinp2 = iinp7
      iinp4 = iinp2+4*nlines
      iinp6 = iinp4+4*nsurel
      iinp7 = iinp6+4*nvolel
      if ( debug ) &
         write(irefwr,1) 'iinp2, iinp4, iinp6, iinp7', &
                          iinp2, iinp4, iinp6, iinp7

      if ( nlines<=0 .and. nsurel<=0 .and. nvolel<=0 ) &
         call errsub ( 121, 0, 0, 0 )

      iinp8 = max(iinp8,iinp7-1)
      if ( iinput(12)/=0 ) then

!     --- Connected elements found

         length = iinp7-1 + 3*iinput(iinp7) + 4*iinput(iinp7+1)
         ndefsf = iinput(iinp7+2)
         do i = 1, ndefsf
            length = length + 5 + max(iinput(length+5),0)
         end do
         iinp8 = max(iinp8,length)

      end if

      iinp8 = max ( iinp10+1, iinp8 )
      if ( iinput(9)<0 ) iinp8 = max ( iinp8, -iinput(9) )
      if ( iinput(14)>0 ) iinp8 = max ( iinp8, iinput(14)-1 )

      nsurvl = 0
      iinp = iinp4
      do i = 1,nsurel
         ifirst = iinput(iinp+2)
         ilast = iinput(iinp+3)
         nsurvl = nsurvl+ilast-ifirst+1
         iinp = iinp+4
      end do
      iinp = iinp6
      do i = 1,nvolel
         ifirst = iinput(iinp+2)
         ilast = iinput(iinp+3)
         nsurvl = nsurvl+ilast-ifirst+1
         iinp = iinp+4
      end do

!     --- Compute nelgrp for standard elements

      if ( debug ) &
         write(irefwr,1) 'nvolel, nsurel, nlines, nelgrp', &
                          nvolel, nsurel, nlines, nelgrp

      if ( nvolel>0 ) then

!     --- nvolel>0

         nelgrp = iinput(iinp6+1+(nvolel-1)*4)

      else if ( nsurel>0 ) then

!     --- nsurel>0, nvolel = 0

         nelgrp = iinput(iinp4+1+(nsurel-1)*4)

      else

!     --- nlines>0, nvolel = nsurel = 0

         nelgrp = iinput(iinp2+1+(nlines-1)*4)

      end if
      if ( debug ) write(irefwr,1) 'nelgrp, nelgrpdum', nelgrp, nelgrpdum
      nelgrp = nelgrp+nelgrpdum

      if ( jmaxin>0 ) iinp8 = max ( iinp8, jmaxin+3 )
      if ( part4 ) irnp3 = irnp3+4*ncurvs
      if ( iinp12>0 ) then

!     --- properties found

         nintprop = iinput(iinp12-2)
         iinp8 = max(iinp8,iinp12+nintprop*nelgrp)
         nrealprop = iinput(iinp12-1)
         irnp5 = irnp3+1
         irnp3 = irnp3+nrealprop*nelgrp

         if ( debug ) &
            write(irefwr,1) 'nintprop, nrealprop, iinp8, irnp3, irnp5', &
                             nintprop, nrealprop, iinp8, irnp3, irnp5

      end if
      if ( iinp13>0 ) then

!     --- change_coordinates found

         nchancoor = iinput(iinp13)
         iinp8 = max(iinp8,iinp13+nchancoor*4)

      end if
      if ( iinp14>0 ) then

!     --- obstacles found

         numobst = iinput(iinp14)
         iinp8 = max(iinp8,iinp14+numobst*2)

      end if
      if ( .not. part4 ) irnp4 = 0

      deallocate ( ipointsurf, stat = error )
      if ( error/=0 ) call erdealloc ( error, 'ipointsurf' )

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,1) 'iinp1, iinp2, iinp3, iinp4, iinp5', &
                          iinp1, iinp2, iinp3, iinp4, iinp5
         write(irefwr,1) 'iinp6, iinp7, iinp8, iinp9, iinp10', &
                          iinp6, iinp7, iinp8, iinp9, iinp10
         write(irefwr,1) 'iinp11, iinp12, iinp13, iinp14, iinp15', &
                          iinp11, iinp12, iinp13, iinp14, iinp15
         write(irefwr,1) 'irnp1, irnp2, irnp3, irnp4, irnp5', &
                          irnp1, irnp2, irnp3, irnp4, irnp5
         call prinin ( iinput, iinp8, 'iinput' )
         write(irefwr,*) 'End msh000'

      end if
      call erclos ( 'msh000' )

      end subroutine msh000
