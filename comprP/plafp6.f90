      subroutine plafp6 ( ichois, x, y )
      use control
! ======================================================================
!
!        programmer     Guus Segal
!        version 3.35   date   02-09-2013 Increase max. number of plots
!        version 3.33   date   02-02-2004 Adaptation for /tmp
!        version 3.32   date   13-06-2002 New call to plpos6
!
!   copyright (c) 1986-2013  "Ingenieursbureau SEPRA"
!   permission to copy or distribute this software or documentation
!   in hard copy or soft copy granted only by written license
!   obtained from "Ingenieursbureau SEPRA".
!   all rights reserved. no part of this publication may be reproduced,
!   stored in a retrieval system ( e.g., in memory, disk, or core)
!   or be transmitted by any means, electronic, mechanical, photocopy,
!   recording, or otherwise, without written permission from the
!   publisher.
!
!
! ********************************************************************
!
!                       DESCRIPTION
!
!    this subroutine contains the machine dependent calls of plotting
!    subroutines, except for texts, special symbols and colours
!    machine dependent subroutine
!    starbase/ hp-gl/ tektronix/ plot10/ apollo/ calcomp/ gks /cgi
!    PC with VGA or EGA
!
! **********************************************************************
!
!                       KEYWORDS
!
!    machine_dependent
!    plot
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulecpack
      use sepmodulemach
      use sepmoduleplot
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer ichois
      double precision x, y

!    ichois   i   The use of the various parameters depends on ichois.
!                 Possibilities:
!                 1 Call plot start subroutine
!                   When X=1 small plotting paper is used, when X=2 wide
!                   plotting paper is used. (not supported in some versions)
!                   The length of the plot is given in Y (centimeters)
!                 2 Plot with pen down to point (x,y)
!                 3 Plot with pen up to point (x,y)
!                 4 Change origin (will disappear in next major revision!)
!                   (x,y) are the values to be added to the current origin
!                 5 Close plot data set
!                   When X=-1 the number for the next (w or u) output file
!                   will be reset to 1 again
!    x        i   See ichois
!    y        i   See ichois
!
! **********************************************************************
!
!                       LOCAL PARAMETERS

      integer lenhome, leng, iref,ireftk, i, lastps
      real sx, sy
      character (len=80) filenm, name
      logical check, isopen
      save leng, isopen, filenm

!    check        Indicates if files exists (true) or not (false)
!    filenm       character field to build filename in
!    iref         Reference number for plot files
!    ireftk       Reference number for Tektronix files
!    isopen       File iref is open
!    leng       Number of characters in filenm or namebn
!    lastps       Number of characters in name
!    lenhome      Number of characters in the name sphome
!    name         character field to build file name in
!    sx           new x coordinate, after addition of new origin
!    sy           new y coordinate, after addition of new origin
! ********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS  Resets old name of previous subroutine of higher level
!     EROPEN  Produces concatenated name of local subroutine
!     ERRCHR  Put character in error message
!     ERRINT  Put integer in error message
!     ERRSUB  Error messages
!     ERRWAR  Warnings
!     PLCOLI  Initialization of color tables
!     PLCOL0  Routine to set default line/marker color
!     PLHGL6  Perform the actions in case of a HPGL plot
!     PLPOS6  Perform the actions in case of a Postscript plot
!     PLSTB6  Perform the actions in case of a StarBase plot
!     PLVGA6  Perform the actions in case of a plot to a VGA screen
!             in combination with the Salford compiler
!     PLTEK0  Routines for Tektronix 4010 plotting
!     PLTEK1  Routines for Tektronix 4010 plotting
!     PLTEK2  Routines for Tektronix 4010 plotting
! ********************************************************************
!
!                       I/O
!
! ********************************************************************
!
!                       ERROR MESSAGES
!
!     294   Picture missing
!    1996   Error in opening plot file
!    1997   Choice not supported
!    1998   Device not supported
!    2285   Plot file is not open
! **********************************************************************
!
!                       PSEUDO CODE
!
!                 Machine-dependent actions
!
!    In order to get the correct package, the user should remove the
!    comments before the calls of the packages
!
!    The following specific comments are used
!
!    chgl   HP-GL plotter version
!    cpos   PostScript output
!    cown   Local version not supported by SEPRAN
! **********************************************************************
!
!                       DATA STATEMENTS
!
!
      data isopen /.false./
! ======================================================================
!
!     --- Initialiations

      iref = abs(irefpl)
      ireftk = iref+3
      if (ichois.eq.1) then

!     --- ichois = 1:  start

         xorig=0.
         yorig=0.
         sx=x
         sy=y

      else if (ichois.eq.4) then

!     --- ichois = 4:  change origin

         xorig=xorig+x
         yorig=yorig+y

      else

!     --- ichois <> 1 or 4:  Compute position

         sx=x+xorig
         sy=y+yorig

      end if

!     --- first reset all colors, init tables etc

      if (ichois.eq.1) call plcoli

!     --- some implementations need explicit setting of the default color

      if ((posdev(1:2).ne.'ps' .and. posdev(1:2).ne.'pl') .and. &
           dcolor.gt.0) call plcol0 (dcolor)

!     --- delete existing files if opening for first call

      if ( ichois.eq.1 .and. &
           (parallel .and. inodenr==1 .or. .not. parallel) ) then
         if ( iplpic.eq.0 .and. ( posdev(1:2).eq.'fa' .or. &
              posdev(1:2).eq.'fb' .or. iplotf.eq.1 ) ) then

!        --- Write file ASCII or binary
!            All files "namepl".*** are deleted

            leng = index ( namepl, ' ' )-1

10          iplpic=iplpic+1

            if ( maxplots.eq.0 ) then
               write (filenm,9000) namepl(1:leng), iplpic
9000           format (a,'.',i4.4)
            else
               write (filenm,9001) namepl(1:leng), iplpic
9001           format (a,'.',i5.5)
            end if

!           --- Check if file exists

            name = filenm
            if ( name(1:8).eq.'/tmp/SEP' ) then
               lastps = index ( name, ' ' )-1
               name = name(21:lastps)
            end if
            inquire ( file=name, exist=check )
            if ( check ) then

!              --- File exists, delete by opening and closing with delete

               open (unit=iref,file=name,status='old')
               close (unit=iref,status='delete')
               isopen = .false.
               go to 10

            end if
            iplpic=0

         else if ( iplpic.eq.0 .and. posdev(1:2).eq.'f1' ) then

!        --- Write one ASCII file
!            The file "namepl" is deleted

            leng = index ( namepl, ' ' )-1
            filenm = namepl(1:leng)
            name = filenm
            if ( name(1:8).eq.'/tmp/SEP' ) then
               lastps = index ( name, ' ' )-1
               name = name(21:lastps)
            end if

!           --- Check if file exists

            inquire ( file=name, exist=check )
            if ( check ) then

!           --- File exists, delete by opening and closing with delete

               open (unit=iref,file=name,status='old')
               close (unit=iref,status='delete')
               isopen = .false.

            end if
            iplpic = 0
            iplopn = 1

         end if

      end if

!     --- Binary (unformatted): open file-set and write all plotdata

      if ( iplotf.eq.1 .and. irefpl.lt.0 .or. posdev(1:2).eq.'fb') then

         if ( ichois.eq.1 ) then

!        --- ichois = 1, open plot

            iplpic=iplpic+1

            leng = index ( namepl, ' ' )-1
            if ( maxplots.eq.0 ) then
               write (filenm,9000) namepl(1:leng), iplpic
            else
               write (filenm,9001) namepl(1:leng), iplpic
            end if

            ! PvK overwrite plote name
            filenm=ourplotname

            open (unit=iref,file=filenm,form='unformatted',err=5000)
            isopen = .true.
            iplopn = 1

         end if

         if ( ichois.ne.4 ) then

!        --- ichois <> 4, write arguments (1 and 5 too!)

            if ( isopen ) then

               write (iref) ichois,sx,sy,0.,0.,dcolor

               if (ichois.eq.5) then

!              --- ichois = 5, close plot

                  close (unit=iref)
                  isopen = .false.
                  if ( x.le.-1d0 .and. x.ge.-1d0 ) iplpic = 0
                  iplopn = 0

               end if

            else if ( ichois.eq.5 ) then
               call eropen ( 'plafp6' )
               call errwar ( 294, 0, 0, 0 )
               call erclos ( 'plafp6' )
               iplopn = 0
            else
               call eropen ( 'plafp6' )
               call errint ( ichois, 1 )
               call errsub ( 2285, 1, 0, 0 )
               call erclos ( 'plafp6' )
               iplopn = 0
            end if

         end if
         if ( iplotf.eq.0 ) return

!        --- Ascii (formatted): open plotfile-set and write all plotdata

      else if ( iplotf.eq.1 .or. posdev(1:2).eq.'fa' ) then

         if (ichois.eq.1) then
            iplpic=iplpic+1
            leng = index ( namepl, ' ' )-1
            if ( maxplots.eq.0 ) then
               write (filenm,9000) namepl(1:leng), iplpic
            else
               write (filenm,9001) namepl(1:leng), iplpic
            end if
            ! PvK overwrite plotname
            filenm=ourplotname
            open (unit=iref,file=filenm,form='formatted',err=5000)
            isopen = .true.
            iplopn = 1
         end if

         if (ichois.ne.4) then

!           --- ichois <> 4, write arguments (1 and 5 too!)

            if ( isopen ) then

               write (iref,9030) ichois,sx,sy,dcolor
9030           format (i1,2(1x,f8.3),' 0 0 ',i3)

               if (ichois.eq.5) then
                  close (unit=iref)
                  isopen = .false.
                  if ( x.le.-1d0 .and. x.ge.-1d0 ) iplpic = 0
                  iplopn = 0
               end if

            else if ( ichois.eq.5 ) then
               call eropen ( 'plafp6' )
               call errwar ( 294, 0, 0, 0 )
               call erclos ( 'plafp6' )
               iplopn = 0
            else
               call eropen ( 'plafp6' )
               call errint ( ichois, 1 )
               call errsub ( 2285, 1, 0, 0 )
               call erclos ( 'plafp6' )
               iplopn = 0
            end if

         end if
         if ( iplotf.eq.0 ) return
      end if

!     --- One Ascii file only (formatted): open plotfile-set if necessary
!         and write all plotdata

      if (posdev(1:2).eq.'f1') then
         if (ichois.eq.1 .and. iplpic.eq.0) then
            iplpic=iplpic+1
            write (filenm,9000) namepl
            open (unit=iref,file=filenm,form='formatted',err=5000)
            isopen = .true.
         end if
         iplopn = 1

         if (ichois.ne.4) then

!        --- ichois <> 4, write arguments (1 and 5 too!)

            if ( isopen ) then
               write (iref,9030) ichois,sx,sy,dcolor
            else if ( ichois.eq.5 ) then
               call eropen ( 'plafp6' )
               call errwar ( 294, 0, 0, 0 )
               call erclos ( 'plafp6' )
               iplopn = 0
            else
               call eropen ( 'plafp6' )
               call errint ( ichois, 1 )
               call errsub ( 2285, 1, 0, 0 )
               call erclos ( 'plafp6' )
               iplopn = 0
            end if

         end if

         go to 1000

!         PostScript plotter output to file:

      else if ( posdev(1:2).eq.'ps' .or. posdev(1:2).eq.'pl' &
           .or. posdev(1:2).eq.'pc' .or. posdev(1:2).eq.'pw' ) then
!
         call plpos6 ( ichois, sx, sy, 0, 0, 560.368, 806.144 )
         return

!     -----------------------------------------------------------------
!         Local version not supported by SEPRAN

      else if ( posdev(1:2).eq.'ow' ) then
!
         call plown6 ( ichois, sx, sy )
         return

!     -----------------------------------------------------------------
!         HP-GL plotter output in file

      else if (posdev(1:1).eq.'p') then
         call plhgl6 ( ichois, sx, sy )
         return


!     -----------------------------------------------------------------
!         Unknown choice, unsupported or commented out of compiled version

      else
         go to 5200
      end if

6000  return

!     -----------------------------------------------------------------
!        Some error messages

5000  call eropen ( 'plafp6' )
      call errsub ( 1996, 0, 0, 0 )
      call erclos ( 'plafp6' )
      stop
!
5100  call eropen ( 'plafp6' )
      call errint ( ichois, 1 )
      call errsub ( 1997, 1, 0, 0 )
      call erclos ( 'plafp6' )
      stop
!
5200  call eropen ( 'plafp6' )
      call errchr ( posdev, 1 )
      call errsub ( 1998, 0, 0, 1 )
      call erclos ( 'plafp6' )
      stop
!
1000  end
