      subroutine to0092 ( isubr, iread, maxmesh1 )
! ======================================================================
!
!        programmer    Guus Segal
!        version  6.1  date 11-12-2017 New storage of intmat
!        version  6.0  date 30-01-2016 New parameter list
!        version  5.32 date 17-12-2015 New call to incommat
!
!   copyright (c) 1993-2017  "Ingenieursbureau SEPRA"
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
!    Undersubroutine of SEPSTR and SEPSTL
!    The actual tasks of SEPSTR/SEPSTL are carried out
!    The mesh is supposed to be made by program SEPMESH and stored in the
!    file meshoutput
!    The input for subroutine PROBDF is stored in array iinput and written
!    to the file sepcomp.inf
!    The files meshoutput, sepcomp.inf and sepcomp.out are opened, &
!    the file meshoutput is closed.
!    The machine dependent reference numbers iref10, iref73 and iref74 get
!    a value
! **********************************************************************
!
!                       KEYWORDS
!
!     read
!     start
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulekmesh
      use sepmodulekprob
      use sepmodulemain
      use sepmoduleintmat
      use sepmodulemat
      use sepmoduleinput
      use sepmodulecpack
      use sepmodulemach
      use mparallel
      implicit none

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer maxmesh1, isubr, iread

!     iread          i    Indicates if all SEPRAN input until end
!                         of file or END_OF_SEPRAN_INPUT is read (1) or that
!                         only the input as described for SEPSTR is read (0)
!     isubr          i    Indicates the type of calling subroutine
!                         Possible values
!                         1:  SEPSTR, nbuf does not have a value
!                         2:  SEPSTL, nbuf has a value
!                         3:  SEPSTN, nbuf may have a value
!                             All SEPRAN input is read in this subroutine
!                             provided iread = 1
!                         4:  See 3, now called by sepcomp or sepfree
!     maxmesh1       i    Maximum number of meshes allowed
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer lniincommat
      parameter ( lniincommat=25 )
      integer iincommt(lniincommat,10), &
              iref, jstart, iseqnr, idummy, ihelp(10), &
              ishift, lastps, lastps1, manag(lenmanag), nextrecord, &
              iseqprob, irefsav, lenhome, maxcall, ipos
      logical check, callcommat
      character (len=100) filename, filename1, name
      type (rmatrix) :: A

!     A              Structure of type matrix containing the necessary
!                    information of the matrix
!     callcommat     If true subroutine commat must be called
!     check          Indication if a file has been opened (true) or not (false)
!     filename       Name of file
!     filename1      Name of file
!     idummy         Dummy parameter
!     ihelp          Work array to store input from the input file
!     iincommt       Input array for subroutine incommat.
!                    Defines the type of matrix to be used. (See User's Manual)
!                    iincommt is an array of length 10, which means that at most
!                    10 problems may be solved at one time
!     ipos           Position of input arrays in sepmoduleinput
!     iref           Local reference number
!     irefsav        Help parameter to save the value of iref
!     iseqnr         Dummy parameter
!     iseqprob       Sequence number of problem input
!     ishift         Shift in mesh sequence number
!     jstart         Absolute value of first parameter in the call of START
!     lastps         Last non-blank position in filename
!     lastps1        Last non-blank position in filename1
!     lenhome        Number of characters in the name of the SEPRAN home
!                    directory
!     lniincommat    Length of first index of iincommat
!     manag          Manager array that sets the parameters in subroutine
!                    SEPSTL
!                    See subroutine readsp
!     maxcall         Maximum number of coefficient inputs read
!     name           Name of file
!     nextrecord     Defines last record read by readsp
!                    Possible values:
!                    0: none
!                    1: PROBLEM
!                    2: POSTPROCESSING
!                    3: MESH
!                    4: READ MESH
!                    5: EXTRA_MESH_INPUT
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
      integer inisetnewref

!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERRCHR         Put character in error message
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     INCOPY         Copy one integer array into another
!     INICHINPSOL    Adapt input array with respect to variables
!     INIFIL         Open files
!     INIMPI         Start MPI
!     INIMPIOPEN     Join parallel application
!     INIRMF         Remove file if it exists
!     INIRMREF       Remove entry iref from the set of reference numbers used
!     INISETNEWREF   Create new reference number
!     INITCB         Initializes SEPRAN machine-independent common blocks
!     INITMD         Set machine-dependent quantities
!     MSH071         Read information of the mesh from the file indicated by
!                    iref
!     MSHEXTRA       Read extra input for mesh and apply the actions required
!                    by this extra input
!     MSHMESHQUAN    Fill contents of kmesh2d_quant
!     MSHRDPARMESH   Read information parallel information about global mesh
!                    from file meshoutput_par.000
!     PRININ         print 1d integer array
!     PROBDFEXT      Reads and stores problem definition
!     READ00         Reads input for subroutine COMMAT from standard input
!                    file.
!     READ74SUB      Read names of vectors from the sepcomp.out file
!     READAL         Read all SEPRAN input for the computational program except
!                    problem definition and preceding information
!     READPROBIINP   Read information concerning the input for probdfbf from
!                    file sepcomp.out and call subroutine probdfbf
!     READSP         Reads information for subroutine SEPSTL
!     SEPCOMMAT      Actual body of COMMAT and MATSTRUC
!     SEPFEMESH      creates a linear element mesh based on the nodes
!                    of a spectral element mesh
!     SEPFILLKMESH   Fill parameters in module sepmodulekmesh
!     SEPGETINPFINF  Extract infor entries of integer and real input arrays from
!                    ISEPIN
!     SEPLOGO        Print copyright statement of SEPRAN
!     SEPMUMPSSLAVE  Perform a repeated call to prmumps on a slave
!     SEPRELKMESH    Make the relation between kmref and the pointers
!                    in sepmodulekmeshp
!     SEPSTART       Standard SEPRAN starting subroutine
!     WRITPROBIINP   Write arrays iinput and rinput corresponding to kprob2d to
!                    file sepcomp.out
! **********************************************************************
!
!                       I/O
!
!    Input is read from the standard input file
!    See the various subroutines that perform this action
!    Furthermore the file meshoutput is read
!    File sepcomp.inf is created
! **********************************************************************
!
!                       ERROR MESSAGES
!
!     719   Too much problems are to be created at one time
!    1669   file does not exist
!    2860   The sepcomp.in and sepcomp.out file have the same name
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
      debug = .false.
  1   format ( a, 1x, (10i6) )

!     --- Start the program

      call initcb
      if ( debug ) write(6,*) 'after initcb'

!     --- Check if parallel computations must be carried out

      call inimpiopen
      inquire ( file = 'sepran_par.input', exist = parallel )
      inquire ( file = 'sepran_mumps.input', exist = usemumps )

      if ( parallel ) then

!     --- File sepran_par.input has been found, i.e. parallel computations
!         initialize mpi

         call inimpi
         cputiming = .false.

      else if ( usemumps ) then

!     --- File sepran_mumps.input has been found, i.e. mumps solver
!         initialize mpi

         call inimpi
         firstmumps = .false.
         cputiming = .false.

      else

!     --- Serial computions

         cputiming = .true.

      end if  ! ( parallel )

      call initmd
      if ( debug ) write(irefwr,*) 'after initmdbf'

!     --- Open read and write files

      call inifil(5)
      if ( parallel .or. usemumps ) then

!     --- Parallel, we rename irefwr to 26 and couple a name to
!         the output file

         irefwr = 26
         name(1:11) = 'sepran_out.'
         write ( name(12:14),100 ) inodenr
100      format(i3.3)
         open ( unit = irefwr, file = name(1:14) )
         write(irefwr,*) 'node number is ',inodenr

      else

!     --- Serial computation, standard

         call inifil(6)

      end if
      irefwrsave = irefwr

!     --- Open standard elements file and error message file

      call inifil(1)
      call inifil(4)
      if ( debug ) write(irefwr,*) 'after inifil'

      call seplogo      ! print copyright statement

      if ( parallel ) then

!     --- File sepran_par.input has been found, i.e. parallel computations
!         First set input file

         open ( unit=irefre, file = 'sepran_par.input' )

      else if ( usemumps ) then

!     --- File sepran_mumps.input has been found, i.e. computations
!         using mumps
!         First set input file
!         In case we are not at the host, prmumpsslave is called
!         This carries out an infinite loop in which only some
!         preparations for each mumps call are made
!         Also it takes care of finishing the process

         open ( unit=irefre, file = 'sepran_mumps.input' )
         if ( inodenr>1 ) call sepmumpsslave (A)

      end if  ! ( parallel )

!     --- Read array manag (defaults)

      manag = 0
      manag(5) = -2
      if ( isubr>=4 .and. isubr/=11 ) then
         manag(5) = manag(5)-1000
         mainsubr = 2
      end if

      call readsp ( manag, lenmanag, nextrecord )
      if ( debug ) then
         write(irefwr,*) 'after readsp'
         write(irefwr,*) 'parallel', parallel
         call prinin ( manag, lenmanag, 'manag' )
      end if  ! ( debug )

!     --- Check if the file $SPHOME/bin/update/numwis exists
!         Set numwis

      lenhome = index(sphome,' ')-1
      filename = sphome(1:lenhome)//'/bin/update/numwis'
      inquire ( file = filename, exist = check )
      numwis = check

      if ( numwis .or. mpi_partrac .or. mpi_parallel ) manag(7) = 0  ! do not write sepcomp.out

      if ( parallel ) then

!     --- create name of general parallel mesh input file

         lastps = index ( name10, ' ' )-1
         filename = name10
         filename(lastps+1:lastps+8) = '_par.000'
         inquire ( file = filename(1:lastps+8), exist = check )
         if ( .not. check ) then

!        --- The file meshoutput_par.000 does not exist

            call errchr ( filename(1:lastps+8), 1 )
            call errsub ( 1669, 0, 0, 1 )
            go to 1000

         end if

!        --- Set name of mesh input file

         write(filename(lastps+1:lastps+8),110) inodenr
110      format('_par.',i3.3)
         name10 = filename(1:lastps+8)

         lastps1 = index ( name74, ' ' )-5
         filename1 = name74
         write(filename1(lastps1+1:lastps1+8),110) inodenr
         name74 = filename1(1:lastps1+8)

      end if

!     --- Standard SEPRAN start subroutine

      jstart = abs(manag(1))
      if ( mod(jstart/100,10)==0 ) jstart = jstart + 100
      if ( manag(1)<0 ) then
         manag(1) = -jstart
      else
         manag(1) = jstart
      end if

      if ( debug ) call prinin ( manag, lenmanag, 'manag' )

      call sepstart ( manag(1), manag(2), manag(3), manag(4) )
      if ( debug ) write(irefwr,*) 'after sepstart'
      if ( isubr==1 ) then
         call eropen( 'sepstr' )
      else if ( isubr==2 ) then
         call eropen( 'sepstl' )
      else
         call eropen( 'sepstn' )
      end if
      call eropen ( 'to0092' )

!     --- Open files meshoutput, sepcomp.out and read mesh

      call inifil ( -10 )
      if ( debug ) write(irefwr,*) 'after inifil ( -10 )'
      if  ( manag(11)==1 ) then

!     --- manag(11)>0, read from sepcomp.xxx

         if ( iref74<0 ) iref73 = -abs(iref73)  ! sepcomp.in has same
                                                ! structure as sepcomp.out
         call inifil ( -73 )

      end if  ! ( manag(11)==1 )

      if ( manag(7)==1 ) then

!     --- manag(7) = 1, fill files 74 (sepcomp.out)

         if ( manag(11)==0 ) then

!        --- manag(11)=0, standard

            call inirmf ( 73 )  ! remove sepcomp.inf
            call inirmf ( 74 )  ! remove sepcomp.out

         end if  ! ( manag(11)==0 )

         if ( manag(7)==1 ) then

!        --- manag(7) = 1, write to sepcomp.out

            call inifil ( 74 )
            if ( manag(11)>0 ) then

!           --- manag(11)>0, in and output files must have
!               different names

               if ( name74==name73 ) then

!              --- sepcomp.in and sepcomp.out have the same name

                  call errchr ( name73, 1 )
                  call errsub ( 2860, 0, 0, 1 )
                  go to 1000

               end if  ! ( namef74==namef73 )

            end if  ! ( manag(11)>0 )

         end if  ! ( manag(7)==1 )

         iref = inisetnewref()
         inquire ( file='sepcomp.freq', exist=check )
         if ( check ) then

            open ( unit = iref, file = 'sepcomp.freq', &
                   status = 'unknown' )
            close ( unit = iref, status = 'delete' )
            call inirmref ( iref )

         end if

      end if  ! ( manag(7)==1 )

      kmesh => kmesh2d(1)
      kprob => kprob2d(1)

      if ( debug ) write(irefwr,*) 'before msh071'
      call msh071 ( manag(5), iref10, ihelp )
      if ( debug ) write(irefwr,*) 'after msh071'
      iref = abs(iref10)
      close ( unit = iref )

      if ( parallel ) then

!     --- parallel is true, hence parallel computing
!         Read file meshoutput_par.000 and some information from meshoutput

         call mshrdparmesh

      end if
      if ( debug ) write(irefwr,*) 'after mshrdparmesh'

!     --- Check if spectral mesh has been found
!         If so call femesh

      if ( isubr==4 .and. kelmu/=0 ) then

!     --- Spectral mesh: fill fem mesh in kmesh2d(.,2)

         iref = inisetnewref()
         call sepfemesh ( iref, 1 )
         call inirmref ( iref )
         close ( iref )
         if ( debug ) write(irefwr,*) 'after femeshbf'

      end if

      if ( nextrecord==5 ) then

!     --- Read extra mesh input and perform corresponding actions

         call mshextra
         if ( debug ) write(irefwr,*) 'after mshextra'

      end if  ! ( nextrecord==5 )

!     --- Compute some mesh dependent quantities and store in kmesh2d_quant

      kmesh => kmesh2d(1)

      call mshmeshquan
      if ( debug ) write(irefwr,*) 'after mshmeshquan'

      if ( isubr==4 .and. kelmu/=0 ) then
         kmesh => kmesh2d(2)
         call seprelkmesh ( kmesh )
      end if  ! ( isubr==4 .and. kelmu/=0 )

!     --- Read problem definition

      if ( debug ) write(irefwr,1) 'manag(11)', manag(11)

      if ( manag(11)<=0 ) then

!     --- manag(11)=0, standard case

         if ( debug ) write(irefwr,1) 'ichprb', manag(6)
         maxcontc = 0
         call probdfext ( manag(6), iseqprob )
         if ( debug ) write(irefwr,*) 'after probdfext'

      else if ( manag(11)>0 ) then

!     --- manag(11)=1, read problem definition

         irefsav = iref74    ! save iref74
         iref74 = iref73     ! replace by iref73
         call readprobiinp

!        --- Read names of vectors

         call read74sub
         iref74 = irefsav    ! reset iref74
         iseqprob = 1

      end if  ! ( manag(11)==0 )
      if ( ierror/=0 ) go to 1000

      if ( debug ) then
         write(irefwr,*) 'after probdfext'
         write(irefwr,200) 'iseqprob', iseqprob
200      format ( a, 1x, 100(i6,1x) )
      end if  ! ( debug )

      if ( manag(11)>=0 .and. manag(7)==1 ) then

!     --- Write array iinput to file sepcomp.out

         call writprobiinp
         if ( debug ) write(irefwr,*) 'after writprobiinp'

      end if
      iseqkp = kprob2d(1)
      nprob = kp(iseqkp)%nprob
      if ( nprob>10 ) then

!     --- nprob>10 not yet implemented

         call errint ( nprob, 1 )
         call errint ( 10, 2 )
         call errsub ( 719, 2, 0, 0 )
         go to 1000

      end if

      if ( manag(11)<0 .or. (manag(11)==0 .and. iread==1) ) then

!     --- manag(11)<0, probdf is skipped as well as commat

         callcommat = .false.  ! commat must not be called

      else

!     --- manag(11)>=0, standard case

         callcommat = .true.  ! default commat must be called

      end if  ! ( manag(11)<0 )
      if ( debug ) then
         write(irefwr,*) 'callcommat/1', callcommat
         write(irefwr,1) 'iread', iread
      end if  ! ( debug )

      kmesh => kmesh2d(1)
      call sepfillkmesh(kmesh)

      if ( iread==1 ) then

!     --- isubr = 1, Read rest of the SEPRAN input

         call readal ( iseqprob, manag )
         if ( debug ) write(irefwr,*) 'after readal'
         if ( ierror/=0 ) go to 1000

!        --- Find first input of sepcommat and store in iincommt

         call sepgetinpfinf ( 3, 1, ipos, maxcall )
         if ( ierror/=0 ) go to 1000

         if ( debug ) write(irefwr,1) 'ipos, maxcall', ipos, maxcall

         if ( in(ipos)%used ) then

!        --- iinput(1)>0, hence input for commat is defined
!            Copy into iincommt

            call incopy ( in(ipos)%iinput, iincommt, nprob*lniincommat )
            callcommat = .true.  ! commat must not be called

         else

!        --- iinput(1)=0, hence input for commat is not defined

            callcommat = .false.  ! commat must not be called

         end if  ! ( iinput(1)>0 )

      else

!     --- Read iincommt and call commat

         iseqnr = 1
         call read00 ( iincommt, nprob, iseqnr, 0, lniincommat )
         if ( debug ) write(irefwr,*) 'after read00'

      end if
      if ( ierror/=0 ) go to 1000

      if ( debug ) write(irefwr,*) 'callcommat/2', callcommat
      if ( callcommat ) then

!     --- callcommat true, call subroutine commat
!         Adapt variables

         call inichinpsol ( iincommt, lniincommat*nprob )
         if ( debug ) write(irefwr,*) 'after inichinpsol'
         if ( ierror/=0 ) go to 1000
         do iprob = 1, nprob
            ishift = iincommt(9,iprob)
            intmat => intmat1d(iprob)
            kmesh => kmesh2d(1+ishift)
            kprob => kprob2d(1)
            call sepcommat ( iincommt(1,iprob) )
            intmt(intmat)%spectral = ishift==1
         end do
         if ( debug ) write(irefwr,*) 'after incommat'

      end if  ! ( callcommat )

1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End to0092'

      end if
      call erclos ( 'to0092' )

      if ( isubr==1 ) then
         call erclos( 'sepstr' )
      else if ( isubr==2 ) then
         call erclos( 'sepstl' )
      else
         call erclos( 'sepstn' )
      end if

      end subroutine to0092
