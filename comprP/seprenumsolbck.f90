      subroutine seprenumsolbck ( solin, ichoice, usolglob )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 17-03-2016
!
!   copyright (c) 2016-2016  "Ingenieursbureau SEPRA"
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
!     Renumber solution vector back to original form
!     Real case
! **********************************************************************
!
!                       KEYWORDS
!
!     real
!     reordering
! **********************************************************************
!
!                       MODULES USED
!
      use sepmodulecpack
      use sepmodulematstr
      use sepmodulekprob
      use sepmodulesolve
      implicit none
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer, intent(in) :: ichoice
      double precision :: solin(nusol)
      double precision, intent(in), optional :: usolglob(*)
      integer :: i=0,ip=0

!     ichoice        i    If > 0 it concerns the solution vector
!     solin         i/o   Solution vector to be transformed
!     usolglob       i    Structure of type vector containing the necessary
!                         information of the global solution vector
!                         i.e. the part that is not positioned in this cpu
!                         but is necessary for the solution at this cpu
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     SEPLOCTRNRL    Compute (back-)transformation in case of local
!                    transformations ( real case )
!     SEPSOLCOPSOL   Send parts of local solution array to neighbor blocks
!                    Receive these solution arrays and store in usolglob
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
      call eropen ( 'seprenumsolbck' )
      debug = .false. .and. ioutp>=0
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from seprenumsolbck'
         write(irefwr,1) 'ichoice, npltra, typecoupl', &
                          ichoice, npltra, typecoupl
  1      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      if ( ichoice==1 .and. parallel .and. typecoupl==1 ) then

!     --- Extend usol to xglob

         call sepsolcopsol ( solin, usolglob )

         if ( ninterface>0 ) &
            solin(nrusol+1:nrusol+ninterface) = &
                usolglob(kprobhglob(nrusol+1:nrusol+ninterface))

      end if  ! ( ichoice==1 .and. parallel )

      npboun = kprobloc(34)
      if ( npboun>0 ) then

!     --- npboun>0; compute effect of periodical boundary conditions

         kprobs => kprobpart%kprobs
         rprobs => kprobpart%rprobs
         solin(nrusol+1:nrusol+npboun) = rprobs(:,2)*solin(kprobs)+rprobs(:,1)

      end if  ! ( npboun>0 )

!     --- Transform essential boundary conditions in case of local transform

      if ( npltra>0 .and. itrans==0 ) then
         indprh = 0
         call seploctrnrl ( 2, solin )
         indprh = kprobloc(50)
      end if  ! ( npltra>0 .and. itrans==0 )

!     ! write out these do loops to avoid segfaults
!     if ( indprh>0 ) solin(kprobhinv) = solin  ! reorder solin back

!     if ( ichoice==1 .and. indprh>0 ) usolbuf(1:nusol) = solin
      if ( indprh>0 .and. ichoice==1 ) then
         do i=1,nusol
            ip=kprobpart%kprobhinv(i)
            usolbuf(ip)=solin(i)
         enddo
      else if ( indprh>0 ) then
         solin(kprobhinv) = solin  ! reorder solin back
         if ( ichoice==1) usolbuf(1:nusol) = solin
      endif


1000  if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End seprenumsolbck'

      end if  ! ( debug )
      call erclos ( 'seprenumsolbck' )

      end subroutine seprenumsolbck
