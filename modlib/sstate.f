      subroutine SSTATE(ar,nr,ai,ni,
     &  xaux,naux,
     &  states,nstatea,mostobs,nsaveset)

*     This is an auxiliary function CALLED BY the selected optimizer (e.g.
*     simflex) via LSQ *only*. Reason: choise of task ('laux') in LSQ.

*     SUBROUTINES NEEDED

*     - Cstates

*     PURPOSE

*     Returns the states in a point specified by 'naux,laux,xaux'

*     COMMON BLOCKS

      include 'common3.inc'

*     ARGUMENTS

      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4      naux,nstatea,mostobs,nsaveset
      real*8         xaux(naux)       ! variables to be optimized
      real*8         states(nstatea,mostobs,nsaveset)

*     LOCAL VARIABLES

      integer*4      i

*     THE ALGORITHM

*     Compute the states: estimate or simulate

      if (debug .gt. 10) write(*,*) ' in sstate'

      if (task(1:3) .eq. 'est' .or. task(1:3) .eq. 'exp' ) then
      if(debug .gt. 10) write(*,*) ' calling cstates, sstate/est'
         call cstates(ar,nr,ai,ni,xaux,naux,ai(plest),
     &                states,nstatea,mostobs,nsaveset)

      elseif (task(1:3) .eq. 'sen') then
       if (optcrit(1:1).eq.'G' .or. optcrit(1:1).eq.'g') then
         call cstates(ar,nr,ai,ni,xaux,naux,ai(plest),
     &                states,nstatea,mostobs,nsaveset)
       else
         call cstates(ar,nr,ai,ni,xaux,naux,ai(plsen),
     &                states,nstatea,mostobs,nsaveset)
       endif
      elseif (task(1:3) .eq. 'opt') then
         call cstates(ar,nr,ai,ni,xaux,naux,ai(plopt),
     &                states,nstatea,mostobs,nsaveset)

      else
         write(*,*) 'task not given in sstate'
         stop
      endif

      return
      end

