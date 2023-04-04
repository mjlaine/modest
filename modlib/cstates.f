      subroutine cstates(ar,nr,ai,ni,xaux,naux,laux,
     &                   states,nstatea,mostobs,nsaveset)

*     This is a MODEST-subroutine CALLED BY the selected objective function
*     (eg. LSQ) or by JACOB

*     SUBROUTINES NEEDED
*
*     - Xparser
*     - Model subroutines

*     PURPOSE
*
*     Calculates the states with given xaux, for a given data set (& observn)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      Integer*4      naux
      Integer*4      laux(4,naux)
      real*8         xaux(naux)    ! variables to be optimized/simulated
      Integer*4      nstatea,mostobs,nsaveset
      real*8         states(nstatea,mostobs,nsaveset)

*     LOCAL VARIABLES

      Integer*4   i,j,k,k1

      include 'common3.inc'

      if(debug .gt. 11) write(6,*) '1. line of cstates'
      if(debug .gt. 10) write(6,*) 'calling xparser '

      call Xparser(ar,nr,ai,ni,'par',xaux,naux,laux)

*     Select the right subroutine to calculate the states

      if (model(1:3) .eq. 'alg') then
          if(debug .gt. 10) write(6,*) 'calling salg'
          call Salg(ar,nr,ai,ni,
     &              states,nstatea,mostobs,nsaveset,
     &              ar(pxdata), dnx,dnsets,
     &              ai(pnobs),
     &              ar(pgpar),max(dngpar,1),
     &              ar(plpar),max(dnlpar,1))

      elseif (model(1:3) .eq. 'imp') then
          if(debug .gt. 10) write(6,*) 'calling simp'
          call Simp(ar,nr,ai,ni,
     &      ar(1), ar(pstates0), states, nstatea,
     &      ar(pyest),dmostyda,mostobs,ai(pnydata),
     &      nsaveset, ar(pxdata), dnx, dnsets, ai(pnobs),
     &      ar(pgpar), max(dngpar,1), ar(plpar), max(dnlpar,1))

      elseif (model(1:3) .eq. 'ode') then
           if(debug .gt. 10) write(6,*) 'calling sode'

          call Sode(ar,nr,ai,ni,
     &      ar(1), ar(pstates0), states, nstatea,
     &      ar(pyest),dmostyda,mostobs,ai(pnydata),
     &      nsaveset, ar(pxdata), dnx, dnsets, ai(pnobs),
     &      ar(pgpar), max(dngpar,1), ar(plpar), max(dnlpar,1))

*      elseif (model(1:3) .eq. 'pde') then
*          call Spde(ar,nr,ai,ni)
*      elseif (model(1:3) .eq. 'dae') then
*          call Sdae(ar,nr,ai,ni)
      else
          write(*,*) 'No such module for states: ',model
      endif


      return
      end

