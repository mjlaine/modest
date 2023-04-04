      subroutine Clmdiff(ar,nr,ai,ni)

      implicit none

*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

*     LOCAL PARAMETERS

      real*8      optim
      real*8      lsq
      external    lsq
* removed by ML 15.8.2002
*     real*8      lsqlmdif
*     external    lsqlmdif

*    COMMON BLOCKS

      include 'common3.inc'

      if (task(1:3).eq.'est' .and. objfun(1:3).eq.'lsq') then

*        Select the arguments to be optimized

         if(debug .gt. 10) write(*,*) 'calling Xparser/est in CLevMar'
         call Xparser(ar,nr,ai,ni,
     &                'sel',ar(pxest),dnest,ai(plest))
         
         if (debug .gt. 10) write(*, *) 'calling LevMar, est/lsq'

         call Lmdiff(ar,nr,ai,ni,
     &        lsq,
     &        ar(pxest),dnest,ai(plest),
     &        ai(pipivlm),ar(pzlm),
     &        ar(pbounds),ai(pbtype),
     &        ai(pnobs),dnsets,
     &        ar(pstates),ar(pstates0),dnstatea,dmostobs,dnsavese,
     &        ar(pxdata),dnx,
     &        ar(pydata),ar(pweight),
     &        ar(pyest),ar(pxjac),ar(pjtj),ar(pjtj1),ar(pxaux),
     &        dmostyda,ai(pnydata),
     &        ar(pgpar),max(dngpar,1),
     &        ar(plpar),max(dnlpar,1))

      else
        write(*,*) 'LevMar minimizer only for '
        write(*,*) 'least squares fitting, objfun = lsq '
        stop
      endif

      return
      end


