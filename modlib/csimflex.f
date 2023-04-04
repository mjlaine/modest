      subroutine csimflex(ar,nr,ai,ni)

      implicit none

*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

*     LOCAL PARAMETERS

      real*8      optim
      integer*4   simpsize
      integer*4   ierror,npen
      real*8      lsq,expdes,opt,R2
* removed by ML 15.8.2002
*     real*8      lsqlmdif
      external    lsq
*     external    lsqlmdif
      external    expdes
      external    opt                       
      external    R2

*    COMMON BLOCKS

      include 'common3.inc'

	

      ierror  = 1                  !non-monitor mode
      if(optmonit .gt. 0) ierror=-1

      if (task(1:3).eq.'est') then

         simpsize = dnest+1

*        Select the arguments to be optimized

         if(debug .gt. 10) write(*,*) 'calling Xparser/est in csimflex'
         call Xparser(ar,nr,ai,ni,
     &                'sel',ar(pxest),dnest,ai(plest))

         npen    = 1                  !n. of penalty param.(even for dummy)??

         if (objfun(1:3).eq.'lsq') then
           if (debug .gt. 10) write(*, *) 'calling Simflex, est/lsq'

           call Simflex(ar,nr,ai,ni,
     &                lsq,
     &                ar(pxest),dnest,
     &                sizes,itmaxs,abstols,reltols,
     &                ar(pbounds),ai(pbtype),
     &                optim,ar(psimplex),simpsize,
     &                ar(pwrksimf),ai(piwrsimf),ierror,optmonit)

         elseif (objfun(1:2).eq.'r2') then
           if (debug .gt. 10) write(*, *) 'calling Simflex, est/R2'

           call Simflex(ar,nr,ai,ni,
     &                R2,
     &                ar(pxest),dnest,
     &                sizes,itmaxs,abstols,reltols,
     &                ar(pbounds),ai(pbtype),
     &                optim,ar(psimplex),simpsize,
     &                ar(pwrksimf),ai(piwrsimf),ierror,optmonit)

         elseif (objfun(1:3).eq.'opt') then
           if (debug .gt. 10) write(*, *) 'calling Simflex, est/opt'

           call Simflex(ar,nr,ai,ni,
     &                opt,
     &                ar(pxest),dnest,
     &                sizes,itmaxs,abstols,reltols,
     &                ar(pbounds),ai(pbtype),
     &                optim,ar(psimplex),simpsize,
     &                ar(pwrksimf),ai(piwrsimf),ierror,optmonit)

         else
           write(*,*) 'no such objfun'
           stop
         endif ! task:est

      elseif (task(1:3).eq.'exp') then

         simpsize = dnexp+1

         if(debug .gt. 10) write(*,*) 'calling Xparser/exp in csimflex'

         call xparser(ar,nr,ai,ni,
     &                'sel',ar(pxexp),dnexp,ai(plexp))

c         npen    = 1                  !n. of penalty param.(even for dummy)??
          ierror  = 1                  !non-monitor mode

         if(optmonit .gt. 0) ierror=-1

c         if (objfun(1:3).eq.'exp') then
           if (debug .gt. 10) write(*, *) 'calling Simflex, exp'

           call Simflex(ar,nr,ai,ni,
     &                expdes,
     &                ar(pxexp),dnexp,
     &                sizes,itmaxs,abstols,reltols,
     &                ar(pbounds),ai(pbtype),
     &                optim,ar(psimplex),simpsize,
     &                ar(pwrksimf),ai(piwrsimf),ierror,optmonit)

c         elseif (objfun(1:3).eq.'opt') then
c           if (debug .gt. 10) write(*, *) 'calling Simflex, exp/opt'
c
c           call Simflex(ar,nr,ai,ni,
c     &                opt,
c     &                ar(pxest),dnest,
c     &                sizes,itmaxs,abstols,reltols,
c     &                ar(pbounds),ai(pbtype),
c     &                optim,ar(psimplex),simpsize,
c     &                ar(pwrksimf),ai(piwrsimf),ierror,optmonit)
c         endif ! task:exp

      elseif (task(1:3).eq.'opt') then

         simpsize = dnopt+1

         if(debug .gt. 10) write(*,*) 'calling Xparser (opt)'

         call xparser(ar,nr,ai,ni,
     &                'sel',ar(pxopt),dnopt,ai(plopt))

         npen    = 1                  !n. of penalty param.(even for dummy)??
         ierror  = 1                  !non-monitor mode

         if(optmonit .gt. 0) ierror=-1

         if (debug .gt. 10) write(*, *) 'calling Simflex, opt'

         call Simflex(ar,nr,ai,ni,
     &                opt,
     &                ar(pxopt),dnopt,
     &                sizes,itmaxs,abstols,reltols,
     &                ar(pbounds),ai(pbtype),
     &                optim,ar(psimplex),simpsize,
     &                ar(pwrksimf),ai(piwrsimf),ierror,optmonit)


      endif

      return
      end


