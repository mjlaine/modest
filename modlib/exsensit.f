
      subroutine exsensit(ar,nr,ai,ni,
     &                  xest,nest,lest,
     &                  xdisc,valdisc,ndiscpoint)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4 nest
      integer*4 ndiscpoint,moststep,ndir
      real*8    xest(nest)
      integer*4 lest(4,nest)
      real*8    xdisc(nest,ndiscpoint)
      real*8    valdisc(ndiscpoint)


*     LOCAL VARIABLES

      integer*4 i,j,k,i1,i2,isen,ii,jj,kk,nresultfile,ndir1,count
      real*8    lsq,expdes,opt,r2
      external  lsq,expdes,opt,r2

      include 'common3.inc'

      data   count  /0/

*      read arbitrary  variable points directly to 'xdisc'
*      from a user given matrix in file 'discfile'

       if (discfile(1:5) .ne. 'dummy') then
         count = count + 1
         if (count.eq.1) then
          if (ndiscpoint.gt.0) then
            open(unit=1,file = discfile, status = 'old')
            do j = 1,ndiscpoint
               read(1,*) (xdisc(i,j), i=1,nest)
            enddo
            close(1)
          else
            write(*,*) 'ndiscpoint not given'
            stop
          endif
         endif

         do jj = 1,ndiscpoint
           if (debug.gt.10) write(*,*) 'calling lsq in exsensit '
           valdisc(jj) = lsq(ar,nr,ai,ni,xdisc(1,jj),nest)
           call Xparser(ar,nr,ai,ni,'par',xest,nest,lest)
         enddo

       endif

      return
      end


