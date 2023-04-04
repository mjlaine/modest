

      subroutine xdirline(iaux,naux,ndir,dir,x0,bounds,nevastep,
     &                   xdir,moststep,coline,ncodir,tlaux,tuaux)

*     PARAMETERS

      implicit none
      integer*4 iaux,naux,ndir,moststep,ncodir
      real*8    dir(naux)
      real*8    x0(naux)
      real*8    bounds(2,naux)
      integer*4 nevastep(ndir)
      real*8    xdir(naux,moststep)
      real*8    coline(moststep)
      real*8    tlaux(naux)
      real*8    tuaux(naux)

*     LOCAL VARIABLES

      integer*4   i,j
      integer*4   idmin,idmax       ! blas-routines
      real*8      deltat            ! size of step
      real*8      tt,tl,tu           ! values of parameter tt corresponding
      real*8      aux1,aux2

****************************************************************
*     compute the points (inside "box") of line that goes
*     through point x0 in the direction 'dir'.
*
*     line is defined as        xdir = x0 + tt*dir
*
*     the coordinates of lines parallel to single variables
*     (iaux <= naux) saved in 'coline'
****************************************************************

*      first find out:  tl,tu,deltat,h

        do i=1,naux
           if (dir(i) .eq. 0) then
              tlaux(i) = -1.0d10
              tuaux(i) = 1.0d10
           else
              aux1 = (bounds(1,i)-x0(i))/dir(i)
              aux2 = (bounds(2,i)-x0(i))/dir(i)
              tuaux(i) = max(aux1,aux2)
              tlaux(i) = min(aux1,aux2)
              if (aux1*aux2.gt.0.0d0) then
                 write(*,*) 'the reference value of parameter',i
                 write(*,*) 'not inside given bounds'
                 stop
              endif
           endif
        enddo

        if (nevastep(iaux).ge.2) then
          tl     = tlaux(idmax(naux,tlaux,1))
          tu     = tuaux(idmin(naux,tuaux,1))
          deltat = (tu-tl)/(nevastep(iaux)-1)
        else
          write(*,*) 'nevastep  = ', nevastep(iaux)
          write(*,*) 'n of simulation steps must be at least 2'
          stop
        endif

*       compute points inside of "box" (and corresponding values of tt) until
*       optimal point

        do i=1,nevastep(iaux)
           tt  = tl+(i-1)*deltat
           do j=1,naux
              xdir(j,i) = x0(j) + tt * dir(j)
           enddo
*          take the values for coordinate lines
           if (iaux .le. ncodir) then
              coline(i) = xdir(iaux,i)
           endif
        enddo


      return
      end


