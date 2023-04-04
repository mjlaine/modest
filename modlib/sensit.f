
      subroutine sensit(ar,nr,ai,ni,
     &                  xsen,nsen,lsen,
     &                  nobs,nsets,mostobs,
     &                  xdata,nx,nydata,
     &                  dir,ndir,nevastep,
     &                  xdir,coline,moststep,nfilepoint,
     &                  optline,optgrid,igrid,ncontour)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4 nsen,nsets,mostobs,nx,ndir
      integer*4 nfilepoint,moststep,ncontour
      real*8    xsen(nsen)
      integer*4 lsen(4,nsen)
      integer*4 nobs(nsets)
      real*8    xdata(nx,mostobs,nsets)
      integer*4 nydata(mostobs,nsets)
      real*8    dir(nsen,ndir)
      real*8    xdir(nsen,moststep)
      real*8    coline(moststep,nsen)
      integer*4 nevastep(ndir)
      real*8    optline(moststep,ndir)
      real*8    optgrid(ncontour,moststep,moststep)
      integer*4 igrid(2,ncontour)


*     LOCAL VARIABLES

      integer*4 i,j,k,i1,i2,isen,ii,jj,kk,nresultfile,ndir1
      logical   openfile
      real*8    lsq,expdes,opt,r2
      external  lsq,expdes,opt,r2

      include 'common3.inc'

*     calculations along colines or in eig/file directions only:


*     open file to save the computed objective function values?

      nresultfile    = 3
      openfile       = task(1:3).eq.'sen'

      if (openfile) then
       open(nresultfile,file=resultfile,status='unknown')

***************************
*      first write some scalars and matrixes (??!)
       write(nresultfile,*)  nsets
       write(nresultfile,*)  nx
c       write(nresultfile,*)  nstatea
       write(nresultfile,*)  dnsaux

       do k = 1,nsets
          write(nresultfile,*) nobs(k)
          do j = 1,nobs(k)
             write(nresultfile,*) nydata(j,k)
          enddo
       enddo

*      the 'x' data

       do k = 1,nsets
          do j = 1,nobs(k)
             do i = 1,nx
                write(nresultfile,*) real(xdata(i,j,k))
             enddo
          enddo
       enddo

***********************
       write(nresultfile,*) nsen
       write(nresultfile,*) ndir
       write(nresultfile,*) dncodir
       write(nresultfile,*) dneigdir
       write(nresultfile,*) nfilepoint
       write(nresultfile,*) ncontour
       do i = 1,ndir
          write(nresultfile,*) nevastep(i)
       enddo

      endif   ! openfile

      ndir1 = dncodir + dneigdir   !only - filepoints below
      if (ndir1.gt.0) call direct(ar,nr,ai,ni,
     &                     ai(pcodir),dncodir,
     &                     ai(peigdir),dneigdir,
     &                     dir,ndir1,nsen,dnest,
     &                     ar(pjtj),ar(pjtj1),ar(peigval),ar(peigvec))

****** compute the obj function on the lines given by dir,xsen and bounds

       do isen = 1,ndir1

         call dirline(isen,nsen,ndir,dir(1,isen),xsen,
     &                ar(psbounds),nevastep,xdir,moststep,
     &                coline(1,min(isen,dncodir)),dncodir,
     &                ar(ptlaux),ar(ptuaux))

         if (isen.le.dncodir .and. (openfile)) then
******    write individual coordinate lines
          do jj = 1,nevastep(isen)+1
           write(nresultfile,*) real(coline(jj,isen))
          enddo
         endif

         do jj = 1,nevastep(isen)+1

******** compute the values for objective functions ****

              if (objfun(1:3).eq.'lsq' .or.
     &           (task(1:3).eq.'exp' .and. optcrit(1:3).eq.'lsq')) then

                if (debug.gt.10) write(*,*) 'calling lsq in sen '
                optline(jj,isen) = lsq(ar,nr,ai,ni,xdir(1,jj),nsen)
                call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

              elseif (objfun(1:3).eq.'exp') then

                if (debug.gt.10) write(*,*) 'calling expdes in sen'
                optline(jj,isen) = expdes(ar,nr,ai,ni,xdir(1,jj),nsen)
                call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

              elseif (objfun(1:3).eq.'opt') then

                if (debug.gt.10) write(*,*) 'calling opt in sen '
                optline(jj,isen) = opt(ar,nr,ai,ni,xdir(1,jj),nsen)
                call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

              elseif (objfun(1:2).eq.'r2') then

                if (debug.gt.10) write(*,*) 'calling opt in sen '
                optline(jj,isen) = r2(ar,nr,ai,ni,xdir(1,jj),nsen)
                call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

              endif  ! objfun
         enddo  ! jj

         if (openfile) then
            do jj = 1,nevastep(isen)+1
                if (isen.gt.dncodir) then
                  do i = 1,nsen
                     write(nresultfile,*) real(xdir(i,jj))
                  enddo
                endif
                write(nresultfile,*) real(optline(jj,isen))
            enddo
         endif

       enddo  ! isen

*      read arbitrary  variable points directly to 'xdir'
*      from a user given matrix in file 'evalfile'

       if (evalfile(1:5) .ne. 'dummy') then
         sencount = sencount + 1
         if (sencount.eq.1) then
          if (nfilepoint.gt.0) then
            open(unit=1,file = evalfile, status = 'old')
            do j = 1,nfilepoint
               read(1,*) (xdir(i,j), i=1,nsen)
            enddo
            close(1)
          else
            write(*,*) 'nfilepoint not given'
            stop
          endif
         endif

         do jj = 1,nfilepoint

******** compute the values for objective functions ****

              if (objfun(1:3).eq.'lsq' .or.
     &           (task(1:3).eq.'exp' .and. optcrit(1:3).eq.'lsq')) then

                if (debug.gt.10) write(*,*) 'calling lsq in sen '
                optline(jj,isen) = lsq(ar,nr,ai,ni,xdir(1,jj),nsen)
                call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

              elseif (objfun(1:3).eq.'exp') then

                if (debug.gt.10) write(*,*) 'calling expdes in sen'
                optline(jj,isen) = expdes(ar,nr,ai,ni,xdir(1,jj),nsen)
                call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

              elseif (objfun(1:3).eq.'opt') then

                if (debug.gt.10) write(*,*) 'calling opt in sen '
                optline(jj,isen) = opt(ar,nr,ai,ni,xdir(1,jj),nsen)
                call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

              elseif (objfun(1:2).eq.'r2') then

                if (debug.gt.10) write(*,*) 'calling opt in sen '
                optline(jj,isen) = r2(ar,nr,ai,ni,xdir(1,jj),nsen)
                call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

              endif  ! objfun

             if (openfile) then
               do i = 1,nsen
                 write(nresultfile,*) real(xdir(i,jj))
               enddo
               write(nresultfile,*) real(optline(jj,isen))
             endif
           call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)
         enddo  !jj, filepoints

       endif ! evalfile

*     compute the object function values on the contour grids (2-D)

      if (ncontour.gt.0 .and. (openfile)) then
       do kk = 1,ncontour
         write(nresultfile,*) igrid(1,kk)
         write(nresultfile,*) igrid(2,kk)
       enddo
      endif

      do kk = 1,ncontour

        call dcopy(nsen,xsen(1),1,xdir(1,1),1)

          i1 = igrid(1,kk)
          i2 = igrid(2,kk)

          do ii = 1,nevastep(i1)+1
            xdir(i1,1) = coline(ii,i1)
            do jj = 1,nevastep(i2)+1
               xdir(i2,1) = coline(nevastep(i2)-jj+2,i2)

c               if (xdir(i1,1) .eq. xsen(i1)) then
c                   optgrid(kk,jj,ii) = optline(jj,i2)
c               elseif (xdir(i2,1) .eq. xsen(i2)) then
c                   optgrid(kk,jj,ii) = optline(ii,i1)
c               else
                 if (objfun(1:3).eq.'lsq') then

                 if (debug.gt.10) write(*,*) 'calling lsq in sen'
                 optgrid(kk,jj,ii) = lsq(ar,nr,ai,ni,xdir,nsen)
                 call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

                 elseif (objfun(1:3).eq.'exp') then

                  if (debug.gt.10) write(*,*) 'calling expdes in sen'
                  optgrid(kk,jj,ii) = expdes(ar,nr,ai,ni,xdir,nsen)
                  call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

                 elseif (objfun(1:3).eq.'opt') then

                  if (debug.gt.10) write(*,*) 'calling opt in sen'
                  optgrid(kk,jj,ii) = opt(ar,nr,ai,ni,xdir,nsen)
                  call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

                 elseif (objfun(1:2).eq.'r2') then

                  if (debug.gt.10) write(*,*) 'calling opt in sen'
                  optgrid(kk,jj,ii) = r2(ar,nr,ai,ni,xdir,nsen)
                  call Xparser(ar,nr,ai,ni,'par',xsen,nsen,lsen)

                 endif  ! objfun

c               endif  ! xdir

            enddo  ! jj
          enddo  ! ii

          if (openfile) then

            do ii = 1,nevastep(i1)+1
            do jj = 1,nevastep(i2)+1
                write(nresultfile,*) real(optgrid(kk,jj,ii))
            enddo
            enddo

          endif

       enddo  ! kk

      if (openfile) then
        do i=1,nsen
          write(nresultfile,*) real(xsen(i))
        enddo
      endif
      if (openfile) close(nresultfile)

      return
      end


