      subroutine simulate(ar,nr,ai,ni,
     &                    xsim,nsim,lsim,
     &                    dir,ndir,nevastep,
     &                    xdir,coline,moststep,nfilepoint,
     &                    nobs,nsets,
     &                    states,states0,nstatea,mostobs,nsaveset,
     &                    xdata,nx,yest,mostydat,nydata,
     &                    gpar,ngpar,
     &                    lpar,nlpar)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4 nsim,ndir,moststep,nfilepoint
      integer*4 nsets,nstatea,mostobs,nsaveset
      integer*4 nx,mostydat,ngpar,nlpar
      real*8    xsim(nsim)
      integer*4 lsim(4,nsim)
      real*8    dir(nsim,ndir)
      integer*4 nevastep(ndir)
      real*8    xdir(nsim,moststep)
      real*8    coline(moststep,nsim)
      integer*4 nobs(nsets)
      real*8    states(nstatea,mostobs,nsaveset)
      real*8    states0(nstatea,nsets)
      real*8    xdata(nx,mostobs,nsets)
      real*8    yest(mostydat,mostobs,nsaveset)
      integer*4 nydata(mostobs,nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)

*     LOCAL VARIABLES

      integer*4 nresultfile,isim,i,j,jj,k,k1,ndir1

      include 'common3.inc'

*     simulations along colines or in eig/file directions only:

*     open file to save the computed states

      if (resultfile(1:5).ne.'dummy') then

       nresultfile    = 3
       open(nresultfile,file=resultfile,status='unknown')

*      first write some scalars and matrixes

       write(nresultfile,*)  nsets
       write(nresultfile,*)  nx
       write(nresultfile,*)  nstatea
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

       write(nresultfile,*) nsim
       write(nresultfile,*) ndir
       write(nresultfile,*) dncodir
       write(nresultfile,*) dneigdir
       write(nresultfile,*) nfilepoint

      if (dnsim.eq.0) then

*       simulation results at one given variable point only

          write(nresultfile,*) jacout

          do k = 1,nsets
             iset = k
             k1   = mod(k,nsaveset) + 1
             iobs = -1                    ! flag: states for all obss

             if(debug .gt. 10) write(6,*) 'calling cstates, write'
             dumpcall = 1
             call cstates(ar,nr,ai,ni,xsim,dnsim,lsim,
     &                    states,nstatea,mostobs,nsaveset)
             dumpcall = 0
            do j = 1,nobs(k)

               if(debug.gt.10)write(6,*)'calling observations,simulate'

               call observations(states(1,j,k1),nstatea-dnsaux,
     &                        yest(1,j,k1),nydata(j,k),
     &                        xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &                        states0(1,k),
     &                        gpar(1),ngpar,
     &                        lpar(1,k),nlpar,
     &                        j,k)

             enddo  !j, nobs

            do j = 1,nobs(k)
                do i = 1,nstatea
                   write(nresultfile,*) real(states(i,j,k1))
                enddo
            enddo

            do j = 1,nobs(k)
               do i = 1,nydata(j,k)
                 write(nresultfile,*) real(yest(i,j,k1))
               enddo  ! i, nydata
            enddo


          enddo   !k, nsets


      elseif (dnsim.gt.0) then   !'target' given

        if (ndir.le.0) then
           write(*,*) 'no simulation directions given'
           stop
        endif

*       simulation results at several variable points, in the directions
*       given by 'xsim', 'bounds' and 'dir'

        do i = 1,ndir
          write(nresultfile,*) nevastep(i)
        enddo

        ndir1 = dncodir + dneigdir  !only - filepoints below
        if (ndir1.gt.0) call direct(ar,nr,ai,ni,
     &             ai(pcodir),dncodir,
     &             ai(peigdir),dneigdir,
     &             dir,ndir1,nsim,dnest,
     &             ar(pjtj),ar(pjtj1),ar(peigval),ar(peigvec))

******  compute the states on the lines given by dir,xsim and bounds

        do isim = 1,ndir1

         call dirline(isim,nsim,ndir,dir(1,isim),xsim,ar(psbounds),
     &               ai(pnsimste),xdir,moststep,
     &               coline(1,min(isim,nsim)),dncodir,
     &               ar(ptlaux),ar(ptuaux))


         if (isim.le.dncodir) then
******    write individual coordinate lines
          do jj = 1,nevastep(isim)+1
           write(nresultfile,*) real(coline(jj,isim))
          enddo

         endif

         do jj = 1,nevastep(isim)+1

         if (isim.gt.dncodir) then
******    write  coordinate vectors at each point of line
          do i = 1,nsim
            write(nresultfile,*) real(xdir(i,jj))
          enddo
         endif

          do k = 1,nsets
             iset = k
             k1   = mod(k,nsaveset) + 1
             iobs = -1                    ! flag: states for all obss

             if(debug .gt. 10) write(6,*) 'calling cstates, write'
             dumpcall = 1
             call cstates(ar,nr,ai,ni,xdir(1,jj),nsim,lsim,
     &                    states,nstatea,mostobs,nsaveset)
             dumpcall = 0
             do j = 1,nobs(k)

               if(debug.gt.10)write(6,*)'calling observations,simulate'

               call observations(states(1,j,k1),nstatea-dnsaux,
     &                        yest(1,j,k1),nydata(j,k),
     &                        xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &                        states0(1,k),
     &                        gpar(1),ngpar,
     &                        lpar(1,k),nlpar,
     &                        j,k)

             enddo  !j, nobs

             do j = 1,nobs(k)
                do i = 1,nstatea
                   write(nresultfile,*) real(states(i,j,k1))
                enddo
             enddo

             do j = 1,nobs(k)
               do i = 1,nydata(j,k)
                 write(nresultfile,*) real(yest(i,j,k1))
               enddo  ! i, nydata
             enddo


          enddo   !k, nsets
         enddo  !jj
        enddo  !isim, ndir1
       endif !isim

       if (nfilepoint.gt.0) then

*      read arbitrary simulation variable points directly from
*      a matrix in file 'evalfile'

         if (evalfile(1:5) .ne. 'dummy') then
            open(unit=1,file = evalfile, status = 'old')
            do j = 1,nfilepoint
               read(1,*) (xdir(i,j), i=1,nsim)
            enddo
            close(1)

         else
            write(*,*) 'no evalfile file name given'
            stop
         endif

         do jj = 1,nfilepoint

******    write  coordinate vectors at each point of line
          do i = 1,nsim
            write(nresultfile,*) real(xdir(i,jj))
          enddo

          do k = 1,nsets
             iset = k
             k1   = mod(k,nsaveset) + 1
             iobs = -1                    ! flag: states for all obss

             if(debug .gt. 10) write(6,*) 'calling cstates, write'
             call cstates(ar,nr,ai,ni,xdir(1,jj),nsim,lsim,
     &                    states,nstatea,mostobs,nsaveset)

             do j = 1,nobs(k)

                if(debug.gt.10)write(6,*)'calling observations,simulate'
                call observations(states(1,j,k1),nstatea-dnsaux,
     &                         yest(1,j,k1),nydata(j,k),
     &                         xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &                         states0(1,k),
     &                         gpar(1),ngpar,
     &                         lpar(1,k),nlpar,
     &                         j,k)

              enddo  !j, nobs

              do j = 1,nobs(k)
                 do i = 1,nstatea
                    write(nresultfile,*) real(states(i,j,k1))
                 enddo
              enddo

              do j = 1,nobs(k)
                do i = 1,nydata(j,k)
                  write(nresultfile,*) real(yest(i,j,k1))
                enddo  ! i, nydata
              enddo

          enddo   !k, nsets

         enddo  !jj

************ states computed & dumped *****************************
       close(nresultfile)

       endif ! nfilepoint

      else
        write(*,*) 'no file for simulation results given'
        stop
      endif !resultfile

      return
      end

