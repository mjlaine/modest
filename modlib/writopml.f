      subroutine writopml(ar,nr,ai,ni,
     &                  xest,nest,lest,
     &                  xexp,nexp,lexp,
     &                  xopt,nopt,lopt,
     &                  xaux,nesopsi,laux,
     &                  nobs,nsets,
     &                  states,states0,nstatea,mostobs,nsaveset,
     &                  xdata,nx,
     &                  ydata,weight,
     &                  yest,xjac,jtj,jtj1,
     &                  mostydat,nydata,
     &                  gpar,ngpar,
     &                  lpar,nlpar)

      implicit none

      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4 nsim,nest,nexp,nsen,nopt,nesopsi
      real*8    xest(nest),xexp(nexp)
      real*8    xopt(nopt),xaux(nesopsi)
      integer*4 lest(4,nest),lexp(4,nexp)
      integer*4 lopt(4,nopt)
      integer*4 laux(4,nesopsi)
      integer*4 nsets,ngpar,nlpar
      integer*4 nstatea,mostobs,nsaveset,mostydat,nx
      integer*4 nobs(nsets)
      real*8    states(nstatea,mostobs,nsaveset)
      real*8    states0(nstatea,nsets)
      real*8    xdata(nx,mostobs,nsets)
      real*8    ydata(mostydat,mostobs,nsets)
      real*8    weight(mostydat,mostobs,nsets)
      real*8    yest(mostydat,mostobs,nsaveset)
      real*8    xjac(mostydat,nest,mostobs,nsaveset)
      real*8    jtj(nest,nest),jtj1(nest,nest)
      integer*4 nydata(mostobs,nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)


*     LOCAL PARAMETERS

      integer*4 i,p,p1,p2,j,k,k1            ! loop indices
      integer*4 naux,nresultfile,noptfile
      logical   casey                       ! false: sim/exp/opt
                                            ! true : est/sen

*     COMMON BLOCS

      include 'common3.inc'

*     write the results to file (all in one column, to be analyzed by
*     Matlab)

      if(debug .gt. 11) write(*,*) ' in writopml'

      if (task(1:3).eq.'est' .or. task(1:3).eq.'sen') then
         casey = .true.
      else
         casey = .false.
      endif

      if (task(1:3).eq.'est') then
         naux = dnest
         do i = 1,naux
            xaux(i) = xest(i)
            call icopy(4,lest(1,i),1,laux(1,i),1)
         enddo
      elseif (task(1:3).eq.'exp') then
         naux = dnexp
         do i = 1,naux
            xaux(i) = xexp(i)
            call icopy(4,lexp(1,i),1,laux(1,i),1)
         enddo
      elseif (task(1:3).eq.'opt') then
         naux = dnopt
         do i = 1,naux
            xaux(i) = xopt(i)
            call icopy(4,lopt(1,i),1,laux(1,i),1)
         enddo
      endif

*     open the files

      nresultfile = 1
      noptfile  = 3

      if (resultfile(1:5).eq.'dummy') then
          open(nresultfile,status='scratch')
      else
          open(nresultfile,file=resultfile,status='unknown')
      endif

      if (optfile(1:5).eq.'dummy') then
         open(noptfile,status='scratch')
      else
         open(noptfile,file=optfile,status='unknown')
      endif

*     first write scalars and matrixes kept in memory

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

*     the 'x' data

      do k = 1,nsets
         do j = 1,nobs(k)
            do i = 1,nx
               write(nresultfile,*) real(xdata(i,j,k))
            enddo
         enddo
      enddo

*     the 'y' data

      if (casey) then
        do k = 1,nsets
           do j = 1,nobs(k)
              do i = 1,nydata(j,k)
                write(nresultfile,*) real(ydata(i,j,k))
              enddo
           enddo
        enddo
      endif

*     write the optimized arguments

      write(nresultfile,*) naux
        do i = 1,naux
           write(nresultfile,*) real(xaux(i))
        enddo
        if (task(1:3).eq.'exp') then
           call xparser(ar,nr,ai,ni,'sel',xest,dnest,lest)
           write(nresultfile,*) dnest
           do i = 1,dnest
              write(nresultfile,*) real(xest(i))
           enddo
      endif


      write(nresultfile,*) jacout
      if (jacout.eq.1) then
         do p1 = 1,nest
            call dcopy(nest,0.0d0,0,jtj(1,p1),1)
         enddo
      endif

*     for IO updates, the estimation result saved independently
      do i = 1,naux
         write(noptfile,*) real(xaux(i))
      enddo

*     then compute & write volatile response variables for each set

      do k = 1,nsets
         iset = k
         k1   = mod(k,nsaveset) + 1
         iobs = -1                    ! flag: states for all obss

         if(debug .gt. 10) write(6,*) 'calling cstates, write'


         if (task(1:3).eq.'est') then
            dumpcall = 1
            call cstates(ar,nr,ai,ni,xest,dnest,lest,
     &                states,nstatea,mostobs,nsaveset)
            dumpcall = 0
         elseif (task(1:3).eq.'exp') then
            dumpcall = 1
            call cstates(ar,nr,ai,ni,xexp,dnexp,lexp,
     &                states,nstatea,mostobs,nsaveset)
            dumpcall = 0
         elseif (task(1:3).eq.'opt') then
            dumpcall = 1
            call cstates(ar,nr,ai,ni,xopt,dnopt,lopt,
     &                states,nstatea,mostobs,nsaveset)
            dumpcall = 0
         endif


         if(debug.gt.10) write(6,*)'calling observations, write'

         do j = 1,nobs(k)
            call observations(states(1,j,k1),nstatea-dnsaux,
     &                        yest(1,j,k1),nydata(j,k),
     &                        xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &                        states0(1,k),
     &                        gpar(1),ngpar,
     &                        lpar(1,k),nlpar,
     &                        j,k)
         enddo

*        computed states

            do j = 1,nobs(k)
               do i = 1,nstatea
                  write(nresultfile,*) real(states(i,j,k1))
               enddo
            enddo

*        computed observables

            do j = 1,nobs(k)
               do i = 1,nydata(j,k)
                  write(nresultfile,*) real(yest(i,j,k1))
               enddo
            enddo

*        the Jacobian

           if (jacout .eq. 1) then

              call jacob(ar,nr,ai,ni,
     &                xest,nest,lest,
     &                xjac,jtj1,xaux,mostydat,
     &                xdata,nx,mostobs,nsets,
     &                weight,
     &                nobs,nydata,
     &                states,nstatea,nsaveset,
     &                yest,states0,
     &                gpar,ngpar,
     &                lpar,nlpar)

*             write the Jacobian

              do p = 1,nest
                  do j = 1,nobs(k)
                    do i = 1,nydata(j,k)
                      write(nresultfile,*) real(xjac(i,p,j,k1))
                    enddo
                  enddo
              enddo

*            update the approx. Hessian

              do p1 = 1, nest
                do p2 = 1, nest
                   jtj(p1,p2) = jtj(p1,p2) + jtj1(p1,p2)
                enddo
              enddo

           endif

      enddo

*     finally, write the cumulatively computed Hessian and
*     standard statistics for parameter estimates

          if (jacout.eq.1) then
             do i = 1,nest
                do j = 1,nest
                 write(nresultfile,*) real(jtj(i,j))
                enddo
             enddo
          endif

*       compute the standard statistics

         if(task(1:3) .eq. 'est' .and. stats .eq. 1) then
            if(debug .gt. 10) write(*,*) ' calling cstats'
            call cstats(ar,nr,ai,ni,
     &                  xest,nest,lest,
     &                  nsets,nobs,nydata,
     &                  ydata,weight,mostydat,mostobs,
     &                  jtj,jtj1,jtj1,ar(peigval),ar(peigvec))

         endif

      return
      end
