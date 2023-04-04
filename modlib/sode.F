      subroutine sode(ar,nr,ai,ni,
     &                s, s0, states, nstatea,
     &                yest,mostydat,mostobs,nydata,
     &                nsaveset, xdata, nx, nsets, nobs,
     &                gpar, ngpar, lpar, nlpar)

*     CALLED BY

*     cstates

*     SUBROUTINES NEEDED

*     - Ode solvers

*     PURPOSE

*     This subroutine calculates the states in the user given points (in x)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4 nstatea,mostydat,mostobs,nsaveset,nx,nsets
      integer*4 ngpar,nlpar
      real*8    s(nstatea)                       ! s is actually the whole ar/VMT 24.11.93
      real*8    s0(nstatea,nsets)
      real*8    states(nstatea,mostobs,nsaveset)
      real*8    yest(mostydat,mostobs,nsaveset)
      integer*4 nydata(mostobs,nsets)
      real*8    xdata(nx,mostobs,nsets)
      integer*4 nobs(nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)

*     LOCAL VARIABLES
      real*8    t0aux                              ! start of integration
      integer*4 k,k1,kk,jj,i,dumpc                 ! loop indices,etc
      integer*4 nstates

      logical   dumpbin

*     COMMON BLOCKS

      include 'common3.inc'
      data dumpc /1/
      data dumpbin /.false./

      if(debug .gt. 11) write(*,*) ' in sode'

      nstates = nstatea - dnsaux                 ! dnsaux comes from common3.inc

*     Variable iset tells the current data set and may
*     be used in inits0 and Fode; it comes from common3.inc.

      k  = iset
      k1 = mod(k,nsaveset)+1

*     Get initial state values and start time for integration
*     in current data set.

      if(debug .gt. 10) write(*,*) ' calling inits0'

      iobs = 1
      call inits0(nstates,t0aux,s,
     &     xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &     s0(1,k),
     &     gpar,max(dngpar,1),
     &     lpar(1,k),max(dnlpar,1),
     &     iobs,iset)

*     Auxiliary states at inital point:
      if (dnsaux.gt.0) call ddrFode(nstates,t0aux,s,states(1,1,k1))

      call dcopy(nstatea,s,1,states(1,1,k1), 1)

      if(dumpcall .eq. 1 .and. dumpfile(1:6) .ne. 'nodump') then
          if (dumpc.eq.1) then
             if (dumpfile(len_trim(dumpfile)-3:len_trim(dumpfile)) 
     &             == '.bin') then
                dumpbin = .true.
             else
                dumpbin = .false.
             end if

             if (dumpbin) then

#ifdef __GFORTRAN__
                open(unit=11, file=dumpfile, status='replace', 
     &               form='unformatted', position='rewind', 
     &               access='stream')
#else
                open(unit=11, file=dumpfile, status='replace', 
     &               form='binary', position='rewind')
#endif
                dumpc = dumpc+1
                
                write(11)  nsets
                write(11)  nstatea
                write(11)  dnsaux
                write(11)  ndumpp
                do kk = 1,nsets
c                  do jj = 1,nobs(kk)
                   write(11) nydata(1,kk)
c                  enddo
                enddo

             else
                open(11,file=dumpfile,status='unknown')
                dumpc = dumpc+1
                
                write(11,*)  nsets
                write(11,*)  nstatea
                write(11,*)  dnsaux
                write(11,*)  ndumpp
                do kk = 1,nsets
c                  do jj = 1,nobs(kk)
                   write(11,*) nydata(1,kk)
c                  enddo
                enddo
             end if
c          else
c             open(11,access='APPEND',file=dumpfile,status='unknown')
          endif

          if(debug.gt.10) write(6,*)'calling Observations'
          call observations(states(1,1,k1),nstates,
     &             yest(1,1,k1),nydata(1,k),
     &             xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &             s0(1,k),
     &             gpar(1),ngpar,
     &             lpar(1,k),nlpar,
     &             1,k)

          if (dumpbin) then
             write(11) t0aux
             do i = 1,nydata(1,k)
                write(11) yest(i,1,k1)
             enddo
             do i = 1,nstatea
                write(11) s(i)
             enddo
          else
             write(11,*) t0aux
             do i = 1,nydata(1,k)
                write(11,*) yest(i,1,k1)
             enddo
             do i = 1,nstatea
                write(11,*) s(i)
             enddo
          endif

      endif

      if (odesolver(1:4) .eq. 'odes') then

*         Odessa is used.
          call odeode(ar,nr,ai,ni,
     &                s,s0,states,nstatea,nstates,
     &                yest,mostydat, mostobs,nydata,
     &                nsaveset, xdata, nx, nsets, nobs, t0aux,
     &                gpar, ngpar, lpar, nlpar,
     &                ar(prworkod),dlrwrkod,ai(piworkod),dliwrkod,
     &                ar(prtolod),dlrtolod,ar(patolod),dlatolod,
     &                ar(pteta0od),max(dnest,1),ar(pxaux))
c added 24.10.1999 ML
#ifdef USE_LSODE
      else if (odesolver(1:4) .eq. 'lsod') then
          call lodeode(ar,nr,ai,ni,
     &                s,s0,states,nstatea,nstates,
     &                yest,mostydat, mostobs,nydata,
     &                nsaveset, xdata, nx, nsets, nobs, t0aux,
     &                gpar, ngpar, lpar, nlpar,
     &                ar(prworkod),dlrwrkod,ai(piworkod),dliwrkod,
     &                ar(prtolod),dlrtolod,ar(patolod),dlatolod,
     &                max(dnest,1),ar(pxaux))
#endif
      else if (odesolver(1:4) .eq. 'ddri') then

*         Ddriv2 is used.
c          call ddrode(s,s0,states,nstatea,nstates,
c     &                yest,mostydat, mostobs,nydata,
c     &                nsaveset, xdata, nx, nsets, nobs, t0aux,
c     &                gpar, ngpar, lpar, nlpar,
c     &                ar(pworkddr),dlrwddr,ai(piwrkddr),dliwddr)
           write(*,*) 'ddriv2 Ode solver not implemented'
           stop

      else if (odesolver(1:4) .eq. 'eule') then

*          Euler is used.
           call eulode(s,s0,states,nstatea,nstates,
     &                 yest,mostydat, mostobs,nydata,
     &                 nsaveset, xdata, nx, nsets, nobs, t0aux,
     &                 gpar, ngpar, lpar, nlpar,
     &                 ar(pworkeul))


      else if (odesolver(1:4) .eq. 'divp') then

c       Divpag /IMSL used
c       If necessary, rebuild from ddriv interface.
        write(*,*) 'divpag is not implemented (template in sodebup.f)'
        stop

      else
         write(*,*) 'The odesolver not specified!'
      end if

      iobs = -1           ! all observations in a set have been processed

c      close(11)

      return
      end

