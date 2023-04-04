C $Id: mdstlsqs.f,v 1.5 2009/02/10 11:38:11 mjlaine Exp $
C ------------------------------------------------------------------------
C mdstmcmc library 
C File: mdstlsqs.f90
C Purpose: sum of squares function, vectorized version
C
C Marko Laine <Marko.Laine@Helsinki.fi>
C ------------------------------------------------------------------------
C $Author: mjlaine $
C $Revision: 1.5 $
C $Date: 2009/02/10 11:38:11 $
C ------------------------------------------------------------------------
C $Log: mdstlsqs.f,v $
C Revision 1.5  2009/02/10 11:38:11  mjlaine
C precision checks on constants
C
C Revision 1.4  2008/11/17 14:34:38  mjlaine
C gfortran -Wall cleanups, array bound in mdstcov.F
C
C Revision 1.3  2007/09/11 19:28:54  mjlaine
C sstrans and sstype
C
C Revision 1.1  2003/02/28 11:22:04  mjlaine
C Initial revision
C
C ------------------------------------------------------------------------
C
C f90 but fixed format
C this is just modlib/lsq.f but vectorised for ny
C
C called in mdstcov.F mcmc.F90
C
      function mdstlsqs(ar,nr,ai,ni,xaux,naux,ny,sstype,sstrans)

*     Interphase routine between optimizer & real 'dimensional' lsq LSQD.

*     ARGUMENTS
      implicit none
      integer*4 ny
      real*8         mdstlsqs(ny)
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)

      integer*4      naux               ! number of variables in Xopt
      real*8         xaux(naux+1)       ! variables to be optimized
c      real*8         lsqds(ny)
c      external       lsqds
      integer*4      sstype
      real*8         sstrans
c      logical        dumpit

      interface
         function lsqds(ar,nr,ai,ni,
     &        xaux,naux,ny,
     &        xdata,nx,mostobs,nsets,
     &        ydata,mostydat,
     &        weight,
     &        nobs,
     &        nydata,
     &        states,nstatea,nsaveset,
     &        yest,
     &        states0,
     &        gpar,ngpar,
     &        lpar,nlpar,sstype,sstrans)

         implicit none
         integer*4    ny
         real*8       lsqds(ny)
         integer*4    nr,ni
         real*8       ar(nr)
         integer*4    ai(ni)
         integer*4    naux,nsets,nx,mostobs,mostydat
         integer*4    nstatea,nsaveset,ngpar,nlpar
         real*8       xaux(naux) ! variables to be optimized
         real*8       xdata(nx,mostobs,nsets)
         real*8       ydata(mostydat,mostobs,nsets)
         real*8       weight(mostydat,mostobs,nsets)
         integer*4    nobs(nsets)
         integer*4    nydata(mostobs,nsets)
         real*8       states(nstatea,mostobs,nsaveset)
         real*8       yest(mostydat,mostobs,nsaveset)
         real*8       states0(nstatea,nsets)
         real*8       gpar(ngpar)
         real*8       lpar(nlpar,nsets)
         integer*4    sstype
         real*8       sstrans
c         logical      dumpit
         end function lsqds
      end interface



*     COMMON BLOCKS

      include 'common3.inc'

c pitääkö tämä tehdä????
      call dcopy(dnest,xaux,1,ar(pxest),1)

      mdstlsqs = LSQDs(ar,nr,ai,ni,xaux,naux,ny,
     & ar(pxdata),dnx,dmostobs,dnsets,
     & ar(pydata),dmostyda,
     & ar(pweight),
     & ai(pnobs),
     & ai(pnydata),
     & ar(pstates),dnstatea,dnsavese,
     & ar(pyest),
     & ar(pstates0),
     & ar(pgpar),max(dngpar,1),
     & ar(plpar),max(dnlpar,1),sstype,sstrans)


      return
      end

      function lsqds(ar,nr,ai,ni,
     & xaux,naux,ny,
     & xdata,nx,mostobs,nsets,
     & ydata,mostydat,
     & weight,
     & nobs,
     & nydata,
     & states,nstatea,nsaveset,
     & yest,
     & states0,
     & gpar,ngpar,
     & lpar,nlpar,sstype,sstrans)

*     This is a MODEST-function CALLED BY the selected optimizer (e.g.
*     simflex). LSQ MUST BE specified as an EXTERNAL function in the
*     calling routine.

*     SUBROUTINES NEEDED

*     - Xparser
*     - Constraints
*     - Penalties
*     - States
*     - Observations

*     PURPOSE

*     Returns the sum of squares of the residuals i.e.

*                    SUM (Yobs - Yest)**2
*                  over all
*                observations

*     PARAMETERS AND COMMON BLOCKS



*     ARGUMENTS
      implicit none
      integer*4    ny
      real*8       lsqds(ny)
      integer*4    nr,ni
      real*8       ar(nr)
      integer*4    ai(ni)

      integer*4    naux,nsets,nx,mostobs,mostydat
      integer*4    nstatea,nsaveset,ngpar,nlpar
      real*8       xaux(naux)       ! variables to be optimized
      real*8       xdata(nx,mostobs,nsets)
      real*8       ydata(mostydat,mostobs,nsets)
      real*8       weight(mostydat,mostobs,nsets)
      integer*4    nobs(nsets)
      integer*4    nydata(mostobs,nsets)
      real*8       states(nstatea,mostobs,nsaveset)
      real*8       yest(mostydat,mostobs,nsaveset)
      real*8       states0(nstatea,nsets)
      real*8       gpar(ngpar)
      real*8       lpar(nlpar,nsets)
      integer*4    sstype
      real*8       sstrans
c      logical      dumpit

*     LOCAL VARIABLES

      real*8         e                  ! residual
c      real*8         penalty            ! penalty
c      real*8         dnrm2,lsq2         ! see BLAS/IMSL
      integer*4      i,j,k,j1,k1        ! loop indices
c      integer*4      ndata              ! n of data in dataset (obs)
      integer*4      icount             ! count for calls

      real*8     aux

      integer*4 nobsall

      include 'common3.inc'
      data  icount   /0/                ! initalization

*     THE ALGORITHM

      lsqds   = 0.0d0
      icount = icount + 1

      if (model(1:3) .eq. 'ode' .and. obss0.eq.0) then
         j1  = 2
      else
         j1  = 1
      endif

c     calculate total number of observations
      nobsall=0
      do k=1,nsets
         nobsall = nobsall + nobs(k)
      enddo

      do k=1,nsets
         iset = k
         k1   = mod(k,nsaveset)+1
         iobs = -1                    !flag: all observations


         call sstate(ar,nr,ai,ni,xaux,naux,
     &               states,nstatea,mostobs,nsaveset)

*        Calculate the estimated observations Yest, given by the user given
*        observation function

         do j=j1,nobs(k)
            iobs = j

            if(debug .gt. 10) write(*,*) ' calling Observations'

            call observations(states(1,j,k1),nstatea-dnsaux,
     &                        yest(1,j,k1),nydata(j,k),
     &                        xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &                        states0(1,k),
     &                        gpar(1),ngpar,
     &                        lpar(1,k),nlpar,
     &                        j,k)


c dump yest if requested
c            if (dumpit) then
c               call dumpyest(yest(1,1,k1),mostydat,nobs(k),iset,nobsall)
c            endif


c just checking
c            write(*,*) 'nydata=',nydata(1,1)
c            open(1,file='yest.dat',status='unknown')
c            open(2,file='w.dat',status='unknown')
c            open(3,file='mxdata.dat',status='unknown')
c            open(4,file='mydata.dat',status='unknown')
c            do j1=1,nobs(1)
c               write(1,*) (yest(i,j1,1),i=1,nydata(j1,1))
c               write(2,*) (weight(i,j1,1),i=1,nydata(j1,1))
c               write(3,*) (xdata(i,j1,1),i=1,nx)
c               write(4,*) (ydata(i,j1,1),i=1,nydata(j1,1))
c            enddo            
c            close(1)
c            close(2)
c            close(3)
c            close(4)
c            stop 'pysahdytaan'

*           Update the sum of squares

C$$$            if (usew.eq.1) then               
C$$$               do i=1,nydata(j,k)
C$$$                  if (sstype.eq.0) then
C$$$                     e = yest(i,j,k1) - ydata(i,j,k)
C$$$                  elseif (sstype.eq.1) then
C$$$                     e = dsqrt(yest(i,j,k1)) - dsqrt(ydata(i,j,k))
C$$$                  elseif (sstype.eq.2) then
C$$$                     e = dlog(yest(i,j,k1)) - dlog(ydata(i,j,k))
C$$$                  else
C$$$                     stop 'STOP wrong sstype'
C$$$                  endif
C$$$                  aux  = weight(i,j,k) * (e*e)
C$$$                  lsqds(i) = lsqds(i) + aux
C$$$               enddo
C$$$            else
C$$$               do i=1,nydata(j,k)
C$$$                  if (sstype.eq.0) then
C$$$                     e = yest(i,j,k1) - ydata(i,j,k)
C$$$                  elseif (sstype.eq.1) then
C$$$                     e = dsqrt(yest(i,j,k1)) - dsqrt(ydata(i,j,k))
C$$$                  elseif (sstype.eq.2) then
C$$$                     e = dlog(yest(i,j,k1)) - dlog(ydata(i,j,k))
C$$$                  else
C$$$                     stop 'STOP wrong sstype'
C$$$                  endif
C$$$                  lsqds(i)   = lsqds(i) + e*e
C$$$               enddo
C$$$            endif


c 15.7.2007: use sstrans, 27.8.2007: Poisson flag -999
c 31.8.2007: now both sstype and sstrans
            do i=1,nydata(j,k)
               
               if (sstype .eq. 1) then ! Gaussian with trans
                  if (sstrans .eq. 1.0d0 ) then
                     e = yest(i,j,k1) - ydata(i,j,k)
                  elseif (sstrans .eq. 0.0d0 ) then ! lognormal
                     e = dlog(yest(i,j,k1)) - dlog(ydata(i,j,k))
                  else
                     e = yest(i,j,k1)**sstrans - ydata(i,j,k)**sstrans
                  endif
                  if (usew.eq.1) then
                     lsqds(i) = lsqds(i) + weight(i,j,k) * (e*e)
                  else
                     lsqds(i) = lsqds(i) + e*e
                  endif

               elseif (sstype .eq. 2) then ! ??
                  ! this is lognormal in mcmc.F90
                  ! sstype == 1 .and. sstrans == 0.0
                  write(*,*) 'sstype == 2'
                  stop

               elseif (sstype .eq. 3) then ! Poisson
                  e = yest(i,j,k1)-ydata(i,j,k1)*dlog(yest(i,j,k1))
                  if (usew.eq.1) then
                     lsqds(i) = lsqds(i) + 2.0d0*weight(i,j,k)*e
                  else
                     lsqds(i) = lsqds(i) + 2.0d0*e
                  endif

               elseif (sstype .eq. 4) then ! Binomial
! how to get n here?
!                  e = ydata(i,j,k1)*dlog(yest(i,j,k1)) 
!                  e = e + (n-ydata(i,j,k1))*dlog(1-yest(i,j,k1))
                  write(*,*) 'Binomial not ready yet'
                  stop
                  if (usew.eq.1) then
                     lsqds(i) = lsqds(i) - 2.0d0*weight(i,j,k)*e
                  else
                     lsqds(i) = lsqds(i) - 2.0d0*e
                  endif

               elseif (sstype .eq. 5) then ! t - with given df=sstrans

                  e = yest(i,j,k1) - ydata(i,j,k)
                  if (usew.eq.1) then
                     aux = weight(i,j,k)*e*e
                  else
                     aux = e*e
                  endif
                  aux = (sstrans+1)*dlog(sstrans+aux)
                  lsqds(i) = lsqds(i) + aux

               elseif (sstype .eq. 6) then ! L1, Laplace

                  e = dabs(yest(i,j,k1) - ydata(i,j,k))
                  if (usew.eq.1) then
                     lsqds(i) = lsqds(i) + 2.0d0*weight(i,j,k)*e
                  else
                     lsqds(i) = lsqds(i) + 2.0d0*e
                  endif         

               else
                  write(*,*) 'sstype not defined:',sstype
                  stop
               endif
               
            enddo

         enddo
      enddo

      if (optmonit.gt.0) then
         if (mod(icount,optmonit).eq.0 .and.
     &        task(1:3).eq.'sen'    .and.
     &        objfun(1:3).eq.'lsq') 
     &        write(*,*) 'count and lsq ', icount,lsqds
      endif
c      lsqd = lsqd + penalty(xopt,nopt,ar(penpar(1)),dnpen)

      return

      end
