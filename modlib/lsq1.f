      real*8 function LSQ(ar,nr,ai,ni,xaux,naux)

*     Interphase routine between optimizer & real 'dimensional' lsq LSQD.

*     ARGUMENTS
      implicit none
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)

      integer*4      naux               ! number of variables in Xopt
      real*8         xaux(naux+1)       ! variables to be optimized
      real*8         lsqd
      external       lsqd

*     COMMON BLOCKS

      include 'common3.inc'

      lsq = LSQD(ar,nr,ai,ni,xaux,naux,
     & ar(pxdata),dnx,dmostobs,dnsets,
     & ar(pydata),dmostyda,
     & ar(pweight),
     & ai(pnobs),
     & ai(pnydata),
     & ar(pstates),dnstatea,dnsavese,
     & ar(pyest),
     & ar(pstates0),
     & ar(pgpar),max(dngpar,1),
     & ar(plpar),max(dnlpar,1))


      return
      end

      real*8 function LSQD(ar,nr,ai,ni,
     & xaux,naux,
     & xdata,nx,mostobs,nsets,
     & ydata,mostydat,
     & weight,
     & nobs,
     & nydata,
     & states,nstatea,nsaveset,
     & yest,
     & states0,
     & gpar,ngpar,
     & lpar,nlpar)

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

*                    SUM abs((Yobs - Yest))
*                  over all
*                observations

*     PARAMETERS AND COMMON BLOCKS



*     ARGUMENTS
      implicit none
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

*     LOCAL VARIABLES

      real*8         e                  ! residual
      real*8         penalty            ! penalty
      real*8         dnrm2,lsq2         ! see BLAS/IMSL
      integer*4      i,j,k,j1,k1        ! loop indices
c      integer*4      ndata              ! n of data in dataset (obs)
      integer*4      icount             ! count for calls
      integer*4      irmeth
      real*8         lpp,c
      real*8     aux,ww

      include 'common3.inc'
      data  icount   /0/                ! initalization

*     THE ALGORITHM

      lsqd   = 0.0d0
      if (icount.eq.0) then
        open(1,file='robust.met',status = 'unknown')
        read(1,*) irmeth,lpp
        close(1)
      endif

      icount = icount + 1

      if (model(1:3) .eq. 'ode' .and. obss0.eq.0) then
         j1  = 2
      else
         j1  = 1
      endif

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

*           Update the sum of squares

          if (irmeth.eq.1) then  ! Lp
            if (usew.eq.1) then
               do i=1,nydata(j,k)
                  e = (yest(i,j,k1) - ydata(i,j,k))*weight(i,j,k)
                  e = dabs(e)
                  aux  = e **lpp
                  lsqd = lsqd + aux
               enddo
            else
               do i=1,nydata(j,k)
                  e = yest(i,j,k1) - ydata(i,j,k)
                  e = dabs(e)
                  lsqd   = lsqd + e**lpp
               enddo
            endif
          elseif (irmeth.eq.2) then !fair
            if (usew.eq.1) then
               do i=1,nydata(j,k)
                  e = (yest(i,j,k1) - ydata(i,j,k))*weight(i,j,k)
                  e = dabs(e)
                  c = lpp
                  aux  = c*c*(e/c - log(1.0d0+e/c))
                  lsqd = lsqd + aux
               enddo
            else
               do i=1,nydata(j,k)
                  e = yest(i,j,k1) - ydata(i,j,k)
                  e = dabs(e)
                  c = lpp
                  aux  = c*c*(e/c - log(1.0d0+e/c))
                  lsqd   = lsqd + aux
               enddo
            endif
          elseif (irmeth.eq.3) then !tukey
            if (usew.eq.1) then
               do i=1,nydata(j,k)
                  e = (yest(i,j,k1) - ydata(i,j,k))*weight(i,j,k)
                  e = dabs(e)
                  c = lpp
                  if (e.le.c) then
                   aux = c*c/6.0*(1.0d0-(1.0d0-(e/c)**2)**3)
                   lsqd = lsqd + aux
                  elseif (e.gt.c) then
                   aux = c*c/6.0d0
                   lsqd = lsqd + aux
                  endif
               enddo
            else
               do i=1,nydata(j,k)
                  e = yest(i,j,k1) - ydata(i,j,k)
                  e = dabs(e)
                  c = lpp
                  if (e.le.c) then
                   aux = c*c/6.0*(1.0d0-(1.0d0-(e/c)**2)**3)
                   lsqd = lsqd + aux
                  elseif (e.gt.c) then
                   aux = c*c/6.0d0
                   lsqd = lsqd + aux
                  endif
               enddo
            endif

          endif !irmeth


         enddo
      enddo

      if (optmonit.gt.0) then
         if (mod(icount,optmonit).eq.0 .and.
     &        task(1:3).eq.'sen'    .and.
     &        objfun(1:3).eq.'lsq') 
     &        write(*,*) 'count and lsq ', icount,lsqd
      endif
c      lsqd = lsqd + penalty(xopt,nopt,ar(penpar(1)),dnpen)

      return

      end

