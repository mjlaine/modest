      real*8 function R2(ar,nr,ai,ni,xaux,naux)

*     Interphase routine between optimizer & real 'dimensional'  R2D.

*     ARGUMENTS
      implicit none
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)

      integer*4      naux               ! number of variables in Xopt
      real*8         xaux(naux+1)       ! variables to be optimized
      real*8         R2d
      external       R2d

*     COMMON BLOCKS

      include 'common3.inc'

      R2 = R2D(ar,nr,ai,ni,xaux,naux,
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

      real*8 function R2D(ar,nr,ai,ni,
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
*     simflex). R2 MUST BE specified as an EXTERNAL function in the
*     calling routine.

*     SUBROUTINES NEEDED

*     - Xparser
*     - Constraints
*     - Penalties
*     - States
*     - Observations
*     - lsq

*     PURPOSE

*     Returns the R2 value of the residuals i.e.

*                  SUM (1 - (Yobs - Yest)**2/(Yobs - Ymean)**2)
*              over all
*          observations

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
      real*8         dnrm2              ! see BLAS/IMSL
      integer*4      i,j,k,j1,jj,k1        ! loop indices
      integer*4      icount


      Real*8       y                              ! ydata element
      Real*8       weigh                          ! corresponding weight
      Real*8       ss                             ! total corrected ss
      Real*8       ysum                           ! weighted sum of ydata
      Real*8       wtot                           ! sum of weight
      Real*8       sse                            ! residual ss
      Real*8       mse                            ! residual ms
      Real*8       stde                           ! std error of residuals
      real*8       tol
      real*8       sum
      Real*8       lsq                            ! see LSQ
      external     lsq


      include 'common3.inc'

*     THE ALGORITHM

*     calculate the total corrected sum of squares for the data
*     first take the weight into account

      if(debug .gt. 10) write(*,*) ' in R2D'

      data  icount   /0/                ! initalization
      icount = icount + 1

      ss   = 0.0d0
      ysum = 0.0d0
      wtot = 0.0d0

      if (model(1:3) .eq. 'ode' .and. obss0.eq.0) then !obss0(j,k)!!
         jj    = 2
      else
         jj    = 1
      endif

         do k=1,nsets
            do j=jj,nobs(k)
               do i = 1,nydata(j,k)

                  y      = ydata(i,j,k)
                  weigh  = weight(i,j,k)
                  wtot   = wtot + weigh
                  ysum   = ysum + weigh *y
                  ss     = ss + weigh *(y**2)

               enddo
            enddo
         enddo

*     correct ss by taking the mean into account

      ss = ss - (ysum ** 2) / wtot

*     calculate the mean sum of squares of the residuals

      sse  = lsq(ar,nr,ai,ni,xaux,naux)
      R2d  = -(1.d0-sse/ss)*100.

      if (optmonit.gt.0) then
        if (mod(icount,optmonit).eq.0) write(*,*) 'R2 ', -R2d
      endif

      return
      end

