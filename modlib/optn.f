      real*8 function opt(ar,nr,ai,ni,xopt,nopt)

*     Interphase routine between optimizer & real 'dimensional' optd.

*     ARGUMENTS
      implicit none
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)

      integer*4      nopt               ! number of variables in Xopt
      real*8         xopt(nopt)       ! variables to be optimized
      real*8         optd
      external       optd

*     COMMON BLOCKS

      include 'common3.inc'

      opt = optd(ar,nr,ai,ni,xopt,nopt,
     & ar(pstates),dnstates,dnsaux,dnstatea,dnsavese,
     & ar(pyest),
     & ar(pxdata),dnx,dmostobs,dnsets,
     & ar(pydata),dmostyda,
     & ar(pweight),
     & ai(pnobs),
     & ai(pnydata),
     & ar(pstates0),
     & ar(pgpar),max(dngpar,1),
     & ar(plpar),max(dnlpar,1))
      return
      end

      real*8 function optd(ar,nr,ai,ni,xopt,nopt,
     & states,nstates,nsaux,nstatea,nsaveset,
     & yest,
     & xdata,nx,mostobs,nsets,
     & ydata,mostydat,
     & weight,
     & nobs,
     & nydata,
     & states0,
     & gpar,ngpar,
     & lpar,nlpar)

*     This is a MODEST-function CALLED BY the selected optimizer (e.g.
*     simflex). opt MUST BE specified as an EXTERNAL function in the
*     calling routine.

*     SUBROUTINES NEEDED

*     - Xparser
*     - Constraints
*     - Penalties
*     - States
*     - Observations

*     PURPOSE

*     Compute the states, observations etc for the user written
*     routine 'optfun'. There the user may define any objective function
*     that can be expressed in terms of the given variables.

*     PARAMETERS AND COMMON BLOCKS

*     ARGUMENTS
      implicit none
      integer*4    nr,ni
      real*8       ar(nr)
      integer*4    ai(ni)

      integer*4    nopt,nsets,nx,mostobs,mostydat
      integer*4    nstates,nsaux,nstatea,nsaveset,ngpar,nlpar
      real*8       xopt(nopt)       ! variables to be optimized
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


*     local variables
      integer*4    j,k,k1
      integer*4    icount             ! count for calls
      real*8       optfun
      external     optfun

      include 'common3.inc'
      data  icount   /0/                ! initalization
      icount = icount + 1

*     Note: the loop below makes sense if nsets>1 and nsaveset=nsets.
*     Typically in task 'opt', however, nset = 1 = nsaveset.

      do k=1,nsets
         iset = k
         k1   = mod(k,nsaveset)+1
         iobs = -1                    !flag: all observations

         if(debug .gt. 10) write(6,*) 'calling cstates, opt'

         call sstate(ar,nr,ai,ni,xopt,nopt,
     &               states,nstatea,mostobs,nsaveset)

*        Calculate the estimated observations Yest, given by the user given
*        observation function

         do j = 1,nobs(k)
            iobs = j

            if(debug .gt. 10) write(*,*) ' calling Observations'

            call observations(states(1,j,k1),nstates,
     &                        yest(1,j,k1),nydata(j,k),
     &                        xdata(1,1,k),nx,nobs(k),nsaux,nstatea,
     &                        states0(1,k),
     &                        gpar(1),ngpar,
     &                        lpar(1,k),nlpar,
     &                        j,k)


         enddo ! j

         if(debug .gt. 10) write(*,*) 'calling the user written optfun'

         optd = optfun(states(1,1,k1),nstates,
     &                 yest(1,1,k1),nydata(1,k),
     &                 xdata(1,1,k),nx,nobs(k),nsaux,nstatea,
     &                 ydata,weight,mostydat,
     &                 states0(1,k),
     &                 gpar(1),ngpar,
     &                 lpar(1,k),nlpar,
     &                 k)


      enddo ! k

      if (optmonit.gt.0) then
         if (mod(icount,optmonit).eq.0 .and.
     &        task(1:3).eq.'sen'    .and.
     &        objfun(1:3).eq.'opt')
     &        write(*,*) 'count and opt ', icount,optd
      endif

      return
      end

