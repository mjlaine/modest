      real*8 function EXPDES(ar,nr,ai,ni,xaux,naux)

*     Interphase routine between optimizer & real 'dimensional' expdesd.

*     ARGUMENTS
      implicit none
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)

      integer*4      naux               ! number of variables in Xopt
      real*8         xaux(naux)       ! variables to be optimized
      real*8         expdesd
      external       expdesd

*     COMMON BLOCKS

      include 'common3.inc'

      expdes = expdesd(ar,nr,ai,ni,xaux,naux,
     &                 ar(pxest),ai(plest),dnest,
     &                           ai(plexp),dnexp,
     &                 ar(pjtj),ar(pjtj0),ar(pjtj1),ar(pjtj2),
     &                 ar(pexpdu),ar(psigma),ar(pexpdv),
     &                 ai(pnewobs),dnnew,
     &                 ai(pexppar),dnexpp,
     &                 ar(pxdata),dnx,dmostobs,dnsets,
     &                 ar(pydata),dmostyda,
     &                 ai(pnobs),
     &                 ai(pnydata),
     &                 ar(pstates),dnstatea,dnsavese,
     &                 ar(pyest),
     &                 ar(pstates0),
     &                 ar(pgpar),max(dngpar,1),
     &                 ar(plpar),max(dnlpar,1),
     &                 ar(pvaldisc),ddiscpoi)


      return
      end

      real*8 function expdesd(ar,nr,ai,ni,xaux,naux,
     &                        xest,lest,nest,
     &                             lexp,nexp,
     &                        jtj,jtj0,jtj1,jtj2,
     &                        u,s,v,
     &                        newobs,nnew,
     &                        exppar,nexpp,
     &                        xdata,nx,mostobs,nsets,
     &                        ydata,mostydat,
     &                        nobs,
     &                        nydata,
     &                        states,nstatea,nsaveset,
     &                        yest,
     &                        states0,
     &                        gpar,ngpar,
     &                        lpar,nlpar,
     &                        valdisc,ndiscpoint)


*     This is a MODEST-function CALLED BY the selected optimizer (e.g.
*     simflex). The function computes the values of the objective function
*     for optimal design of experiments

*     ARGUMENTS
      implicit none
      integer*4    nr,ni
      real*8       ar(nr)
      integer*4    ai(ni)

      integer*4    naux,nest,nexp,nnew,nexpp
      integer*4    nx,mostydat,mostobs,nsets,nstatea,nsaveset,
     &             ndiscpoint,ngpar,nlpar
      real*8       xaux(naux)       ! variables to be optimized
      real*8       xest(nest)
      integer*4    lest(4,nest),lexp(4,nexp)
      real*8       jtj(nest,nest),jtj0(nest,nest),jtj1(nest,nest)
      real*8       jtj2(nexpp,nexpp)
      real*8       u(nest,nest),s(nest),v(nest,nest)
      integer*4    newobs(nnew,2)
      integer*4    exppar(nexpp)
      real*8       xdata(nx,mostobs,nsets)
      real*8       ydata(mostydat,mostobs,nsets)
      integer*4    nobs(nsets)
      integer*4    nydata(mostobs,nsets)
      real*8       states(nstatea,mostobs,nsaveset)
      real*8       yest(mostydat,mostobs,nsaveset)
      real*8       states0(nstatea,nsets)
      real*8       gpar(ngpar)
      real*8       lpar(nlpar,nsets)
      real*8       valdisc(ndiscpoint)


*     local variables

      integer*4 i,j,p1,p2,k,k1,isen,jj,icount
      real*8    tol
      integer*4 ipath,irank,info
      real*8    sum, var, logdet
      logical   globcrit
      real*8    x0,x1,sig,Desir,R2mean,R2dev

      include 'common3.inc'
      data  icount   /0/                ! initalization

      icount   = icount + 1
      globcrit = optcrit(1:1).eq.'G' .or. optcrit(1:1).eq.'g'

      if(debug .gt. 10) write(*,*) 'calling Xparser (exp)'

      call xparser(ar,nr,ai,ni,'par',xaux,naux,lexp)

      if(debug .gt. 10) write(*,*) 'calling Xparser (est)'

      call xparser(ar,nr,ai,ni,'sel',xest,nest,lest)

      jcount = jcount + 1                !jcount initialized in IO

      if (.not.globcrit) then

        if (jcount.eq.1) then

*     * Initialize the local criteria  ************************
         if(debug .gt. 10) write(*,*) 'calling hess'

         call hess(ar,nr,ai,ni,
     &             jtj,u,
     &             xest,nest,lest)

         if(debug .gt. 10) write(*,*) 'calling hessup'

         call hessup(ar,nr,ai,ni,
     &               jtj1,u,
     &               xest,nest,lest,
     &               newobs,nnew)


         do i = 1,nest
            do j = 1,nest
               jtj0(i,j) = jtj(i,j) - jtj1(i,j)
            enddo
         enddo

       endif ! jcount

      elseif (globcrit) then
                                          ! compute the reference 'data'
          do k = 1,nsets
             iset = k
             k1   = mod(k,nsaveset) + 1
             iobs = -1                    ! flag: states for all obss

             if(debug .gt. 10) write(6,*) 'calling cstates, expdes'
             call cstates(ar,nr,ai,ni,xest,dnest,lest,
     &                    states,nstatea,mostobs,nsaveset)
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
               do i = 1,nydata(j,k)
                 ydata(i,j,k) = yest(i,j,k1)
               enddo  ! i, nydata
            enddo

          enddo   !k, nsets

      endif

*     Initializations done ***********************************

       if (globcrit) then

        call exsensit(ar,nr,ai,ni,
     &             xest,nest,lest,
     &             ar(pxdisc),valdisc,ndiscpoint)

       if (optcrit(4:6).eq.'lsq') then

          expdesd = 0.0d0
          do jj   = 1,ndiscpoint
             expdesd = expdesd + valdisc(jj)
          enddo

       elseif (optcrit(4:6).eq.'opt') then
***** test only *****
          expdesd = 1.0d0
          do jj   = 1,ndiscpoint
             expdesd = expdesd*valdisc(jj)
          enddo
          expdesd = expdesd**0.25d0

       elseif (optcrit(4:5).eq.'r2') then
******** test only *******
         R2mean = gpar(7)
         R2dev  = gpar(8)
         expdesd = 1.0d0
         do jj   = 1,ndiscpoint
            x1 = -valdisc(jj)    ! the R2 value of test param
            x0 = R2mean
            sig= R2dev
            Desir = 1.0d0 - 1.0d0/(1.0d0 + dexp(-(x1-x0)/sig))
            expdesd = expdesd*Desir
         enddo
         expdesd = expdesd**0.25d0

       endif

       expdesd = -expdesd

      else !  local optcrits, with jcount > 1

        if(debug .gt. 10) write(*,*) 'calling hessup'

        call hessup(ar,nr,ai,ni,
     &               jtj,u,
     &               xest,nest,lest,
     &               newobs,nnew)

        do i = 1,nest
         do j = 1,nest
            jtj(i,j) = jtj0(i,j) + jtj(i,j)
         enddo
        enddo

*       Several 'alphabetic' optimal design criteria computed by
*       using the SVD of J'J.

*       Compute the SVD:

        ipath = 20
        tol   = 1d-16
        call dlsvrr(nest, nest, jtj, nest, ipath, info, tol, irank,
     &              s, u, nest, v, nest, jtj1,ar(pwrk2svd),ar(pwrksvd))


        if (irank.lt.nest .or. s(nest).lt.1.0d-30) then
            if (optmonit.gt.0) print *,'singular jtj, irank',irank
            expdesd = 1.0d+30
            return

        elseif (optcrit(1:1) .eq. 'd') then

*       minimize the generalized variance, i.e. maximize
*       logdet(jtj)

            logdet = 0.0d0
            do i = 1,nest
                logdet = logdet + log(s(i))
            enddo
            expdesd = - logdet

        elseif (optcrit(1:1) .eq. 'c') then

*       minimize the condition number

            expdesd = s(1)/s(nest)

        elseif (optcrit(1:1) .eq. 't') then

*       minimize the trace

            sum = 0.0d0
            do i = 1,nest
                sum = sum + s(i)
            enddo
            expdesd = sum

        elseif (optcrit(1:1) .eq. 'v') then

*       minimize the variance of a single parameter, j = expparam(1), i.e.
*       the diagonal (j,j) of inv(J' * J):

               j   = exppar(1)
               var = 0.0d0
               do i = 1,nest
                  var = var + u(j,i)**2/s(i)
               enddo
               expdesd = var

        elseif (optcrit(1:1) .eq. 's') then
*
*       minimize the generalized variance (D-optimum) of a subset
*       of the parameters given originally in the vector 'expparam'
*       (cf. the routine 'subset')

            do p1 = 1,nexpp
                do p2 = 1,nexpp
                    jtj2(p1,p2) = jtj(exppar(p1),exppar(p2))
                enddo
            enddo

*       compute first '-logdet(jtj)':

            logdet = 0.0d0
            do i = 1,nest
                logdet = logdet - log(s(i))
            enddo

*       Compute the SVD of jtj2:

            ipath = 20
            tol   = 1d-16
            call dlsvrr(nexpp,nexpp,jtj2,nexpp,ipath,info,tol,irank,
     &              s,u,nexpp,v,nexpp,jtj1,ar(pwrk2svd),ar(pwrksvd))

            if (irank.lt.nexpp.or.s(nexpp).lt.1.0d-30) then
               if (optmonit.eq.1) print *,'singular jtj2, irank',irank
               expdesd = 1.0d+30
               return
            endif

*       compute  '-logdet(jtj) + logdet(jtj2)'

            do i = 1,nexpp
                logdet = logdet + log(s(i))
            enddo
            expdesd = logdet

        else

            write(*,*) 'The Design criteria not implemented'
            stop

        endif

       endif  ! optcrit

      if (optmonit.gt.0) then
         if (mod(icount,optmonit).eq.0 .and.
     &        task(1:3).eq.'sen'    .and.
     &        objfun(1:3).eq.'exp')
     &        write(*,*) 'count and exp ', icount,expdesd

      endif

      return
      end

