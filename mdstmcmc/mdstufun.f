C $Id: mdstufun.f,v 1.3 2007/09/12 08:14:51 mjlaine Exp $
C ------------------------------------------------------------------------
C mdstmcmc library 
C File: mdstufun.f
C Purpose: Interface to user defined function usrfun
C 
C Marko Laine <Marko.Laine@Helsinki.fi>
C ------------------------------------------------------------------------
C $Author: mjlaine $
C $Revision: 1.3 $
C $Date: 2007/09/12 08:14:51 $
C ------------------------------------------------------------------------
C $Log: mdstufun.f,v $
C Revision 1.3  2007/09/12 08:14:51  mjlaine
C more header fixes
C
C Revision 1.2  2007/09/12 08:12:24  mjlaine
C headers added
C 
C ------------------------------------------------------------------------
C
C
      function mdstufun(ar,nr,ai,ni,xopt,nopt,len)

C     ARGUMENTS
      implicit none
      integer*4 len
      real*8 mdstufun(len)
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)

      integer*4      nopt       ! number of variables in Xopt
      real*8         xopt(nopt) ! variables to be optimized

c     real*8         mdstufund(len)
c     external       mdstufund

      interface
         function mdstufund(ar,nr,ai,ni,xopt,nopt,
     &        states,nstates,nsaux,nstatea,nsaveset,
     &        yest,
     &        xdata,nx,mostobs,nsets,
     &        ydata,mostydat,
     &        weight,
     &        nobs,
     &        nydata,
     &        states0,
     &        gpar,ngpar,
     &        lpar,nlpar,len)
         
         implicit none
         integer*4    len
         real*8       mdstufund(len)
         integer*4    nr,ni
         real*8       ar(nr)
         integer*4    ai(ni)

         integer*4    nopt,nsets,nx,mostobs,mostydat
         integer*4    nstates,nsaux,nstatea,nsaveset,ngpar,nlpar
         real*8       xopt(nopt) ! variables to be optimized
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
         end function mdstufund
      end interface


C     COMMON BLOCKS

      include 'common3.inc'

      mdstufun = mdstufund(ar,nr,ai,ni,xopt,nopt,
     &     ar(pstates),dnstates,dnsaux,dnstatea,dnsavese,
     &     ar(pyest),
     &     ar(pxdata),dnx,dmostobs,dnsets,
     &     ar(pydata),dmostyda,
     &     ar(pweight),
     &     ai(pnobs),
     &     ai(pnydata),
     &     ar(pstates0),
     &     ar(pgpar),max(dngpar,1),
     &     ar(plpar),max(dnlpar,1),len)
      return
      end

      function mdstufund(ar,nr,ai,ni,xopt,nopt,
     &     states,nstates,nsaux,nstatea,nsaveset,
     &     yest,
     &     xdata,nx,mostobs,nsets,
     &     ydata,mostydat,
     &     weight,
     &     nobs,
     &     nydata,
     &     states0,
     &     gpar,ngpar,
     &     lpar,nlpar,len)

      implicit none
      integer*4    len
      real*8 mdstufund(len)
      integer*4    nr,ni
      real*8       ar(nr)
      integer*4    ai(ni)

      integer*4    nopt,nsets,nx,mostobs,mostydat
      integer*4    nstates,nsaux,nstatea,nsaveset,ngpar,nlpar
      real*8       xopt(nopt)   ! variables to be optimized
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


C     local variables
      integer*4    j,k,k1
      integer*4    icount             ! count for calls
 

      interface
         function usrfun(len,s,ns,yest,ny, 
     &        xdata,nx,nobs, 
     &        nsaux,nstatea, 
     &        ydata,weight,mostydat, 
     &        states0, 
     &        gpar,ngpar, 
     &        lpar,nlpar, 
     &        iset)

         implicit none
         integer*4 len
         real*8 usrfun(len)
         integer*4 ns,nx,nobs,nsaux,nstatea,mostydat,ngpar,nlpar,iset
         integer*4 ny(nobs)
         real*8    s(nstatea,nobs)
         real*8    yest(mostydat,nobs)
         real*8    xdata(nx,nobs)
         real*8    ydata(mostydat,nobs)
         real*8    weight(mostydat,nobs)
         real*8    states0(nstatea)
         real*8    gpar(ngpar)
         real*8    lpar(nlpar)
         end function usrfun
      end interface

      include 'common3.inc'
      data  icount   /0/        ! initalization
      icount = icount + 1

C     Note: the loop below makes sense if nsets>1 and nsaveset=nsets.
C     Typically in task 'opt', however, nset = 1 = nsaveset.

      do k=1,nsets
         iset = k
         k1   = mod(k,nsaveset)+1
         iobs = -1              !flag: all observations

c        if(debug .gt. 10) write(6,*) 'calling cstates, opt'
c
c need to do this again
         call sstate(ar,nr,ai,ni,xopt,nopt,
     &        states,nstatea,mostobs,nsaveset)

C        Calculate the estimated observations Yest, given by the user given
C        observation function

         do j = 1,nobs(k)
            iobs = j
c
c            if(debug .gt. 10) write(*,*) ' calling Observations'
c

            call observations(states(1,j,k1),nstates,
     &           yest(1,j,k1),nydata(j,k),
     &           xdata(1,1,k),nx,nobs(k),nsaux,nstatea,
     &           states0(1,k),
     &           gpar(1),ngpar,
     &           lpar(1,k),nlpar,
     &           j,k)
c
c
         enddo                  ! j

c        if(debug .gt. 10) write(*,*) 'calling the user written optfun'

         mdstufund = usrfun(len,states(1,1,k1),nstates,
     &        yest(1,1,k1),nydata(1,k),
     &        xdata(1,1,k),nx,nobs(k),nsaux,nstatea,
     &        ydata,weight,mostydat,
     &        states0(1,k),
     &        gpar(1),ngpar,
     &        lpar(1,k),nlpar,
     &        k)


      enddo                     ! k

c      if (optmonit.gt.0) then
c         if (mod(icount,optmonit).eq.0 .and.
c     &        task(1:3).eq.'sen'    .and.
c     &        objfun(1:3).eq.'opt')
c     &        write(*,*) 'count and opt ', icount,optd
c      endif

      return
      end
