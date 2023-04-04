C $Id: mdstdump.f,v 1.2 2007/05/07 16:56:56 mjlaine Exp $
C ------------------------------------------------------------------------
C mdstmcmc library 
C File: mdstdump.f
C Purpose: Subroutine to dump yest to mat-file
C
C Marko Laine <Marko.Laine@Helsinki.fi>
C ------------------------------------------------------------------------
C $Author: mjlaine $
C $Revision: 1.2 $
C $Date: 2007/05/07 16:56:56 $
C ------------------------------------------------------------------------
C $Log: mdstdump.f,v $
C Revision 1.2  2007/05/07 16:56:56  mjlaine
C header comments added
C
C ------------------------------------------------------------------------
C
C
      subroutine mdstdump(ar,nr,ai,ni,xopt,nopt)

C     ARGUMENTS
      implicit none
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)


      integer*4      nopt       ! number of variables in Xopt
      real*8         xopt(nopt) ! variables to be optimized



      interface
         subroutine mdstdumpd(ar,nr,ai,ni,xopt,nopt,
     &        states,nstates,nsaux,nstatea,nsaveset,
     &        yest,
     &        xdata,nx,mostobs,nsets,
     &        ydata,mostydat,
     &        weight,
     &        nobs,
     &        nydata,
     &        states0,
     &        gpar,ngpar,
     &        lpar,nlpar)
         
         implicit none
         integer*4    nr,ni
         real*8       ar(nr)
         integer*4    ai(ni)

         integer*4      nopt    ! number of variables in Xopt
         real*8         xopt(nopt) ! variables to be optimized

         integer*4    nsets,nx,mostobs,mostydat
         integer*4    nstates,nsaux,nstatea,nsaveset,ngpar,nlpar
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
         end subroutine mdstdumpd
      end interface


C     COMMON BLOCKS

      include 'common3.inc'

      call mdstdumpd(ar,nr,ai,ni,xopt,nopt,
     &     ar(pstates),dnstates,dnsaux,dnstatea,dnsavese,
     &     ar(pyest),
     &     ar(pxdata),dnx,dmostobs,dnsets,
     &     ar(pydata),dmostyda,
     &     ar(pweight),
     &     ai(pnobs),
     &     ai(pnydata),
     &     ar(pstates0),
     &     ar(pgpar),max(dngpar,1),
     &     ar(plpar),max(dnlpar,1))
      return
      end

      subroutine mdstdumpd(ar,nr,ai,ni,xopt,nopt,
     &     states,nstates,nsaux,nstatea,nsaveset,
     &     yest,
     &     xdata,nx,mostobs,nsets,
     &     ydata,mostydat,
     &     weight,
     &     nobs,
     &     nydata,
     &     states0,
     &     gpar,ngpar,
     &     lpar,nlpar)

      implicit none
      integer*4    nr,ni
      real*8       ar(nr)
      integer*4    ai(ni)

      integer*4      nopt       ! number of variables in Xopt
      real*8         xopt(nopt) ! variables to be optimized


      integer*4    nsets,nx,mostobs,mostydat
      integer*4    nstates,nsaux,nstatea,nsaveset,ngpar,nlpar
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
      integer*4 nobsall

      interface
         subroutine dumpyest(yest,mostydat,nobs, iset, nobsall)
         implicit none
         integer*4 nobs,mostydat,nobsall,iset
         real*8    yest(mostydat,nobs)
         end subroutine dumpyest

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

c     calculate total number of observations
      nobsall=0
      do k=1,nsets
         nobsall = nobsall + nobs(k)
      enddo

c      write(*,*) 'Calling mdstdump'
      do k=1,nsets
         iset = k
         k1   = mod(k,nsaveset)+1

ccc just a call to states -> cstates with dumpcall = 1
         dumpcall = 1
cccccc need to call states and observations to get yest
         call sstate(ar,nr,ai,ni,xopt,nopt,
     &        states,nstatea,mostobs,nsaveset)
         dumpcall = 0

ccc with dumpcall = 1, no need to call observations or dumpyest
C$$$         do j = 1,nobs(k)
C$$$            call observations(states(1,j,k1),nstates,
C$$$     &           yest(1,j,k1),nydata(j,k),
C$$$     &           xdata(1,1,k),nx,nobs(k),nsaux,nstatea,
C$$$     &           states0(1,k),
C$$$     &           gpar(1),ngpar,
C$$$     &           lpar(1,k),nlpar,
C$$$     &           j,k)
C$$$         enddo
cccccc
c        dumpyest dumps yest to mat-file
C$$$         call dumpyest(yest(1,1,k1),mostydat,nobs(k),iset,nobsall)
      enddo

      return
      end
