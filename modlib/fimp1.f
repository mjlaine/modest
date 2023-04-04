      subroutine Fimp1(ns, s, f, iflg,
     &           xdata,nx,mostobs,nsets,nsaux,
     &           rnobs,
     &           gpar,ngpar,
     &           lpar,nlpar,
     &           iobs,iset)

*    Another aux routine: rnobs --> nobs(iset),
*    and to get xdata & lpar without data set dimension

*     ARGUMENTS

      integer*4 ns,nsaux          ! n of y variables (states or equations)
      real*8    s(ns+nsaux)       ! y variables (states or equations)
      real*8    f(ns)
      integer*4 iflg,nx,mostobs,ngpar,nlpar,nsets
      real*8    rnobs(nsets)
      real*8    xdata(nx,mostobs,nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)
      integer*4 iobs,iset
      integer*4 nstatea

      nstatea = ns + nsaux
      call  Fimp(ns, s, f, iflg,
     &           xdata(1,1,iset),nx,int(rnobs(iset)),
     &           nsaux,nstatea,
     &           gpar,ngpar,
     &           lpar(1,iset),nlpar,
     &           iobs,iset)

      return
      end
