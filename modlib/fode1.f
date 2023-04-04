      subroutine Fode1(ns, t, s, ds,
     &           xdata,nx,mostobs,nsets,nsaux,
     &           rnobs,
     &           states0,
     &           gpar,ngpar,
     &           lpar,nlpar,
     &           iobs,iset)

*    Another aux routine: to get xdata & lpar without data set dimension

*     ARGUMENTS

      integer*4 ns,nsaux          ! n of y variables (states or equations)
      real*8    t                 ! time
      real*8    s(ns+nsaux)       ! y variables (states or equations)
      real*8    ds(ns)            ! derivatives
      integer*4 nx,mostobs,ngpar,nlpar,nsets
      real*8    rnobs(nsets)
      real*8    xdata(nx,mostobs,nsets)
      real*8    states0(ns+nsaux,nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)
      integer*4 iobs,iset

*     local variable
      integer*4 nstatea

      nstatea = ns + nsaux
       ! if (nstatea.ne.1457 .or.iset.ne.1) then
       !   write(*,*) 'fode1:iset,ns,nsaux,nstatea',iset,ns,nsaux,nstatea
       ! endif

      call Fode(ns, t, s, ds,
     &          xdata(1,1,iset),nx,int(rnobs(iset)),
     &          nsaux,nstatea,
     &          states0(1,iset),
     &          gpar,ngpar,
     &          lpar(1,iset),nlpar,
     &          iobs,iset)

      return
      end
