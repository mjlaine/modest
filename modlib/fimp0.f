      subroutine Fimp(ns, s, f, iflg,
     &           xdata,nx,nobs,
     &           nsaux,nstatea,
     &           gpar,ngpar,
     &           lpar,nlpar,
     &           iobs,iset)

*    The user function for implicit algebraic systems

*     ARGUMENTS

      integer*4 ns,nsaux,nstatea   ! n of y variables (states or equations)
      real*8    s(nstatea)         ! y variables (states or equations)
      real*8    f(ns)
      integer*4 iflg,nx,nobs,ngpar,nlpar
      real*8    xdata(nx,nobs)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar)
      integer*4 iobs,iset

      if (1.eq.1) then
        write(*,*) 'dummy routine for Fimp'
        stop
      endif

      return
      end


