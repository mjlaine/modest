            subroutine Inits0(ns, t, s,
     &          xdata, nx, nobs,
     &          nsaux,nstatea,
     &          states0,
     &          gpar,ngpar,
     &          lpar,nlpar,
     &          iobs,iset )

*     test model for modest/lsq/sode
      implicit none

*     ARGUMENTS

      integer*4 ns,nsaux,nstatea  ! n of state variables
      real*8    t                 ! time
      real*8    s(nstatea)       ! state variables

      integer*4 nx,nobs,ngpar,nlpar
      real*8    xdata(nx,nobs)
      real*8    states0(nstatea)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar)
      integer*4 iobs,iset

*     LOCAL VARIABLES

      if (iobs.eq.iobs) then
        write(*,*) 'dummy routine for inits0'
        stop
      endif

      return
      end


