      subroutine Fode(ns,t, s, ds,
     &          xdata, nx, nobs,
     &          nsaux,nstatea,
     &          states0,
     &          gpar,ngpar,
     &          lpar,nlpar,
     &          iobs,iset )

*     test model for modest/lsq/sode
      implicit none

*     ARGUMENTS

      integer*4 ns,nsaux,nstatea  ! n of state variables (states & aux var)
      real*8    t                 ! time
      real*8    s(nstatea)        ! state variables
      real*8    ds(ns)            ! derivatives
      integer*4 nx,nobs,ngpar,nlpar
      real*8    xdata(nx,nobs)
      real*8    states0(nstatea)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar)
      integer*4 iobs,iset

*     LOCAL VARIABLES (All user defined arguments must be declared here !)

      if (iobs.eq.iobs) then
        write(*,*) 'dummy routine for Fode'
        stop
      endif

      return
      end



