      subroutine  Falg(ns,s,
     &          xdata, nx, nobs,
     &          gpar,ngpar,
     &          lpar,nlpar,
     &          iobs,iset )

*     test model for modest/lsq/salg
      implicit none

*     ARGUMENTS

      integer*4 ns                ! n of state variables
      real*8    s(ns)             ! state variables

      integer*4 nx,nobs,ngpar,nlpar
      real*8    xdata(nx,nobs)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar)
      integer*4 iobs,iset

*     LOCAL VARIABLES (All user defined arguments must be declares here !!!)


*     dummy model for modest/ib

*     ARGUMENTS

      write(6,*) 'dummy falg for modest.olb'
      if (1.eq.1) stop

      return
      end

