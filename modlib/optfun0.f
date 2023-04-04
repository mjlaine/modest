      real*8 function optfun(s,ns,yest,ny,
     &                 xdata,nx,nobs,
     &                 nsaux,nstatea,
     &                 ydata,weight,mostydat,
     &                 states0,
     &                 gpar,ngpar,
     &                 lpar,nlpar,
     &                 iset)

*     test model for desirability function
      implicit none

*     ARGUMENTS

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

****** Template, an example *********************
*     LOCAL VARIABLES (All user defined arguments must be declared here !)
c      integer   i,ii,j,last
c      real*8    x1,x2,x3,x4,x5,x0,sig,Desir(10),s(6,20)
c
c      last   = nobs
c      x1     = s(3,last)/s(1,1)
c      x0     = 0.75
c      sig    = 0.05
c      Desir(1) = 1.0d0/(1.0d0 + dexp(-(x1-x0)/sig))
c
c      optfun = 1.0d0
c      do i = 1,1
c           optfun = optfun*Desir(i)
c      enddo
c      if (optfun.gt.0.0d0) optfun = -optfun**1.0d0
c
      optfun = 0.0d0
      if (1.eq.1) then
         write(*,*) 'dummy optfun'
         stop
      endif
      return
      end


