      subroutine fode(ns,t,s,ds,
     &           xdata,nx,nobs,
     &           nsaux,nstatea,
     &           states0,
     &           gpar,ngpar,
     &           lpar,nlpar,
     &           iobs,iset)
 
      implicit none
 
c     arguments
 
       integer*4 ns,nsaux,nstatea    !n of state variables
       real*8    t           !time
       real*8    s(nstatea) !state variables
       real*8    ds(ns)      !derivatives
       integer*4 nx,nobs,ngpar,nlpar
       real*8    xdata(nx,nobs)
       real*8    states0(nstatea)
       real*8    gpar(ngpar)
       real*8    lpar(nlpar)
       integer*4 iobs,iset
 
c     local variables (all user defined variables must be declared here!)
      integer*4 i
      real*8 A,B,R,z,k1,k2

      include 'boxo.inc'
 
c     user code:

      A = s(1)
      B = s(2)

      e1 = e1*1.0d+6
      e2 = e2*1.0d+6

      R  =  8.314
      z  = 1.0d0/R * ( 1.0d0/Temp - 1.0d0/Tmean)
      k1 = k1mean * dexp(-e1*z)
      k2 = k2mean * dexp(-e2*z)

      ds(1) = - k1 * A
      ds(2) =   k1 * A - k2 * B

      do i = 3,ns        !just to test large ns
         ds(i) = 0.0d0
      enddo

c      s(ns+1) = ds(1)    !just to demostrate auxiliary states
c      s(ns+2) = ds(2)
c      s(ns+3) = t
 
      return
      end
