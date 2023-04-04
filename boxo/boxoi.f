      subroutine inits0(ns,t,s,
     &           xdata,nx,nobs,
     &           nsaux,nstatea,
     &           states0,
     &           gpar,ngpar,
     &           lpar,nlpar,
     &           iobs,iset)
 
      implicit none
 
c     arguments
 
       integer*4 ns,nsaux,nstatea   !n of state variables
       real*8    t           !time
       real*8    s(nstatea) !state variables
 
       integer*4 nx,nobs,ngpar,nlpar
       real*8    xdata(nx,nobs)
       real*8    states0(ns)
       real*8    gpar(ngpar)
       real*8    lpar(nlpar)
       integer*4 iobs,iset
 
* local variables (all user defined variables must be declared here!)
 
       integer i
       include 'boxo.inc'
 
* user code:


       do i = 1,ns
         s(i) = states0(i)
       enddo

       t = xdata(1,1)
 
       return
       end
