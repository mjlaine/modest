      subroutine observations(s,ns,yest,ny,
     &           xdata,nx,nobs,
     &           nsaux,nstatea,
     &           states0,
     &           gpar,ngpar,
     &           lpar,nlpar,
     &           iobs,iset)
 
      implicit none
 
c     arguments
 
       integer*4 ns,nsaux,nstatea    !n of state variables
       real*8    s(nstatea) !state variables
       integer*4 ny          !n of obs. vars
       real*8    yest(ny)    !observed variables
 
       integer*4 nx,nobs,ngpar,nlpar
       real*8    xdata(nx,nobs)
       real*8    states0(nstatea)
       real*8    gpar(ngpar)
       real*8    lpar(nlpar)
       integer*4 iobs,iset
 
* local variables (all user defined variables must be declared here!)

*   template:
c         do i=1,ny
c             yest(i) = s(i)
c         enddo
       if (1.eq.1) then
          write(*,*) 'dummy observat function'
          stop
       endif

       return
       end

