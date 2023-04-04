C $Id: mdstcov.F,v 1.6 2009/02/10 11:38:11 mjlaine Exp $
C ------------------------------------------------------------------------
C mdstmcmc library 
C File: mdstcov.F
C Purpose: Calculates covariance matrix from modest
C
C          returns also the mse and nobs
C          F77 version
C          uses dcopy (BLAS), dpotrf, dpotri, dgesvd (LAPACK)
C
C Marko Laine <Marko.Laine@Helsinki.fi>
C ------------------------------------------------------------------------
C $Author: mjlaine $
C $Revision: 1.6 $
C $Date: 2009/02/10 11:38:11 $
C ------------------------------------------------------------------------
C $Log: mdstcov.F,v $
C Revision 1.6  2009/02/10 11:38:11  mjlaine
C precision checks on constants
C
C Revision 1.5  2008/11/17 14:34:38  mjlaine
C gfortran -Wall cleanups, array bound in mdstcov.F
C
C Revision 1.4  2008/05/05 21:07:01  mjlaine
C max(1,nobs-dnest) for mse calculations
C
C Revision 1.3  2007/09/11 19:30:24  mjlaine
C sstype & sstrans
C
C Revision 1.2  2007/05/07 16:56:56  mjlaine
C header comments added
C
C ------------------------------------------------------------------------
C
C 
C 
      subroutine mdstcov(ar,nr,ai,ni,cmat,npar,msevec,nobsvec,nycol,
     &     condmax,sstype,sstrans)
      implicit none
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)
      integer*4      npar    
      real*8         cmat(npar,npar) 
      integer*4      nycol
      real*8         msevec(nycol)
      integer*4      nobsvec(nycol)
      real*8         condmax
      integer*4      sstype
      real*8      sstrans

      real*8         jtj(npar,npar) 
      real*8         jtj1(npar,npar)
c      real*8         lest(4,npar)
c      integer*4      p1,p2 

      real*8    u(npar,npar),s(npar),work(npar*5)
      integer*4 lwork
      real*8 sum, tol
      integer*4 jj      
        
      real*8 lsq
      external dcopy, dpotrf, dgesvd, lsq

      real*8 sse
      integer*4 info, i, j

#include "common3.inc"

      real*8         mse
      integer*4      nobs
      real*8 sw(dmostyda)
      real*8 lsqvec(dmostyda)

      interface
         function mdstlsqs(ar,nr,ai,ni,xaux,naux,ny,sstype,sstrans)
         implicit none
         integer*4 ny
         real*8         mdstlsqs(ny)
         integer*4      nr,ni
         real*8         ar(nr)
         integer*4      ai(ni)
         integer*4      naux
         real*8         xaux(naux+1)
         integer*4      sstype
         real*8      sstrans
c        logical        dumpit
         end function mdstlsqs
      end interface

c     ! par ~ ar(pxest), npar=dnest
c     ! jac ~ ar(pjtj)  dnest*dnest

      if (dmostyda .ne. nycol) then
         write(*,*) 'error with dmostyda and nycol', dmostyda, nycol
         stop 
      endif

c calculate jtj with Hessian
      if (jacout .eq. 0) then
         cmat(1,1) = -1.0d0
         nobs = dntot
         return
c        stop 'jacout .eq. 0!'
      end if

      call hess(ar,nr,ai,ni,
     $     jtj,jtj1,
     $     ar(pxest),dnest,ai(plest))

c residual ss ja mse
      sse  = lsq(ar,nr,ai,ni,ar(pxest),dnest)
      mse  = sse/(dntot-dnest)
c     rmse  = dsqrt(mse)
      nobs = dntot ! OBS dntot is really nobs*nycol

      if (dnest .ne. npar) then
         stop 'mdstcov: illegel input'
      end if
c copy jtj to cmat
      call dcopy(dnest*dnest,jtj,1,cmat,1)

      lwork = 5*npar

!      open(101, file="cmat1.dat", status='replace')
!      do i=1,npar
!         write(101,*) cmat(i,:)
!      enddo
!      close(101)

c cmat contains now j'j
c calculate svd of symmetric matrix cmat = u*s*u'
c!!!  call dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
      call dgesvd('A','N',npar,npar,cmat,npar,s,u,npar,u,
     &     npar,work,lwork,info)
      if (s(1) == 0.0d0) then
         write(*,*) 'DGESVD: s(1)=0, can not go on'
         stop
      end if
      if (info > 0) then
         write(*,*) 'dgesvd info = ', info
      end if
      if (condmax>0.0d0) then
         tol = s(1)/condmax
         if (s(npar) <= tol) then
            write(*,*) 'NOTE: adjusting svd in mdstcov'
            do i=1,npar
               if (s(i) < tol) then
                  s(i) = tol
               end if
            end do
            info = 0
         end if
      elseif (s(npar) .eq. 0.0d0) then
         info = -1
      endif

      if (info .ne. 0) then
         cmat(1,1) = -1.0d0
         write(*,*) 'Modest cmat singular'
         return
      end if


c invert u*s*u' and save it in cmat
      do i = 1,npar
         call dscal(npar,1.0d0/dsqrt(s(i)),u(1,i),1)
      enddo
      do i = 1,npar
         do j = 1,i
            sum = 0.0d0
            do jj = 1,npar
               sum = sum + u(i,jj)*u(j,jj)
            enddo
            cmat(i,j) = sum
            cmat(j,i) = sum
         enddo
      enddo
c
c cmat is now inv(jtj)
c scale with mse
      do i=1,npar
         do j=i,npar
            cmat(i,j) = cmat(i,j) * mse
            if (i /= j) cmat(j,i) = cmat(i,j)
         end do
      end do

c call sumofweights to get corrected nobs
      call sumofweights(ar(pweight),dmostyda,dmostobs,dnsets,
     &     ai(pnobs),ai(pnydata),sw)

c      write(*,*) 'old nobs, mse, most', nobs, mse, dmostyda
c fix nobs mse
c we now really return msevec and nobsvec so this code should be cleaned
      nobs = sw(1)
      mse = sse/max(1.0d0,dble(nobs-dnest))
      lsqvec =  mdstlsqs(ar,nr,ai,ni,ar(pxest),dnest,dmostyda,
     &     sstype,sstrans)
      do i=1,dmostyda
         nobsvec(i) = sw(i)
         msevec(i) = lsqvec(i)/max(1.0d0,dble(nobsvec(i)-dnest))
c         write(*,*) 'mse: ', i,msevec(i)
      enddo
c      write(*,*) 'new nobs, mse', nobs, mse

      return
      end
c
c subroutine to calculate sum of weights
c
      subroutine sumofweights(weight,mostydat,mostobs,nsets,
     &     nobs,nydata,sw)

      implicit none

      include 'common3.inc'

      integer*4 mostydat,mostobs,nsets
      real*8    weight(mostydat,mostobs,nsets)
      integer*4    nobs(nsets)
      integer*4    nydata(mostobs,nsets)

      real*8 sw(dmostyda)
      integer*4 i, j1, j, k

      do i=1,dmostyda !   was dnest
         sw(i) = 0.0d0
      enddo

      if (model(1:3) .eq. 'ode' .and. obss0.eq.0) then
         j1  = 2
      else
         j1  = 1
      endif

      do k=1,nsets
         do j=j1,nobs(k)
            do i=1,nydata(j,k)
               if (usew.eq.1) then
                  sw(i)  = sw(i) + weight(i,j,k)
               else
                  sw(i)  = sw(i) + 1.0d0
               endif
            enddo
         enddo
      enddo

      end
