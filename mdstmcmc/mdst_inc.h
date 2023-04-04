!!! $Id: mdst_inc.h,v 1.1 2011/10/28 07:23:52 mjlaine Exp $ -*- mode: F90; -*-
!!! ------------------------------------------------------------------------
!!! mdstmcmc library
!!! File: mdst_inc.f90
!!! Purpose: interface functions to get information on the current model 
!!!          from modest
!!!
!!! Marko Laine <Marko.Laine@Helsinki.fi>
!!! ------------------------------------------------------------------------
!!! $Author: mjlaine $
!!! $Revision: 1.1 $
!!! $Date: 2011/10/28 07:23:52 $
!!! ------------------------------------------------------------------------
!!!
!!!

!DEC$ IF .NOT. DEFINED(_mdst_inc)
!DEC$ DEFINE _mdst_inc

interface

   ! modestin lsq
   function lsq(ar,nr,ai,ni,xaux,naux)
     real*8         lsq
     integer*4      nr,ni
     real*8         ar(nr)
     integer*4      ai(ni)
     integer*4      naux           ! number of variables in Xopt
     real*8         xaux(naux+1)   ! variables to be optimized
   end function lsq

   function mdstlsq(ar,nr,ai,ni,xaux,naux)
     real*8         mdstlsq
     integer*4      nr,ni
     real*8         ar(nr)
     integer*4      ai(ni)
     integer*4      naux
     real*8         xaux(naux+1)
   end function mdstlsq

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
     real*8         sstrans
!    logical        dumpit
   end function mdstlsqs

   function mdstbnd(ar,nr,ai,ni,par,npar)
     logical         mdstbnd
     integer(4)      nr,ni
     real(8)         ar(nr)
     integer(4)      ai(ni)
     integer(4)      npar
     real(8)         par(npar)
   end function mdstbnd

   subroutine mdstcov(ar,nr,ai,ni,cmat,npar,mse,nobs,nycol,condmax, &
        sstype,sstrans)
     integer*4      nr,ni
     real*8         ar(nr)
     integer*4      ai(ni)
     integer*4      npar
     real*8         cmat(npar,npar)
     integer*4      nycol
     real*8         mse(nycol)
     integer*4      nobs(nycol)
     real*8      condmax
     integer*4   sstype
     real*8   sstrans
   end subroutine mdstcov

   subroutine mdstpar(ar,nr,ai,ni,xaux,naux)
     implicit none
     integer*4      nr,ni
     real*8         ar(nr)
     integer*4      ai(ni)
     integer*4      naux
     real*8         xaux(naux+1)
   end subroutine mdstpar

   subroutine mdstnpar(ar,nr,ai,ni,naux, nycol)
     implicit none
     integer*4      nr,ni
     real*8         ar(nr)
     integer*4      ai(ni)
     integer*4      naux, nycol
!     real*8         xaux(naux+1)
   end subroutine mdstnpar

   function mdstufun(ar,nr,ai,ni,xopt,nopt,len)
      implicit none
      integer*4      len
      real*8         mdstufun(len)
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)
      integer*4      nopt
      real*8         xopt(nopt)
   end function mdstufun

   subroutine mdstdump(ar,nr,ai,ni,xopt,nopt)
     implicit none
     integer*4      nr,ni
     real*8         ar(nr)
     integer*4      ai(ni)
      integer*4      nopt
      real*8         xopt(nopt)
   end subroutine mdstdump

   function mdstexptarget(ar,nr,ai,ni,xaux,naux)
     implicit none
     real*8         mdstexptarget
     integer*4      nr,ni
     real*8         ar(nr)
     integer*4      ai(ni)
     integer*4      naux
     real*8         xaux(naux+1)
   end function mdstexptarget

   subroutine mdstexpnpar(ar,nr,ai,ni,naux)
     implicit none
     integer*4      nr,ni
     real*8         ar(nr)
     integer*4      ai(ni)
     integer*4      naux           
     real*8         xaux(naux+1) 
   end subroutine mdstexpnpar

   subroutine mdstexppar(ar,nr,ai,ni,xaux,naux,xl,xu)
     implicit none
     integer*4      nr,ni
     real*8         ar(nr)
     integer*4      ai(ni)
     integer*4      naux
     real*8         xaux(naux+1), xl(naux), xu(naux)
   end subroutine mdstexppar

end interface

!DEC$ ENDIF
