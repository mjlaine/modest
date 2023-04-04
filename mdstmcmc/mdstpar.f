C $Id: mdstpar.f,v 1.4 2008/11/17 14:34:38 mjlaine Exp $
C ------------------------------------------------------------------------
C mdstmcmc library 
C File: mdstpar.f
C Purpose: copy modest parameters from ar(pxest) to xaux
C   returns the current (optimized) parameter vector from modest to xaux
C   F77 version, uses dcopy (BLAS)
C 
C Marko Laine <Marko.Laine@Helsinki.fi>
C ------------------------------------------------------------------------
C $Author: mjlaine $
C $Revision: 1.4 $
C $Date: 2008/11/17 14:34:38 $
C ------------------------------------------------------------------------
C $Log: mdstpar.f,v $
C Revision 1.4  2008/11/17 14:34:38  mjlaine
C gfortran -Wall cleanups, array bound in mdstcov.F
C
C Revision 1.3  2007/09/12 08:14:51  mjlaine
C more header fixes
C
C Revision 1.2  2007/09/12 08:12:24  mjlaine
C headers added
C
C ------------------------------------------------------------------------
C
      subroutine mdstpar(ar,nr,ai,ni,xaux,naux)
      implicit none
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)
      integer*4      naux           
      real*8         xaux(naux+1) 
        
      external dcopy

      include 'common3.inc'

      ! par ~ ar(pxest), npar=dnest
      ! jac ~ ar(pjtj) dnest,dnest
        
      if (naux .ne. dnest) then
         stop 'mdstpar: illegal parameters'
      end if

      call dcopy(dnest,ar(pxest),1,xaux,1)

      return
      end
