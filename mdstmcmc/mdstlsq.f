C
CCCCC NOTE THIS IS NOT USED ANYMORE CCCCCC
C
C $Id: mdstlsq.f,v 1.3 2008/11/17 14:34:38 mjlaine Exp $
C ------------------------------------------------------------------------
C mdstmcmc library 
C File: mdstlsq.f
C Purpose: interface from mcmc code to modest lsq function
C
C Marko Laine <Marko.Laine@Helsinki.fi>
C ------------------------------------------------------------------------
C $Author: mjlaine $
C $Revision: 1.3 $
C $Date: 2008/11/17 14:34:38 $
C ------------------------------------------------------------------------
C $Log: mdstlsq.f,v $
C Revision 1.3  2008/11/17 14:34:38  mjlaine
C gfortran -Wall cleanups, array bound in mdstcov.F
C
C Revision 1.2  2007/09/11 19:23:36  mjlaine
C file not used
C
C Revision 1.1  2002/08/28 04:57:10  mjlaine
C Initial revision
C
C ------------------------------------------------------------------------
C
C F77 version
C uses dcopy (BLAS)
      function mdstlsq(ar,nr,ai,ni,xaux,naux)
      implicit none
      real*8 mdstlsq
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)
      integer*4      naux           
      real*8         xaux(naux+1) 
        
      real*8   lsq
      external lsq
      external dcopy

      real*8 ss

      include 'common3.inc'
        
      if (naux .ne. dnest) then
         write(*,*) 'mcmc-code thinks npar=',naux
         write(*,*) 'modest thinks npar=',dnest
         stop 'STOP in mdstlsq: illegal parameters'
      end if

c pitääkö tämä tehdä????
      call dcopy(dnest,xaux,1,ar(pxest),1)
        
      ss = lsq(ar,nr,ai,ni,ar(pxest),dnest)
      mdstlsq = ss

      return
      end
