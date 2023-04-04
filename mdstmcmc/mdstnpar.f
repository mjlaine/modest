C $Id: mdstnpar.f,v 1.3 2008/11/17 14:34:38 mjlaine Exp $
C ------------------------------------------------------------------------
C mdstmcmc library 
C File: mdstnpar.f
C Purpose: copy modest dnest to naux
C
C Marko Laine <Marko.Laine@Helsinki.fi>
C ------------------------------------------------------------------------
C $Author: mjlaine $
C $Revision: 1.3 $
C $Date: 2008/11/17 14:34:38 $
C ------------------------------------------------------------------------
C $Log: mdstnpar.f,v $
C Revision 1.3  2008/11/17 14:34:38  mjlaine
C gfortran -Wall cleanups, array bound in mdstcov.F
C
C Revision 1.2  2007/09/12 08:14:51  mjlaine
C more header fixes
C
C Revision 1.1.1.1  2003/02/28 11:24:20  mjlaine
C MCMC f90 code import
C
C ------------------------------------------------------------------------
C
C F77 version
C 
      subroutine mdstnpar(ar,nr,ai,ni,naux,nycol)
      implicit none
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)
      integer*4      naux, nycol
c     real*8         xaux(naux+1) 
        
      real*8   lsq
      external lsq

      include 'common3.inc'

      naux = dnest
c nycol is the n:o of columns in the first observation of the first batch
c but this is assumed to be same over the whole data
      nycol = ai(pnydata)

      return
      end
