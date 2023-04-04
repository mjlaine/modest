C $Id: mdstexptarget.f,v 1.2 2008/11/17 14:34:38 mjlaine Exp $
C ------------------------------------------------------------------------
C mdstmcmc library 
C File: mdstexptarget.f
C Purpose: Interface between modlib expdes and mcmclib
C          Returns the target function for annealing.
C
C Marko Laine <Marko.Laine@Helsinki.fi>
C ------------------------------------------------------------------------
C $Author: mjlaine $
C $Revision: 1.2 $
C $Date: 2008/11/17 14:34:38 $
C ------------------------------------------------------------------------
C $Log: mdstexptarget.f,v $
C Revision 1.2  2008/11/17 14:34:38  mjlaine
C gfortran -Wall cleanups, array bound in mdstcov.F
C
C Revision 1.1.1.1  2002/08/28 05:00:26  mjlaine
C MCMC f90 code import
C
C Revision 1.2  2002/08/28 05:00:16  mjlaine
C *** empty log message ***
C
C Revision 1.1  2002/05/17 11:12:11  mjlaine
C Initial revision
C
C ------------------------------------------------------------------------
C
C F77 version
C uses dcopy (BLAS)
      function mdstexptarget(ar,nr,ai,ni,xaux,naux)
      implicit none
      real*8 mdstexptarget
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)
      integer*4      naux           
      real*8         xaux(naux+1) 
        
      real*8   expdes
      external expdes
      external dcopy

      real*8 ee

      include 'common3.inc'
        
      if (naux .ne. dnexp) then
         write(*,*) 'mcmc-code thinks nxpar=',naux
         write(*,*) 'modest thinks nxpar=',dnexp
         stop 'STOP in mdstexptarget: illegal parameters'
      end if

c pitääkö tämä tehdä????
c     call dcopy(naux,xaux,1,ar(pxexp),1)
c     ee = expdes(ar,nr,ai,ni,ar(pxexp),dnexp)
      ee = expdes(ar,nr,ai,ni,xaux,naux)
      mdstexptarget = -ee

      return
      end
