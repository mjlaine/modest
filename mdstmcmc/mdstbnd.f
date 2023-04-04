C $Id: mdstbnd.f,v 1.3 2008/11/17 14:34:38 mjlaine Exp $
C ------------------------------------------------------------------------
C mdstmcmc library 
C File: mdstbnd.f
C Purpose: Check parameter bounds
C          returns .false. if par is out of bounds
C
C Marko Laine <Marko.Laine@Helsinki.fi>
C ------------------------------------------------------------------------
C $Author: mjlaine $
C $Revision: 1.3 $
C $Date: 2008/11/17 14:34:38 $
C ------------------------------------------------------------------------
C $Log: mdstbnd.f,v $
C Revision 1.3  2008/11/17 14:34:38  mjlaine
C gfortran -Wall cleanups, array bound in mdstcov.F
C
C Revision 1.2  2007/05/07 16:56:56  mjlaine
C header comments added
C
C ------------------------------------------------------------------------
C
C mdstbnd
      logical function mdstbnd(ar,nr,ai,ni,par,npar)
      implicit none
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)
      integer*4 npar
      real*8 par(npar)

      logical inside  
        
      integer i,j,k

      include 'common3.inc'

      integer*4 ibound(dnest)
      real*8 bounds(2,dnest)

c       ar(pbounds)
      call dcopy(dnest*2,ar(pbounds),1,bounds,1)
        
c       ai(pbtype)
      call icopy(dnest,ai(pbtype),1,ibound,1)
        
      inside = .true.
        
      do i=1,npar

         if (ibound(i) == -1) then
        
            if (par(i) .le. bounds(1,i)) then
               inside = .false.
            end if

         else if (ibound(i) == 1) then
            if (par(i) .ge. bounds(2,i)) then
               inside = .false.
            end if

         else if (ibound(i) == 2) then

            if (par(i) .le. bounds(1,i)) then
               inside=.false.
            elseif (par(i) .ge. bounds(2,i)) then
               inside=.false.
            end if
         end if
      end do

      mdstbnd = inside

      return
      end
