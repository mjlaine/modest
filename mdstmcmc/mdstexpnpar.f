C F77 version
C copy modest dnexp to naux
      subroutine mdstexpnpar(ar,nr,ai,ni,naux)
      implicit none
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)
      integer*4      naux           
      real*8         xaux(naux+1) 
        

      include 'common3.inc'

      naux = dnexp

      return
      end
