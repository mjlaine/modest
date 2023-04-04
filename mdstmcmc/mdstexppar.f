C F77 version
C uses dcopy (BLAS)
C copy modest exp parameters from ar(pxexp) to xaux
      subroutine mdstexppar(ar,nr,ai,ni,xaux,naux,xl,xu)
      implicit none
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)
      integer*4      naux
      real*8         xaux(naux+1), xl(naux), xu(naux)
        
      real*8   lsq
      external lsq
      external dcopy
      integer i

      include 'common3.inc'
        
      if (naux .ne. dnexp) then
         stop 'mdstexppar: illegal parameters'
      end if

      call dcopy(dnexp,ar(pxexp),1,xaux,1)
      do i=1,naux
         xl(i) = ar(pbounds+(i-1)*2)
         xu(i) = ar(pbounds+(i-1)*2+1)
      end do

      return
      end
