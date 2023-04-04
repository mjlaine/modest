      subroutine odefode(neq, t, y, par, dy)

c     an auxiliary routine between Odessa and Fode

      integer*4 neq
      real*8    t
      real*8    y(*)
      real*8    par(*)
      real*8    dy(*)

      include 'common3.inc'

      if (iset.eq.0) write(*,*) 'odefode: iset',iset

c      call Xparser('par',max(dnest,1),y(plest),par)
      call Fode1(neq, t, y, dy,
     &          y(pxdata),dnx,dmostobs,dnsets,dnsaux,
     &          y(prnobs),
     &          y(pstates0),
     &          y(pgpar),max(dngpar,1),
     &          y(plpar),max(dnlpar,1),
     &          iobs,iset)

      return
      end


