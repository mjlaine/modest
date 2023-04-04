      subroutine ddrfode(neq, t, y, dy)

c     an auxiliary routine between ddriv2 and Fode

      integer*4 neq, nneq
      real*8    t
      real*8    y(*)
      real*8    dy(*)

      include 'common3.inc'

      call Fode1(neq, t, y, dy,
     &          y(pxdata),dnx,dmostobs,dnsets,dnsaux,
     &          y(prnobs),
     &          y(pstates0),
     &          y(pgpar),max(dngpar,1),
     &          y(plpar),max(dnlpar,1),
     &          iobs,iset)

      return
      end

