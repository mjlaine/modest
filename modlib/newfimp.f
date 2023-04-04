      subroutine newfimp(neq, s, f, iflg)

c     an auxiliary routine between newrap and Fimp

      integer*4 neq,iflg
      real*8    s(*)
      real*8    f(*)

      include 'common3.inc'

      call Fimp1(neq, s, f, iflg,
     &          s(pxdata),max(dnx,1),dmostobs,dnsets,dnsaux,
     &          s(prnobs),
     &          s(pgpar),max(dngpar,1),
     &          s(plpar),max(dnlpar,1),
     &          iobs,iset)

      return
      end

