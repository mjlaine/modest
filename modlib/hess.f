      subroutine hess(ar,nr,ai,ni,
     &                jtj,wrkjtj,
     &                xest,nest,lest)

*     calculates the 'Hessian' J'*J for EXPDES

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      Integer*4      nest
      real*8         jtj(nest,nest)   ! the approx Hessian
      real*8         wrkjtj(nest,nest)  ! work area,the volatile part of jtj
      Integer*4      lest(4,nest)
      real*8         xest(nest)

*     LOCAL VARIABLES

      integer*4      k,p1,p2      ! loop indices

*     COMMONS

      include 'common3.inc'

*     THE ALGORITHM

*     initialize 'jtj'

      if(debug .gt. 11) write(6,*) ' in hess'

      call dset(nest*nest,0.0d0,jtj(1,1),1) ! initialize to 0

      do k = 1,dnsets
         iset = k
         iobs = -1                    ! flag: states for all obss
         call jacob(ar,nr,ai,ni,
     &             xest,nest,lest,
     &             ar(pxjac),wrkjtj,ar(pxaux),dmostyda,
     &             ar(pxdata),dnx,dmostobs,dnsets,
     &             ar(pweight),
     &             ai(pnobs),ai(pnydata),
     &             ar(pstates),dnstatea,dnsavese,
     &             ar(pyest),ar(pstates0),
     &             ar(pgpar),max(dngpar,1),
     &             ar(plpar),max(dnlpar,1))

*        update the approx. Hessian

              do p1 = 1, nest
                do p2 = 1, nest
                   jtj(p1,p2) = jtj(p1,p2) + wrkjtj(p1,p2)
                enddo
              enddo


      enddo

      return
      end
