      subroutine hessup(ar,nr,ai,ni,
     &                  jtj,wrkjtj,
     &                  xest,nest,lest,
     &                  newobs,nnew)

*     calculates the update part of the 'Hessian' J'*J for expdes

*     ARGUMENTS
      implicit none
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      Integer*4      nest,nnew
      real*8         jtj(nest,nest)  ! the update part of J'J
      real*8         wrkjtj(nest,nest) ! work area,the volatile part of J'J
      real*8         xest(nest)
      Integer*4      lest(4,nest)
      integer*4      newobs(nnew,2)

*     LOCAL VARIABLES

      integer*4      j,k,k1,p1,p2,ii   ! loop indices
      integer*4      prevk             ! previous value of k

*     PARAMETERS & COMMONS

      include 'common3.inc'

*     THE ALGORITHM

      if(debug .gt. 11) write(6,*) ' in hessup'

      call dset(nest*nest,0.0d0,jtj(1,1),1) ! initialize to 0

      prevk = -1

      do ii = 1,nnew

*        compute the update jtj piecewise from wrkjtj

         k = newobs(ii,2)
         j = newobs(ii,1)

         if (model(1:3).eq.'ode') then

*          always the whole set computed;
*          check whether already done

           if (k .ne. prevk) then
              prevk = k
              iset  = k
              iobs  = -1
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

*             update the update part for approx. Hessian
              do p1 = 1, nest
                 do p2 = 1, nest
                     jtj(p1,p2) = jtj(p1,p2) + wrkjtj(p1,p2)
                 enddo
              enddo

           endif

         elseif (model(1:3).eq.'alg'.or.model(1:3).eq.'imp') then

*          alg or imp model, observations one by one
           iset = k
           iobs = j
           call jacob(ar,nr,ai,ni,
     &          xest,nest,lest,
     &          ar(pxjac),wrkjtj,ar(pxaux),dmostyda,
     &          ar(pxdata),dnx,dmostobs,dnsets,
     &          ar(pweight),
     &          ai(pnobs),ai(pnydata),
     &          ar(pstates),dnstatea,dnsavese,
     &          ar(pyest),ar(pstates0),
     &          ar(pgpar),max(dngpar,1),
     &          ar(plpar),max(dnlpar,1))

*           update the update part for approx. Hessian
            do p1 = 1, nest
               do p2 = 1, nest
                   jtj(p1,p2) = jtj(p1,p2) + wrkjtj(p1,p2)
               enddo
            enddo

         endif    ! model

      enddo  ! ii

      return
      end

