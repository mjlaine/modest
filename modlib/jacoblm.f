      subroutine jacoblm(ar,nr,ai,ni,
     &                 xest,nest,lest,
     &                 xjac,jtj,xaux,mostydat,
     &                 xdata,nx,mostobs,nsets,
     &                 weight,bounds,ibound,
     &                 nobs,nydata,
     &                 states,nstatea,nsaveset,
     &                 yest,states0,
     &                 gpar,ngpar,
     &                 lpar,nlpar)

*     PURPOSE
*     calculates (part of) the jacobian J of the response 'yest', and
*     (the respective part of) J'J.

*     ARGUMENTS

      implicit none
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4 nest,nsets,nx,mostobs,mostydat
      integer*4 nstatea,nsaveset,ngpar,nlpar
      Integer*4 lest(4,nest)
      real*8    xest(nest),xaux(nest)
      real*8    xjac(mostydat,nest,mostobs,nsaveset)! jacobian
      real*8    jtj(nest,nest)                      ! approx. hessian
      real*8    xdata(nx,mostobs,nsets)
      real*8    weight(mostydat,mostobs,nsets)
      integer*4 ibound(nest)
      real*8    bounds(2,nest)
      integer*4 nobs(nsets)
      integer*4 nydata(mostobs,nsets)
      real*8    states(nstatea,mostobs,nsaveset)
      real*8    yest(mostydat,mostobs,nsaveset)
      real*8    states0(nstatea,nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)


*     LOCAL VARIABLES

      real*8    dx,aux               ! auxiliary variables
      integer*4 i,j,k,k1,p,p1,p2,kk  ! loop indices

*     COMMON BLOCKS

      include 'common3.inc'

*     THE ALGORITHM

       call dset(nest*nest,0.0d0,jtj(1,1),1) ! initialize to 0
       k    = iset
       k1   = mod(k,nsaveset) + 1

         do p = 1,nest
            dx = dstep*xest(p)

            do kk = 1,2

               xest(p) = xest(p) + (-2*kk+3)*dx
               call cstates(ar,nr,ai,ni,xest,nest,lest,
     &                      states,nstatea,mostobs,nsaveset)

               do j = 1,nobs(k)
                  call observations(states(1,j,k1),nstatea-dnsaux,
     &                        yest(1,j,k1),nydata(j,k),
     &                        xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &                        states0(1,k),
     &                        gpar(1),ngpar,
     &                        lpar(1,k),nlpar,
     &                        j,k)

                 do i = 1,nydata(j,k)
                   if (kk.eq.1) then
                      xjac(i,p,j,k1) = yest(i,j,k1)/dx
                   else
                      aux=(bounds(2,p)-bounds(1,p))/2.0d0*dcos(xaux(p))
                      xjac(i,p,j,k1) = (xjac(i,p,j,k1)-yest(i,j,k1)/dx)*
     &                                    weight(i,j,k)*aux
                   endif
                 enddo ! i

               enddo   ! j
            enddo      ! kk
         enddo         ! p

         do p1 = 1,nest
            do p2 = p1,nest
               do j = 1,nobs(k)
                  do i = 1,nydata(j,k)
                     jtj(p1,p2) = jtj(p1,p2) +
     &                           xjac(i,p1,j,k1) * xjac(i,p2,j,k1)
                  enddo   ! i
               enddo      ! j
            enddo         ! p2
         enddo            ! p1

       do p1 = 1,nest
         do p2 = p1,nest
            jtj(p2,p1) = jtj(p1,p2)
         enddo
       enddo

       if(debug .gt. 10) write(6,*) 'calling Xparser'
       call xparser(ar,nr,ai,ni,'par',xest,nest,lest) ! restore params

      return
      end

