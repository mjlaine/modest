      subroutine jacob(ar,nr,ai,ni,
     &                 xaux,naux,laux,
     &                 xjac,jtj,xorg,mostydat,
     &                 xdata,nx,mostobs,nsets,
     &                 weight,
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

      integer*4 naux,nsets,nx,mostobs,mostydat
      integer*4 nstatea,nsaveset,ngpar,nlpar
      Integer*4 laux(4,naux)
      real*8    xaux(naux)
      real*8    xjac(mostydat,naux,mostobs,nsaveset)! jacobian
      real*8    jtj(naux,naux)                      ! approx. hessian
      real*8    xorg(naux)                          ! orig xaux saved
      real*8    xdata(nx,mostobs,nsets)
      real*8    weight(mostydat,mostobs,nsets)
      integer*4 nobs(nsets)
      integer*4 nydata(mostobs,nsets)
      real*8    states(nstatea,mostobs,nsaveset)
      real*8    yest(mostydat,mostobs,nsaveset)
      real*8    states0(nstatea,nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)


*     LOCAL VARIABLES

      real*8    dx                   ! an auxiliary variable
      integer*4 i,j,k,k1,p,p1,p2,kk  ! loop indices

*     COMMON BLOCKS

      include 'common3.inc'

*     THE ALGORITHM


      call dcopy(naux,xaux,1,xorg(1),1)     ! keep the original xaux
      call dset(naux*naux,0.0d0,jtj(1,1),1) ! initialize to 0

*     calculation of the jacobian

      if (iobs .le. 0) then

         do p=1,naux                      ! derivative for each parameter
            dx=dstep*xorg(p)              ! xorg contains the original parameters
            do kk=1,2                     ! central differences

               call dcopy(naux,xorg(1),1,xaux,1)
               xaux(p)=xaux(p)+(2.0*kk-3.0)*dx

               k    = iset
               k1   = mod(k,nsaveset) + 1

               if(debug .gt. 10) write(6,*) 'calling cstates, jacob'
               call cstates(ar,nr,ai,ni,xaux,naux,laux,
     &                      states,nstatea,mostobs,nsaveset)


                  do j=1,nobs(k)

                     if(debug.gt.10) write(6,*)'calling Observations'
                     call observations(states(1,j,k1),nstatea-dnsaux,
     &                        yest(1,j,k1),nydata(j,k),
     &                        xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &                        states0(1,k),
     &                        gpar(1),ngpar,
     &                        lpar(1,k),nlpar,
     &                        j,k)

                     do i=1,nydata(j,k)
                        if (kk.eq.1) then
                           xjac(i,p,j,k1)=(kk-1.5d0)/dx*yest(i,j,k1)
                        else
                           xjac(i,p,j,k1)=(xjac(i,p,j,k1)+(kk-1.5d0)/
     &                         dx*yest(i,j,k1))*weight(i,j,k)
                        endif

                     enddo                ! i
                  enddo                   ! j
            enddo                         ! kk
         enddo                            ! p

         do p1 = 1,naux
            do p2 = p1,naux
               do j = 1,nobs(k)
                  do i = 1,nydata(j,k)
                     jtj(p1,p2) = jtj(p1,p2) +
     &                           xjac(i,p1,j,k1) * xjac(i,p2,j,k1)
                  enddo   ! i
               enddo      ! j
            enddo         ! p2
         enddo            ! p1


      else                       ! iobs .gt. 0

         j = iobs

         do p=1,naux                      ! derivative for each parameter
            dx=dstep*xorg(p)              ! xorg contains the original parameters
            do kk=1,2                     ! central differences

               call dcopy(naux,xorg(1),1,xaux,1)
               xaux(p)=xaux(p)+(2.0*kk-3.0)*dx

               k    = iset
               k1   = mod(k,nsaveset) + 1

               if(debug .gt. 10) write(6,*) 'calling cstates, jacob'
               call cstates(ar,nr,ai,ni,xaux,naux,laux,
     &                      states,nstatea,mostobs,nsaveset)


              if(debug.gt.10) write(6,*)'calling Observations'
              call observations(states(1,j,k1),nstatea-dnsaux,
     &                        yest(1,j,k1),nydata(j,k),
     &                        xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &                        states0(1,k),
     &                        gpar(1),ngpar,
     &                        lpar(1,k),nlpar,
     &                        j,k)

              do i=1,nydata(j,k)
                 if (kk.eq.1) then
                    xjac(i,p,j,k1)=(kk-1.5d0)/dx*yest(i,j,k1)
                 else
                    xjac(i,p,j,k1)=(xjac(i,p,j,k1)+(kk-1.5d0)/
     &              dx*yest(i,j,k1))*weight(i,j,k)
                 endif

              enddo   !i
            enddo     !kk
         enddo        !p

         do p1 = 1,naux
            do p2 = p1,naux
                  do i = 1,nydata(j,k)
                     jtj(p1,p2) = jtj(p1,p2) +
     &                           xjac(i,p1,j,k1) * xjac(i,p2,j,k1)
                  enddo   ! i
            enddo         ! p2
         enddo            ! p1

      endif    ! iobs

      do p1 = 1,naux
         do p2 = p1,naux
            jtj(p2,p1) = jtj(p1,p2)
         enddo
      enddo

         call dcopy(naux,xorg(1),1,xaux,1)              ! restore xaux
         if(debug .gt. 10) write(6,*) 'calling Xparser'
         call xparser(ar,nr,ai,ni,'par',xaux,naux,laux) ! restore params

      return
      end

