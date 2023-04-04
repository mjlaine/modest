      subroutine Lmdiff (ar,nr,ai,ni,
     &                   f,
     &                   xest,nest,lest,
     &                   ipiv,z,
     &                   bounds,ibound,
     &                   nobs,nsets,
     &                   states,states0,nstatea,mostobs,nsaveset,
     &                   xdata,nx,
     &                   ydata,weight,
     &                   yest,xjac,jtj,jtj1,xaux,
     &                   mostydat,nydata,
     &                   gpar,ngpar,
     &                   lpar,nlpar)


       implicit none

*     ARGUMENTS

      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4 nest           ! number of arguments
      real*8    f              ! The lsq objective function of the form
                               ! f(ar,nr,ai,ni,xest,nest)
      real*8    xest(nest)     ! The initial values for the arguments
      integer*4 lest(4,nest)

      integer*4 ipiv(nest)     ! pivot         for dgefa/dgesl
      real*8    z(nest,3)      ! 1 col: 'z',   right hand for "
                               ! 2 col: 'z' -> 'h' (overwri by gesl)
                               ! 3 col: 'xtest',unbounded
      real*8    bounds(2,nest)    ! bounds for constraints
      integer*4 ibound(nest)      ! types of the bounds
      integer*4 nsets,ngpar,nlpar
      integer*4 nstatea,mostobs,nsaveset,mostydat,nx
      integer*4 nobs(nsets)
      real*8    states(nstatea,mostobs,nsaveset)
      real*8    states0(nstatea,nsets)
      real*8    xdata(nx,mostobs,nsets)
      real*8    ydata(mostydat,mostobs,nsets)
      real*8    weight(mostydat,mostobs,nsets)
      real*8    yest(mostydat,mostobs,nsaveset)
      real*8    xjac(mostydat,nest,mostobs,nsaveset)
      real*8    jtj(nest,nest),jtj1(nest,nest),xaux(nest)
      integer*4 nydata(mostobs,nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)


*     local variables
      real*8    na(4),b,dnrm2,norm0,norm,tol
      integer*4 p,p1,p2,it,nn,k,i,j,k1,info,imax,imin
      integer*4 idmax,idmin
c      real*8    xtest(4),zz(4,3)

      external f,dnrm2,idmax,idmin

      include 'common3.inc'

      call dcopy(nest,xest(1),1,xaux(1),1)         !xest:original
      call constraints(xaux,nest,ibound,bounds,-1) !xaux:unbounded

      norm0 = f(ar,nr,ai,ni,xest,nest)
      na(1) = norm0
      it    = 0

*     Default itmax
      if (itmaxlm.eq.0) then
          itmaxlm = iterdflm
          nn      = 0
      endif

5     continue
       it = it+1

       if (optmonit.gt.0) then
          if (mod(it,optmonit).eq.0) then
             write(*,*) 'LevMar: count, norm and arguments  ', it,norm0
             write(*,*) (xest(i),i=1,nest)
          endif
       endif

       do p = 1,nest
          call dcopy(nest,0.0d0,0,jtj(1,p),1)
          z(p,1) = 0.0d0
       enddo

       do k = 1,nsets
         iset = k
         k1   = mod(k,nsaveset) + 1
         iobs = -1                    ! flag: states for all obss

            call jacoblm(ar,nr,ai,ni,
     &                 xest,nest,lest,
     &                 xjac,jtj1,xaux,mostydat,
     &                 xdata,nx,mostobs,nsets,
     &                 weight,bounds,ibound,
     &                 nobs,nydata,
     &                 states,nstatea,nsaveset,
     &                 yest,states0,
     &                 gpar,ngpar,
     &                 lpar,nlpar)

*            update J'* (yest-ydata)

             do p = 1,nest
               do j = 1,nobs(k)
                 do i = 1,nydata(j,k)
                   z(p,1) = z(p,1) +
     &                      xjac(i,p,j,k1)*(-yest(i,j,k1)+ydata(i,j,k))
                 enddo
               enddo
             enddo

*            update the approx. Hessian

              do p1 = 1, nest
                do p2 = 1, nest
                   jtj(p1,p2) = jtj(p1,p2) + jtj1(p1,p2)
                enddo
              enddo

       enddo  ! k

*     Default alfa
      if (it.eq.1 .and. alfalm.eq.0.0d0) then
          alfalm = 0.01d0 * dnrm2(nest,jtj(1,1),nest+1)
      endif

       do p1 = 1,nest
          jtj1(p1,p1) = jtj(p1,p1)
          jtj (p1,p1) = jtj(p1,p1) + alfalm
       enddo

       call dcopy(nest,z(1,1),1,z(1,2),1)            !orig. z
       call dgefa(jtj,nest,nest,ipiv,info)
       call dgesl(jtj,nest,nest,ipiv,z(1,2),0)  !h

       do p=1,nest
          z(p,3) = xaux(p) + z(p,2)
       enddo

       call dcopy(nest,z(1,3),1,xest(1),1)
       call constraints(xest,nest,ibound,bounds,1)
       norm = f(ar,nr,ai,ni,xest,nest)

       do while (norm.gt.norm0) 
          alfalm = 10.0d0*alfalm
          do p = 1,nest
             jtj(p,p) = jtj1(p,p) + alfalm
          enddo
          call dcopy(nest,z(1,1),1,z(1,2),1)            !orig. z
          call dgefa(jtj,nest,nest,ipiv,info)
          call dgesl(jtj,nest,nest,ipiv,z(1,2),0)
          do p=1,nest
             z(p,3) = xaux(p) + z(p,2)
          enddo
          call dcopy(nest,z(1,3),1,xest(1),1)
          call constraints(xest,nest,ibound,bounds,1)
          norm = f(ar,nr,ai,ni,xest,nest)
       enddo

       call dcopy(nest,z(1,3),1,xaux(1),1)
       alfalm = 0.3d0*alfalm
       norm0  = norm
       if (nn.eq.0) then  ! stop by tol
           do i = 1,3
              na(i) = na(i+1)
           enddo
           na(4) = norm0
           imax  = idmax(4,na,1)
           imin  = idmin(4,na,1)
           b     = na(imax) - na(imin)
           tol   = b/(norm0 + 1d-15)
           if (tol.lt.reltollm) then
              itmaxlm = it-1
              it      = it+1
           endif
       endif !nn

      if (it.lt.itmaxlm) goto 5

      return
      end

