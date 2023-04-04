      subroutine simflex(ar,dar,ai,dai,
     &                   f,x0,m,a,itmax,tol1,tol2,bounds,ibound,
     &                   y,s,is,work,iwork,ierror,optmonit)

*    This program passes the parameters to simflex_exec that does
*    the optimization

*      implicit none

*     ARGUMENTS

      integer*4 dar,dai
      real*8    ar(dar)
      integer*4 ai(dai)

      integer*4 m              ! number of arguments
      integer*4 itmax          ! maximum number of iterations
      integer*4 is             ! row dimension of s
      integer*4 iwork(2*m+2)   ! A work vector
      integer*4 ierror         ! returns the error conditions, 0: tol1 and
                               ! tol2 satisfied, 1: tol1 satisfied, 2: tol2
                               ! satisfied, 3: itmax exceeded, ON INPUT:
                               ! >0: printing, <0: no printing
      integer*4 optmonit       ! monitoring at each optmonit:th iteration
      integer*4 ibound(m)      ! types of the bounds
      real*8 f                 ! The objective function ot the form f(x,nx)
*                                Now f(ar,dar,ai,dai,x,nx)
      real*8 x0(m)             ! The initial values for the arguments
      real*8 a                 ! The relative size of the simplex
      real*8 tol1              ! A variance tolerance for the values of f
      real*8 tol2              ! A relative var. tolerance for the values of f
      real*8 y                 ! The final value of f
      real*8 s(is,m+1)         ! The simplex
      real*8 work(5*m+2)       ! A work vector
      real*8 bounds(2,m)       ! bounds for constraints
      external f


*     the initial values are transformed to unbounded ones

      call constraints(x0,m,ibound,bounds,-1)

      call simflex_exec(ar,dar,ai,dai,
     &    f,x0,m,a,itmax,tol1,tol2,bounds,ibound,
     &    y,s,is,work(1),work(m+1),work(2*m+1),work(3*m+1),
     &    iwork(1),iwork(m+2),ierror,optmonit)

      return
      end


      subroutine simflex_exec(ar,dar,ai,dai,
     &           f,x0,m,a,itmax,tol1,tol2,bounds,
     &           ibound,y,s,is,xreflection,xexpansion,
     &           xcontraction,xcentre,iwork,ni,ierror,optmonit)

      implicit none

*     ARGUMENTS

      integer*4 dar,dai
      real*8    ar(dar)
      integer*4 ai(dai)

      integer*4 m              ! number of arguments
      integer*4 itmax          ! maximum number of iterations
      integer*4 is             ! row dimension of s
      integer*4 iwork(2*m+2)   ! A work vector
      integer*4 ierror         ! returns the error conditions, 0: tol1 and
                               ! tol2 satisfied, 1: tol1 satisfied, 2: tol2
                               ! satisfied, 3: itmax exceeded, ON INPUT:
                               ! >0: printing, <0: no printing
      integer*4 optmonit
      integer*4 ibound(m)      ! types of the bounds
      integer*4 start, iterate ! adresses
      integer*4 iterations     ! n of iterations
      real*8 f                 ! The objective function ot the form f(x,nx)
*                                Now f(ar,dar,ai,dai,x,nx)
      real*8 x0(m)             ! The initial values for the arguments
      real*8 a                 ! The relative size of the simplex
      real*8 tol1              ! A variance tolerance for the values of f
      real*8 tol2              ! A relative var. tolerance for the values of f
      real*8 y                 ! The final value of f
      real*8 s(is,m+1)         ! The simplex
      real*8 bounds(2,m)       ! bounds for constraints
      real*8 xcentre(m)        !
      real*8 xreflection(m)    !
      real*8 xexpansion(m)     !
      real*8 xcontraction(m)   !
      real*8 reflection        !
      real*8 expansion         !
      real*8 contraction       !
      real*8 yreflection       !
      real*8 yexpansion        !
      real*8 ycontraction      !
      real*8 ymean             !
      real*8 ysquare           !
      real*8 test1             !
      real*8 test2             !
      real*8 p, p1, t, r       ! auxiliary variables

      integer*4 i,j,kpos(1),nk !iopt(4) ! auxiliary indices

      external f

      parameter (reflection=1.d0,expansion=2.d0,contraction=.5d0)

      integer   ni(m+1)

      assign 1 to start
      assign 2 to iterate

*     maaritetaan aloitussimplex


1     t=dsqrt(2.d0)   !  start
      r=m

      p1=a*(dsqrt(r+1.d0)+r-1.d0)/(r*t)
      p =a*(dsqrt(r+1.d0)-1.d0)/(r*t)

      do j=1,m
          s(j,1)=x0(j)   !  koordinaatit sarakkeittain
      enddo

      do i=2,m+1
          do j=1,m

          if(x0(j) .ne. 0.d0) then
              s(j,i)=x0(j)+p*x0(j)
          else
              s(j,i)=x0(j)+p
          endif

          enddo

          if(x0(i-1) .ne. 0.d0) then
              s(i-1,i)=x0(i-1)+p1*x0(i-1)
          else
              s(i-1,i)=x0(i-1)+p1
          endif

      enddo

*     sijoitetaan funktion f arvot simpleksin s sarmissa
*     s:n riville m+1

      if(itmax.eq.0) then

          call constraints(x0,m,ibound,bounds,1)
          y=f(ar,dar,ai,dai,x0,m)
          return

      endif

      do i=1,m+1
          call dcopy(m,s(1,i),1,x0,1)
          call constraints(x0,m,ibound,bounds,1)
          s(m+1,i)=f(ar,dar,ai,dai,x0,m)
      enddo

*     jarjestetaan s:n sarakkeet rivin m+1 mukaan kasvavaan jarjestykseen

      kpos(1)=m+1
      nk=1
      iterations=0

2     iterations=iterations+1    ! iterate

*     ylla xcentre on tyovektorina

      call dscolr(ar,dar,ai,dai,
     &            m+1,m+1,s,is,0,0,0,nk,kpos,iwork,m+1,ni)

*     konvergenssitestit

      ymean=0.d0
      ysquare=0.d0

      do i=1,m+1

          ymean=ymean+s(m+1,i)
          ysquare=ysquare+s(m+1,i)**2

      enddo

      ymean=ymean/(m+1)

      if(ierror .lt. 0) then
          if (optmonit.gt.0) then
              if (mod(iterations,optmonit).eq.0)
     &        write(6,*) iterations,s(m+1,m+1),'   ',s(m+1,1)
          endif
          call dcopy(m,s(1,m+1),1,x0,1)
          call constraints(x0,m,ibound,bounds,1)
          call dcopy(m,s(1,1),1,x0,1)
          call constraints(x0,m,ibound,bounds,1)
          if (optmonit.gt.0) then
           if (mod(iterations,optmonit).eq.0) then
             write(6,*) (x0(j),j=1,m)
             write(6,*)
           endif
          endif
      endif

      if((ysquare-(m+1)*ymean**2)/(m+1).le.0.d0) then
          test1=0.d0
      else
          test1=dsqrt((ysquare-(m+1)*ymean**2)/(m+1))
      endif

      if (ymean.eq.0.0d0) then
         write(*,*) 'zero response in Simflex'
         stop
      endif

      test2=dabs(test1/ymean)

      if(test1.lt.tol1.or.test2.lt.tol2) then

          ierror=0

          if(test2.ge.tol2) then

              ierror=1

          else if(test1.ge.tol1) then

              ierror=2

          endif

          call dcopy(m,s(1,1),1,x0,1)
          call constraints(x0,m,ibound,bounds,1)
          y=f(ar,dar,ai,dai,x0,m)

          return

      endif

      if(iterations.ge.itmax) then

          ierror=3
          call dcopy(m,s(1,1),1,x0,1)
          call constraints(x0,m,ibound,bounds,1)
          y=f(ar,dar,ai,dai,x0,m)

          return

      endif

*     lasketaan keskipiste ilman huonointa (so. m+1:nnetta)

      do j=1,m
          xcentre(j)=0.d0
      enddo

      do j=1,m
          do i=1,m

              xcentre(j)=xcentre(j)+s(j,i)

          enddo
          xcentre(j)=xcentre(j)/m

      enddo

*     peilaus; huom! huonoin <-> s(.,m+1) ja paras <-> s(.,1)

      do j=1,m

          xreflection(j)=(1.d0+reflection)*xcentre(j)
     &                  -reflection*s(j,m+1)

      enddo

      call dcopy(m,xreflection,1,x0,1)
      call constraints(x0,m,ibound,bounds,1)
      yreflection=f(ar,dar,ai,dai,x0,m)

      if(yreflection.le.s(m+1,1)) then  ! peilaus paras

*     jatko (=expansion)

          do j=1,m

              xexpansion(j)=expansion*xreflection(j)
     &                     +(1.d0-expansion)*xcentre(j)

          enddo

          call dcopy(m,xexpansion,1,x0,1)
          call constraints(x0,m,ibound,bounds,1)
          yexpansion=f(ar,dar,ai,dai,x0,m)

          if(yexpansion.le.yreflection) then  ! jatko paras

              do j=1,m

                  s(j,m+1)=xexpansion(j)

              enddo

              s(m+1,m+1)=yexpansion

          else   ! peilaus parempi kuin jatko

              do j=1,m

                  s(j,m+1)=xreflection(j)

              enddo

              s(m+1,m+1)=yreflection

          endif

      else if(yreflection.ge.s(m+1,m)) then  ! peilaus vah. 2. huonoin

*     valitaan parempi kahdesta huonoimmasta

          if(yreflection.lt.s(m+1,m+1)) then  ! peilaus 2. huonoin

              do j=1,m

                  s(j,m+1)=xreflection(j)

              enddo

              s(m+1,m+1)=yreflection

          endif

*         kutistus (=contraction)

          do j=1,m

              xcontraction(j)=contraction*s(j,m+1)
     &                       +(1.d0-contraction)*xcentre(j)

          enddo

          call dcopy(m,xcontraction,1,x0,1)
          call constraints(x0,m,ibound,bounds,1)
          ycontraction=f(ar,dar,ai,dai,x0,m)

          if(ycontraction.gt.s(m+1,m+1)) then ! kutistus huonoin

*         tayskutistus (= total contraction)

              do i=2,m+1
                  do j=1,m

                      s(j,i)=.5d0*(s(j,i)+s(j,1))

                  enddo

                  call dcopy(m,s(1,i),1,x0,1)
                  call constraints(x0,m,ibound,bounds,1)
                  s(m+1,i)=f(ar,dar,ai,dai,x0,m)!kutistetun simpleksin y:t

              enddo

          else          !  kutistus parempi kuin huonoin

              do j=1,m

                  s(j,m+1)=xcontraction(j)

              enddo

              s(m+1,m+1)=ycontraction

          endif

      else      !  peilaus parhaan ja 2. huonoimman valissa

          do j=1,m

              s(j,m+1)=xreflection(j)

          enddo

          s(m+1,m+1)=yreflection

      endif
      goto iterate
      end
