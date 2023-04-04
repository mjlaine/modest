      subroutine lodeode(ar,nr,ai,ni,
     &     y,s0,states, nstatea, nstates,
     &     yest,mostydat,mostobs, nydata,
     &     nsaveset, xdata, nx, nsets, nobs, t0aux,
     &     gpar, ngpar, lpar, nlpar,
     &     rworkod,lrwrkod,iworkod,liwrkod,
     &     rtolod,lrtolod,atolod,latolod,
     &     nest,xaux)
*     CALLED BY

*     - SODE

*     SUBROUTINES NEEDED

*     - Fode   (ODE containing the model), DF, Jac (Optional)

*     PURPOSE

*     This subroutine calculates the states in the user given points (in x)

*     ARGUMENTS
      implicit none
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4 nstatea,nstates,mostydat,mostobs,nsaveset
      integer*4 nx,nsets,ngpar,nlpar
      integer*4 lrwrkod,liwrkod,lrtolod,latolod,nest
c
      real*8    y(nstatea,nest)         ! y is actually the whole ar
      real*8    s0(nstatea,nsets)
      real*8    states(nstatea,mostobs,nsaveset)
      real*8    yest(mostydat,mostobs,nsaveset)
      integer*4 nydata(mostobs,nsets)
      real*8    xdata(nx,mostobs,nsets)
      integer*4 nobs(nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)
      real*8    t0aux
      real*8    rworkod(lrwrkod),rtolod(lrtolod),atolod(latolod)
c     real*8    teta0(nest)
      real*8    xaux(nest)
      integer*4 iworkod(liwrkod)

*     LOCAL VARIABLES
      real*8    tjint,tj1,step,stepj,tstep ! interpolation for graphics output
      integer*4 j,k,k1,ii,nj,njsum ! loop indices

      real*8   swork(1000)

      logical dumpbin

*     ODESSA VARIABLES  (see odessa)
      integer*4    neq(2),iopt(3)
      real*8       t
      integer*4    istate

      real*8    tj              ! point at which solution
                                ! is desired
*     Externals in solvers

c     external odefode, DF, Jac
      external ddrfode, Jac

*     COMMON BLOCKS
      include 'common3.inc'

*     For Odessa input, the parameters saved in vector 'par'. If
*     sensitivities computed, the values of 'par' are changed
*     by Odessa, so save temporarily in 'xaux'.

      if (iset.eq.0) write(*,*) 'lodeode: iset',iset

      if (isoptod.eq.1) then
         call Xparser(ar,nr,ai,ni,'sel',xaux,nest,ai(plest))
c       call dcopy(nest,xaux,1,teta0(1),1)
      endif

      neq(1)  = nstates
      neq(2)  = dnest
      istate  = 1
      iopt(1) = ioptod
      iopt(2) = isoptod
      iopt(3) = idfod
      iworkod(1) = mlod
      iworkod(2) = muod
      iworkod(6) = mxstepod
      rworkod(1) = tcritod
      if (itolod.eq.1) then
         rtolod(1)=srtolod
         atolod(1)=satolod
      endif

      k  = iset
      k1 = mod(k,nsaveset)+1
      njsum = 0
      step  = xdata(1,nobs(k),k)-t0aux

      do j = 1,nobs(k)-1
         iobs = j

*         Set the end point for integration. 
         tj   = xdata(1,j+1,k)
         t    = t0aux

*         Integrate. Copy results to states.

         if(dumpcall .eq. 1 .and. dumpfile(1:6) .ne. 'nodump')
     &        then

            if (dumpfile(len_trim(dumpfile)-3:len_trim(dumpfile)) 
     &           == '.bin') then
               dumpbin = .true.
            else
               dumpbin = .false.
            end if

              
            if (j+1.lt.nobs(k)) then
               tj1   = xdata(1,j,k)
               stepj = tj - tj1
               nj    = max(aint(stepj/step*ndumpp),1.0d0)
               njsum = njsum + nj
               tstep = stepj/float(nj)
            else
               tj1   = xdata(1,j,k)
               stepj = tj - tj1
               nj    = max(ndumpp-njsum-1,1)
               tstep = stepj/float(nj)
            endif

            call xsetf(1) ! allow lsode messages

            tjint = xdata(1,j,k)+tstep
            do while(tjint. lt. tj-tj*1.0d-12)
               if(debug .gt. 10) write(*,*) ' calling lsodes, tjint'

ccc ioptd p.o. iopt !!
c                 call odessa(odefode,DF,neq,y,teta0,t0aux,tjint,itolod,
c     &                rtolod,atolod,itaskod,istate,ioptod,
c     &                rworkod,lrwrkod,iworkod,liwrkod,Jac,mfod)


c              write(*,*) itaskod, lrwrkod, liwrkod, mfod
c
c HARD CODED FOR NOW
               mfod = 222

               call dlsodes(ddrfode,neq(1),y,t0aux,tjint,itolod,
     &              rtolod(1),atolod,itaskod,istate,ioptod,
     &              rworkod,lrwrkod,iworkod,liwrkod,Jac,mfod)
               
               if(istate .le. 0) then
                  if (istate .lt. -2) then
                     write(*,*) 'Something is wrong in the LSODES'
                     write(*,*) 'istate = ', istate
                     stop
                  else
                     istate = 3    ! continue with istate = 3
                     call xsetf(0) ! no more messages
                  endif
               endif

*                Auxiliary states:
               if (dnsaux.gt.0) then
                  call ddrfode(nstates,tjint,y,swork)
c                 call odefode(nstates,tjint,y,teta0,swork)
               endif
               if(debug.gt.10) write(6,*)'calling Observations'
               call observations(y,nstates,
     &              yest(1,j,k1),nydata(j,k),
     &              xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &              s0(1,k),
     &              gpar(1),ngpar,
     &              lpar(1,k),nlpar,
     &              j,k)

               if (dumpbin) then
                  write(11) tjint
                  do ii = 1,nydata(j,k)
                     write(11) yest(ii,j,k1)
                  enddo
                  do ii = 1,nstatea
                     write(11) y(ii,1)
                  enddo
               else
                  write(11,*) tjint
                  do ii = 1,nydata(j,k)
                     write(11,*) yest(ii,j,k1)
                  enddo
                  do ii = 1,nstatea
                     write(11,*) y(ii,1)
                  enddo
               end if
c
               if(nx .gt. 1) then
*                 Each interval is a separate integration
*                 if there are more than one x variable.
c
                  istate = 1
                  t0aux=tjint
               endif
               tjint = tjint+tstep
            enddo
         endif

         if(debug .gt. 10) write(*, *) ' calling lsodes'
c          call odessa(odefode, DF, neq, y, teta0, t0aux, tj, itolod,
c     &              rtolod, atolod, itaskod, istate, ioptod,
c     &              rworkod, lrwrkod, iworkod, liwrkod, Jac, mfod)


         mfod = 222

c          write(*,*) itaskod, itolod, lrwrkod, liwrkod, mfod, ioptod
c          write(*,*) istate
c          write(*,*)  neq(1), y(1,1), y(70,1)
c          write(*,*)  'from - to ' , t0aux,tj
c          write(*,*) iworkod(1), iworkod(liwrkod)
         call dlsodes(ddrfode, neq, y, t0aux, tj, itolod,
     &        rtolod, atolod, itaskod, istate, ioptod,
     &        rworkod, lrwrkod, iworkod, liwrkod, Jac, mfod)
c
         if(istate .le. 0) then
            if (istate .lt. -2) then
               write(*,*) 'Something is wrong in the LSODES'
               write(*,*) 'istate = ', istate
               stop
            else
               istate = 3       ! continue with istate = 3
               call xsetf(0)    ! no more messages
            endif
         endif

*         Auxiliary states:
         if (dnsaux.gt.0) then
c           call odefode(nstates,tj,y,teta0,swork)
            call ddrfode(nstates,tj,y,swork)
         endif
         call dcopy(nstatea, y(1,1), 1, states(1,j+1,k1),1)
c          call dcopy(nstatea, y(1), 1, states(1,j+1,k1),1)

         if(dumpcall .eq. 1 .and. dumpfile(1:6) .ne. 'nodump')
     &        then
            if(debug.gt.10) write(6,*)'calling Observations'
            call observations(y,nstates,
     &           yest(1,j,k1),nydata(j,k),
     &           xdata(1,1,k),nx,nobs(k),
     &           dnsaux,nstatea,
     &           s0(1,k),
     &           gpar(1),ngpar,
     &           lpar(1,k),nlpar,
     &           j,k)


            if (dumpbin) then
               write(11) tj
               do ii = 1,nydata(j,k)
                  write(11) yest(ii,j,k1)
               enddo
               do ii = 1,nstatea
                  write(11) y(ii,1)
               enddo
            else
               write(11,*) tj
               do ii = 1,nydata(j,k)
                  write(11,*) yest(ii,j,k1)
               enddo
               do ii = 1,nstatea
                  write(11,*) y(ii,1)
               enddo
            end if

         endif

         if(nx .gt. 1) then
*             Each interval is a separate integration
*             if there are more than one x variable.
            istate = 1
            t0aux  = xdata(1,j+1,k)
         endif

      enddo                     ! j
c      
*     Parse the original values back
c     if(ioptod(2).eq.1)
c       call Xparser(ar,nr,ai,ni,'par',xaux,nest,ai(plest))
*************************
      return
      end
