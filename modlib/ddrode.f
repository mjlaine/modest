      subroutine ddrode(s,s0,states, nstatea, nstates,
     &                  yest,mostydat,mostobs, nydata,
     &                  nsaveset, xdata, nx, nsets, nobs, t0aux,
     &                  gpar, ngpar, lpar, nlpar,
     &                  workddr,lrwddr, iwrkddr,liwddr)

*     CALLED BY

*     - SODE

*     SUBROUTINES NEEDED

*     - Fode   (ODE containing the model)

*     PURPOSE

*     This subroutine calculates the states in the user given points (in x)

      implicit none
*     ARGUMENTS

      integer*4 nstatea,nstates,mostydat,mostobs,nsaveset
      integer*4 nx,nsets,ngpar,nlpar
      integer*4 lrwddr,liwddr
      real*8    s(nstatea)                       ! s is actually the whole ar/VMT 24.11.93
      real*8    s0(nstatea,nsets)
      real*8    states(nstatea,mostobs,nsaveset)
      real*8    yest(mostydat,mostobs,nsaveset)
      integer*4 nydata(mostobs,nsets)
      real*8    xdata(nx,mostobs,nsets)
      integer*4 nobs(nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)
      real*8    t0aux
      real*8    workddr(lrwddr)
      integer*4 iwrkddr(liwddr)


*     LOCAL VARIABLES
      real*8    tjint,tj1,step,stepj,tstep        ! interpolation for graphics output
      real*8    dummyr
      integer*4 j,k,k1,ii,nj,njsum                ! loop indices
      integer*4 iroot

      real*8 swork(1000)

*     DDRIV2 VARIABLES

      real*8    tj                                ! point at which solution
                                                  ! is desired
      integer*4 mstate                            ! status of integration
      real*8    ddrRoot

*     Externals in solvers

      external ddrFode, ddrRoot                      ! ddriv2,

*     COMMON BLOCKS

      include 'common3.inc'

*     Ddriv2 either integrates past given end point and
*     interpolates the result or adjusts its step size
*     to reach given end point.

      if (interpdr2 .eq. 1) then
          mstate = 1
      else
          mstate = -1
      end if

      k  = iset
      k1 = mod(k,nsaveset)+1
      njsum = 0
      step  = xdata(1,nobs(k),k)-t0aux

      do j = 1,nobs(k)-1

*         Variable iobs tells the current observation and may
*         be used in inits0 and Fode.

          iobs = j

*         Set the end point for integration.

          tj = xdata(1,j+1,k)

*         Integrate. Copy results to states.

          if(dumpcall .eq. 1 .and. dumpfile(1:6) .ne. 'nodump')
     &    then
              
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

              tjint = xdata(1,j,k)+tstep
              do while(tjint. lt. tj-tj*1.0d-12)
                  if(debug .gt. 10)
     &                write(*,*) ' calling ddriv2, tjint'
                  call ddriv2(nstates,t0aux,s,ddrFode,
     &               tjint,mstate, nrootdr2,epsdr2,ewtdr2,
     &               mintdr2,workddr,dlrwddr,
     &               iwrkddr,dliwddr,ddrRoot)
c          if (abs(mstate).eq.7) then
          if (abs(mstate).eq.5) then
c             write(*,*) 'a root found at ', t0aux,tjint
c             do ii=j,nobs(k)
c               xdata(2,ii,k)=601.0d0
c             enddo
             mstate = 1
             dummyr =  ddrRoot(nstates,t0aux,s,-1)
          endif
*                 Auxiliary states:
                  if (dnsaux.gt.0) then
                     call ddrFode(nstates,tjint,s,swork)
                  endif
                  if(debug.gt.10) write(6,*)'calling Observations'
                  call observations(s,nstates,
     &                     yest(1,j,k1),nydata(j,k),
     &                     xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &                     s0(1,k),
     &                     gpar(1),ngpar,
     &                     lpar(1,k),nlpar,
     &                     j,k)
                  write(11,*) tjint
                  do ii = 1,nydata(j,k)
                      write(11,*) yest(ii,j,k1)
                  enddo
                  do ii = 1,nstatea
                      write(11,*) s(ii)
                  enddo

                 if(nx .gt. 1) then

*                 Each interval is a separate integration
*                 if there are more than one x variable.

              if (abs(mstate).ne.5) t0aux=tjint
                  if (interpdr2 .eq. 1) then
                    mstate = 1
                  else
                    mstate = -1
                  end if

                  endif
                  tjint = tjint+tstep
              enddo
          endif

          if(debug .gt. 10)
     &        write(*,*) ' calling ddriv2, tj'
          call ddriv2(nstates, t0aux, s, ddrFode, tj,
     &       mstate, nrootdr2, epsdr2, ewtdr2,
     &       mintdr2, workddr, dlrwddr,
     &       iwrkddr, dliwddr, ddrRoot)

**********
          if (abs(mstate).eq.5) then
             write(*,*) 'a root found at ', t0aux,tj
             dummyr =  ddrRoot(nstates,t0aux,s,-1)
c             do ii=j,nobs(k)
c               xdata(2,ii,k)=601.0d0
c             enddo
             mstate = 1
          endif
*         Auxiliary states:
          if (dnsaux.gt.0) then
             call ddrFode(nstates,tj,s,swork)
          endif

          call dcopy(nstatea,s,1,states(1,j+1,k1),1)

          if(dumpcall .eq. 1 .and. dumpfile(1:6) .ne. 'nodump')
     &    then
              if(debug.gt.10) write(6,*)'calling Observations'
              call observations(s,nstates,
     &                     yest(1,j,k1),nydata(j,k),
     &                     xdata(1,1,k),nx,nobs(k),dnsaux,nstatea,
     &                     s0(1,k),
     &                     gpar(1),ngpar,
     &                     lpar(1,k),nlpar,
     &                     j,k)
              write(11,*) tj
              do ii = 1,nydata(j,k)
                  write(11,*) yest(ii,j,k1)
              enddo
              do ii = 1,nstatea
                   write(11,*) s(ii)
              enddo

           endif

          if(nx .gt. 1) then

*             Each interval is a separate integration
*             if there are more than one x variable.

              if (abs(mstate).ne.5) t0aux=xdata(1,j+1,k)
              if (interpdr2 .eq. 1) then
                  mstate = 1
              else
                  mstate = -1
              end if
              

          endif
      enddo

      return
      end


