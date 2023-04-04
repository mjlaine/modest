      subroutine eulode(s,s0,states,nstatea,nstates,
     &                  yest,mostydat, mostobs,nydata,
     &                  nsaveset, xdata, nx, nsets, nobs, t0aux,
     &                  gpar, ngpar, lpar, nlpar,
     &                  workeul)

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
c      real*8    workeul(nstatea)
      real*8    workeul(2*nstatea)

*     LOCAL VARIABLES
      real*8    tj,tjint,tj1,step,stepj,tstep     ! interpolation for graphics output
      integer*4 j,k,k1,ii,nj,njsum                ! loop indices

*     Externals in solvers

      external ddrFode

*     COMMON BLOCKS

      include 'common3.inc'

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

              tjint = t0aux+tstep
              do while(tjint. lt. tj-1.0d-12)

                  if(debug .gt. 10)
     &                write(*,*) ' calling deuler'
                  call deuler(nstates,t0aux,s,ddrFode,
     &                        tjint,stpeul,workeul)

*                 Auxiliary states:
                  if (dnsaux.gt.0) then
                     call ddrFode(nstates,tjint,s,workeul)
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

                  t0aux = tjint
                  tjint = tjint+tstep
              enddo
          endif                                                   !VMT 18.11.    93 end

          if(debug .gt. 10)
     &        write(*,*) ' calling deuler'
              call deuler(nstates,t0aux,s,ddrFode,
     &                    tj,stpeul,workeul)

*         Auxiliary states:
          if (dnsaux.gt.0) then
             call ddrFode(nstates,tj,s,workeul)
          endif

          call dcopy(nstatea,s,1,states(1,j+1,k1),1)

          if(dumpcall .eq. 1 .and. dumpfile(1:6) .ne. 'nodump')
     &    then
              if(debug.gt.10) write(6,*)'calling Observations'
              call observations(s,nstates,
     &                     yest(1,j,k1),nydata(j,k),
     &                     xdata(1,1,k),nx,nobs(k), dnsaux,nstatea,
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

          t0aux=xdata(1,j+1,k)

      enddo

      return
      end

