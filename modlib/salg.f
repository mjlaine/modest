      subroutine Salg(ar,nr,ai,ni,
     &                states,nstatea,mostobs,nsaveset,
     &                xdata, nx, nsets,
     &                nobs,
     &                gpar,ngpar,
     &                lpar,nlpar)

*     CALLED BY

*     - the chosen objective function (e.g. LSQ)

*     SUBROUTINES NEEDED

*     - Falg   (algebraic function containing the model)

*     PURPOSE

*     This subroutine calculates the states in the user given points (in x)

*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)
      real*8    states(nstatea,mostobs,nsaveset)
      real*8    xdata(nx,mostobs,nsets)
      integer*4 nobs(nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)

*     LOCAL VARIABLES

      integer*4 j,k,k1                            ! loop indices
      external Falg                               ! algebraic function
                                                  ! containing the model
      include 'common3.inc'

      if(debug .gt. 11) write(*,*) ' in salg'

      ncount=ncount+1
      k  = iset
      k1 = mod(k,nsaveset) + 1

      if (iobs .le. 0) then     ! all observations in the data set

         do j=1,nobs(k)

*           Variable iobs tells the current observation and may
*           be used in Falg.

            iobs = j

*           Calculate state values.

            call Falg(nstatea,states(1,j,k1),
     &                xdata(1,1,iset),nx,nobs(iset),
     &                gpar(1),ngpar,
     &                lpar(1,iset),nlpar,
     &                iobs,iset )

         enddo                ! j
         iobs = -1

      else       ! the  observation 'iobs' only

            call Falg(nstatea,states(1,iobs,k1),
     &                xdata(1,1,iset),nx,nobs(iset),
     &                gpar(1),ngpar,
     &                lpar(1,iset),nlpar,
     &                iobs,iset )


      endif
      
      return
      end

