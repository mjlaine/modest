      subroutine cstats(ar,nr,ai,ni,
     &                  xaux,naux,laux,
     &                  nsets,nobs,nydata,
     &                  ydata,weight,mostydat,mostobs,
     &                  jtj,jtj1,u,s,v)

*     Calculates the approximate covariances and standard errors of the
*     estimated parameters

*     ARGUMENTS
      implicit none
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4 naux,nsets,nx,mostobs,mostydat
      real*8    xaux(naux)        ! xest or xopt
      Integer*4 laux(4,naux)      ! lest or lopt
      integer*4 nobs(nsets)
      integer*4 nydata(mostobs,nsets)
      real*8    ydata(mostydat,mostobs,nsets)
      real*8    weight(mostydat,mostobs,nsets)
      real*8    jtj(naux,naux),jtj1(naux,naux)
      real*8    u(naux,naux),s(naux),v(naux,naux)


*     LOCAL VARIABLES

      Integer*4    i,j,k,jj                       ! loop indices
      integer*4    ipath,info,irank
      Real*8       y                              ! ydata element
      Real*8       weigh                          ! corresponding weight
      Real*8       ss                             ! total corrected ss
      Real*8       ysum                           ! weighted sum of ydata
      Real*8       wtot                           ! sum of weight
      Real*8       sse                            ! residual ss
      Real*8       mse                            ! residual ms
      Real*8       stde                           ! std error of residuals
      real*8       tol
      real*8       sum
      Real*8       lsq                            ! see LSQ
      external     lsq

*     PARAMETERS AND COMMON BLOCKS

      include 'common3.inc'

*     calculate the total corrected sum of squares for the data
*     first take the weight into account

      if(debug .gt. 10) write(*,*) ' in cstats'

      ss   = 0.0d0
      ysum = 0.0d0
      wtot = 0.0d0

      if (model(1:3) .eq. 'ode' .and. obss0.eq.0) then !obss0(j,k)!!
         jj    = 2
      else
         jj    = 1
      endif

         do k=1,nsets
            do j=jj,nobs(k)
               do i = 1,nydata(j,k)

                  y      = ydata(i,j,k)
                  weigh  = weight(i,j,k)
                  wtot   = wtot + weigh
                  ysum   = ysum + weigh *y
                  ss     = ss + weigh *(y**2)

               enddo
            enddo
         enddo

*     correct ss by taking the mean into account

      ss = ss - (ysum ** 2) / wtot

*     calculate the mean sum of squares of the residuals

      sse  = lsq(ar,nr,ai,ni,xaux,naux)
      mse  = sse/(dntot-naux)
      stde = dsqrt(mse)

*     test for singularity & calculation of the covariances

      open(2,file=reportfile,status='unknown')
      write(2,formats) ' Total SS (corrected for means)',ss
      write(2,formats) ' Residual SS                   ',sse
      write(2,formats) ' Std. Error of estimate        ',stde

      write(2,'(//,1x,a,f7.2)') ' Explained (%):',(1.d0-sse/ss)*100.

      if (jacout.eq.0) then

          write(2,'(///)')

          write(2,*)
     &' Estimated  Parameters '

          write(2,*)

          do i=1,naux
              write(2,formatp) xaux(i)
          enddo


      elseif (jacout.eq.1) then

      write(2,*)
      write(2,*) ' The Hessian:'
      write(2,*)

      do i=1,naux
          write(2,formatm) (jtj(i,j),j=1,naux)
      enddo

*     compute svd of J'J = u*s*v'.

      ipath = 11
      tol   = 1d-16
      call dlsvrr(naux, naux, jtj, naux, ipath, info, tol, irank,
     &            s, u, naux, v, naux, jtj1,ar(pwrk2svd),ar(pwrksvd))

*    compute the appr Hessian mse * inv(J' * J) = mse * v*(1./s)*u'
*    for covariance interpretation

      do i = 1,naux
        call dscal(naux,1.0d0/s(i),v(1,i),1)
      enddo
      do i = 1,naux
      do j = 1,i
        sum = 0.0d0
        do jj = 1,naux
           sum = sum + v(i,jj)*u(j,jj)
        enddo
        jtj(i,j) = mse * sum
        jtj(j,i) = jtj(i,j)
      enddo
      enddo

          write(2,'(///)')

          write(2,*)
     &' Estimated   Estimated  Est. Relative Parameter/'

          write(2,*)
     &' Parameters  Std Error  Std Error (%) Std. Error'

          write(2,*)

          do i=1,naux
             if (jtj(i,i) .gt. 0.0d0) then
                write(2,formatp) xaux(i),dsqrt(jtj(i,i)),
     &               dabs(100.*dsqrt(jtj(i,i))/xaux(i)),
     &               dabs(xaux(i)/dsqrt(jtj(i,i)))
             else
                write(2,formatp) xaux(i),0.0,
     &               0.0,
     &               1.0d99
             endif
          enddo

          write(2,'(//)')

          write(2,*) ' The covariances of the parameters:'
          write(2,*)

          do i=1,naux
              write(2,formatm) (jtj(i,j),j=1,i)
          enddo

          write(2,*)
          write(2,*) ' The correlation matrix of the parameters:'
          write(2,*)

          do i=2,naux
              do j=1,i-1
                 if (jtj(i,i)*jtj(j,j) .gt. 0.0) then ! ML 2006-11-09
                    jtj(i,j) = jtj(i,j)/dsqrt(jtj(i,i)*jtj(j,j))
                 else
                    jtj(i,j) = 1.0d99  !
                 endif
                 jtj(j,i)=jtj(i,j)
              enddo
          enddo
          do i=1,naux
              jtj(i,i)=1.d0
          enddo

          do i=1,naux
              write(2,'(13f7.3)') (jtj(i,j),j=1,i)
          enddo

*     compute the eigenvalues by svd of correlation jtj

          ipath = 11
          tol   = 1d-16
          call dlsvrr(naux, naux, jtj, naux, ipath, info, tol, irank,
     &                s,u,naux,v,naux,jtj1,ar(pwrk2svd),ar(pwrksvd))

              write(2,'(//,1x,a,//)')
     &    ' The eigenvalues of the correlation matrix:'

          write(2,'(13f7.3)') (s(i),i=1,naux)

          write(2,'(//,1x,a,//)') ' The scaled principal components:'

          do i=1,naux
              write(2,'(13f7.3)') (u(i,j),j=1,naux)
          enddo

          close(2)

      endif ! jacout

      return

      end
