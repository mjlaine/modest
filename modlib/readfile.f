
      subroutine Readfile(ar,nr,ai,ni,
     & xdata,nx,mostobs,nsets,
     & ydata,mostydat,
     & weight,
     & t0,
     & tf,
     & nobs,
     & nydata,
     & ncolxy,
     & indx,
     & indy,
     & xyaux,mostxy,
     & states,nstatea,nsaveset,
     & states0,s0file,
     & gpar,ngpar,
     & lpar,nlpar)

*     This subroutine reads the input data from datafile.
*     The routine is CALLED BY Readata when the task is 'estimate'

*     DATA FILES NEEDED

*     - datafile           data or weight or both

*     ARGUMENTS
      implicit none
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4    nx,mostobs,nsets,mostydat,mostxy,nstatea,nsaveset
      integer*4    ngpar,nlpar

      real*8       xdata(nx,mostobs,nsets)
      real*8       ydata(mostydat,mostobs,nsets)
      real*8       weight(mostydat,mostobs,nsets)
      real*8       t0(nsets)
      real*8       tf(mostobs,nsets)
      integer*4    nobs(nsets)
      integer*4    nydata(mostobs,nsets)
      integer*4    ncolxy(nsets)
      integer*4    indx(nx,nsets)
      integer*4    indy(mostydat,nsets)
      real*8       xyaux(mostxy)
      real*8       states(nstatea,mostobs,nsaveset)
      real*8       states0(nstatea,nsets)
      integer*4    s0file(nstatea,nsets)
c      real*8       yest(mostydat,mostobs,nsaveset)
      real*8       gpar(ngpar)
      real*8       lpar(nlpar,nsets)

*     LOCAL VARIABLES

      integer*4     i,j,k,k1,i1,j1,i2,ifile    ! counters
      integer*4     ncolaux                 ! n of col in data files
      integer*4     nstates
      real*8        t0aux
      character*100 line                    ! a line in datafile
      character*100 blankol                 ! a blanco line
      logical       casey                   ! false: sim/exp/opt
                                            ! true : est/sen

*     COMMON BLOCS

      include 'common3.inc'

      if(debug .gt. 11) write(*, *) ' in readfile'

      k = iset
      if (task(1:3).eq.'est' .or.
     &  (task(1:3).eq.'sen' .and.
     &  (objfun(1:3).eq.'lsq' .or. objfun(1:2).eq.'r2'))) then
         casey = .true.
      else
         casey = .false.
      endif
      blankol = '                                                      '
      nstates = nstatea-dnsaux

      if (casey) then
          ncolaux = ncolxy(k)
      else
          ncolaux = 0
          do i = 1,nx
             ncolaux = max(ncolaux,indx(i,k))
          enddo
      endif

      if(echo.eq.1) then
         write(6,*) 'dataset, nobs, ncolxy'
         write(6,*)  k, nobs(k), ncolxy(k)
         write(6,*) '(indx)'
         write(6,*)  (indx(i1,k), i1=1,nx)
         if (casey) then
           write(6,*) 'obs, nydata, (indy)'
           do j = 1,nobs(k)
              write(6,*)j,nydata(j,k),(indy(i1,k),i1=1,nydata(j,k))
           enddo
         endif
      endif

*     Read the data from the datafile and place it in the
*     x- and y-matrixes.

      if(model(1:3) .eq. 'alg'.or. model(1:3) .eq. 'imp') then

         j1 = 1

         do while (j1 .le. nobs(k))

            read(2,'(a)') line
            if (line(1:1) .ne. '!' .and. line(1:1) .ne. '%' .and.
     &          line(1:1) .ne. '*' .and. line(1:1) .ne. 'c' .and.
     &          line(1:1) .ne. 'C' .and. line .ne. blankol) then

                backspace(2)
                read(2,*) (xyaux(i1), i1=1,ncolaux)

                do i2 = 1,nx
                   xdata(i2,j1,k) = xyaux(indx(i2,k))
                enddo

                if (casey) then
                  do i2 = 1,nydata(j1,k)
                     ydata(i2,j1,k) = xyaux(indy(i2,k))
                  enddo
                endif

                j1 = j1 + 1
            endif
         enddo

       else if (model(1:3) .eq. 'ode') then

*      see whether initial values in file

         ifile = 0                 !n of s0 to be found in file k
         do i =1,nstatea
            if (s0file(i,k).eq.1) then
               ifile = ifile + s0file(i,k)
               s0file(ifile,k) = i !use the used s0file places as index
            endif
         enddo
         if (ifile.ne.0 .and. ifile .lt. nydata(1,k)) then
          write(*,*) 'the number of y columns is greater than '
          write(*,*) 'the number of initial values to be found in file'
          write(*,*) 'in data set  ',k
          stop
c          do i = 1,nydata(1,k)
c             s0file(i,k) = i
c          enddo
         elseif (ifile.ne.0 .and. ifile .gt. nydata(1,k)) then
          write(*,*) 'the number of y columns is less than '
          write(*,*) 'the number of initial values to be found in file'
          write(*,*) 'in data set  ',k
          stop
         elseif (ifile.eq.0 .and. nx.gt.1) then
          write(*,*) 'All values for control variables must be given'
          write(*,*) 'in data file. So give the starting point of '
          write(*,*) 'integration and (some) initial values in file'
          stop
         endif

         if(ifile.gt.0) then

*        Initial values first (j1 = 1)

            j1 = 1

            do while (j1 .eq. 1)

                read(2,'(a)') line

                if (line(1:1) .ne. '!' .and. line(1:1) .ne. '%' .and.
     &              line(1:1) .ne. '*' .and. line(1:1) .ne. 'c' .and.
     &              line(1:1) .ne. 'C' .and. line .ne. blankol) then

                   backspace(2)
                   read(2,*) (xyaux(i1), i1=1,ncolxy(k))

                    do i2 = 1,nx
                       xdata(i2,1,k) = xyaux(indx(i2,k))
                    enddo

c                    write(*,*) 'k',k,nydata(1,k),mostydat
c                    write(*,*) indy(100,1),indy(100,2),indy(100,3)

                    do i2 = 1,nydata(1,k)
                       states0(s0file(i2,k),k) = xyaux(indy(i2,k))
                    enddo

                    j1 = j1 + 1

               endif
            enddo

            if (casey) then
              do i2 = 1,nydata(1,k)
                 ydata(i2,1,k) = xyaux(indy(i2,k))
              enddo
            endif

         elseif(ifile.eq.0) then

            xdata(1,1,k) = t0(k)

            if(debug .gt. 10) write(*, *) ' calling observations'

*           Just to substitute 'ydata=yest' for plotting

*****************************************
      if(debug .gt. 10) write(*,*) ' calling inits0'
            k1   = mod(k,nsaveset)+1
            call inits0(nstates,t0aux,states(1,1,k1),
     &       xdata(1,1,k),nx,nobs(k),dnsaux,nstates+dnsaux,
     &       states0(1,k),
     &       gpar,max(dngpar,1),
     &       lpar(1,k),max(dnlpar,1),
     &       1,k)
*****************************************
            call observations(states(1,1,k1),nstates,
     &                        ydata(1,1,k),nydata(1,k),
     &                        xdata(1,1,k),nx,nobs(k),
     &                        dnsaux,nstates+dnsaux,
     &                        states0(1,k),
     &                        gpar(1),max(dngpar,1),
     &                        lpar(1,k),max(dnlpar,1),
     &                        1,k)

c            if (casey) then
               do j = nobs(k),1,-1
                 nydata(j+1,k) = nydata(j,k)
               enddo
                 nydata(1,k) = nydata(2,k)
c            endif

            nobs(k) = nobs(k) + 1    ! n of obs. altogether

         endif

*        the rest then (j1 = 2,...,nobs(k))

         j1 = 2

         do while (j1 .le. nobs(k))
            read(2,'(a)') line
            if (line(1:1) .ne. '!' .and. line(1:1) .ne. '%' .and.
     &          line(1:1) .ne. '*' .and. line(1:1) .ne. 'c' .and.
     &          line(1:1) .ne. 'C' .and. line .ne. blankol) then

                   backspace(2)
                   read(2,*) (xyaux(i1), i1=1,ncolaux)

                do i2 = 1,nx
                   xdata(i2,j1,k) = xyaux(indx(i2,k))
                enddo

                if (casey) then
                   do i2 = 1,nydata(j1,k)
                     ydata(i2,j1,k) = xyaux(indy(i2,k))
                   enddo
                endif

                j1 = j1 + 1

            endif
         enddo

      endif

      if(echodata .eq. 1) then

         if(model(1:3) .eq. 'ode') then

            write(*,*) 'initial values'
            write(*,*) (states0(i1,k),i1=1,nstates)

         endif

            write(*,*) 'x variables'
            do j=1,nobs(k)
               write(*,*) (xdata(i1,j,k),i1=1,nx)
            enddo

         if (casey) then
            write(*,*) 'data'
            do j=1,nobs(k)
               write(*,*) (ydata(i1,j,k),i1=1,nydata(j,k))
            enddo
         endif

      endif

      return
      end

