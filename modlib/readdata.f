      subroutine Readdata (ar,nr,ai,ni,
     & xdata,nx,mostobs,nsets,
     & ydata,mostydat,
     & states0,nstatea,
     & t0,
     & tf,
     & nobs,
     & nydata,
     & weight)

*     This is the second input subroutine in the main program of MODEST.
*     The routine reads all the data not contained in the namelist file.

*     DATA FILES NEEDED

*     - datafile           data or weight or both
*     - weightfile              possible weight

*     SUBROUTINES NEEDED

*     - Readest            To read the data when the task is to estimate
*     - Readexp            To read the data when the task is to simulate or
*                          experimental design.
*     - Readsim            To choose between Readest and Readexp

      implicit none

*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4    nsets,nx,mostobs,mostydat,nstatea
      real*8       xdata(nx,mostobs,nsets)
      real*8       ydata(mostydat,mostobs,nsets)
      real*8       states0(nstatea,nsets)
      real*8       t0(nsets)
      real*8       tf(mostobs,nsets)
      integer*4    nobs(nsets)
      integer*4    nydata(mostobs,nsets)
      real*8       weight(mostydat,mostobs,nsets)



*     LOCAL PARAMETERS

      character*80  auxfile,auxfile1        ! aux name for datafile
      character*100 line                    ! a line in datafile
      character*100 blankol                 ! a blanko line

      integer*4     j,k,i1,j1,ntot          ! counters
      integer*2     ios                     ! for option error output
      logical       opened

      include 'common3.inc'

      blankol = '                    '
      blankol = blankol//blankol//blankol//blankol

*     Read the data from the datafile

      auxfile1 = datafile(1)
      opened = .false.
      do k=1,nsets

          if (combined.eq.1) then                           
             if (auxfile1(1:5).ne.'dummy' .and. k.eq.1) then
                open(unit=2,file = auxfile1, status = 'old')
                opened = .true.
             endif
          elseif (combined.eq.0) then
             opened = .false.
             if (auxfile1(1:5).ne.'dummy') then
                auxfile = datafile(k)
                if (auxfile(1:5).ne.'dummy') then
                   open(unit=2,file = auxfile, status = 'old')
                   opened = .true.
                endif
             endif
          endif

          iset = k

C         time points & possible data given in data sets
                                  
          if (opened) then
             if(debug .gt. 10) write(*,*) ' calling readfile'
             call readfile(ar,nr,ai,ni,
     &                    xdata,nx,mostobs,nsets,
     &                    ydata,mostydat,
     &                    weight,
     &                    t0,
     &                    tf,
     &                    nobs,
     &                    nydata,
     &                    ai(pncolxy),
     &                    ai(pindx),
     &                    ai(pindy),
     &                    ar(pxyaux),dmostxy,
     &                    ar(pstates),dnstatea,dnsavese,
     &                    states0,
     &                    ai(ps0file),
     &                    ar(pgpar),max(dngpar,1),
     &                    ar(plpar),max(dnlpar,1))

          else

*            all time points given in the namelist, no datafiles read

             xdata(1,1,k) = t0(k)
             do j = 2,nobs(k)
                if (tf(2,k).eq.0.0d0) then
*                  only t0 and tfinal given in namelist, set
*                  uniform time points between the first and last ones
                   xdata(1,j,k) = t0(k)+float(j-1)*tf(1,k)/
     &                                  float(nobs(k)-1)
                else
*                  all time points given in namelist
                   xdata(1,j,k) = tf(j-1,k)
                endif
             enddo

          endif

         if(combined .eq. 0) close(2)

      enddo


*     Read the weight for the components of y in each observation

      if(usew .eq.1) then

         do k = 1,nsets

         if(combined. eq. 0) open(unit=2,file = weightfile(k),
     &                            status = 'old')
            j1 = 1
            do while (j1 .le. nobs(k))

               read(2,'(a)') line

               if (line(1:1) .ne. '!' .and. line(1:1) .ne. '%' .and.
     &             line(1:1) .ne. '*' .and. line(1:1) .ne. 'c' .and.
     &             line(1:1) .ne. 'C' .and. line .ne. blankol) then

                   backspace(2)
                   read(2,*) (weight(i1,j1,k),
     &                          i1=1,nydata(j1,k))

                 j1 = j1 + 1
               endif
            enddo

          if(combined .eq. 0) close(2)

         enddo

      else                           ! usew = 0
         do k = 1,nsets
            do j = 1,nobs(k)
               call dset(nydata(j,k),1.0d0,weight(1,j,k),1)
            enddo
         enddo

      endif

      close(2)


      if(modelfile(1:5).ne.'dummy')
     &       call Fset(modelfile)

      return
      end


      subroutine Fset(txt1)

      character*80 txt1
      return
      end



