************************************************************************
*                                                                      *
*     MODEST/Xparser                                                   *
*                                                                      *
*     Connections:                                                     *
*                                                                      *
*         Called by:   ...                                             *
*         Calls    :   parse, select                                   *
*                                                                      *
*     Description:                                                     *
*                                                                      *
*         This routine handles the selection and return of the         *
*         parameters between the optimizer and the tasks.              *
*                                                                      *
************************************************************************

      Subroutine Xparser(ar,nr,ai,ni,direction,xaux,naux,laux)


      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)


      character*3 direction          ! 'par' = parse, 'sel' = select
      Integer*4   naux
      Integer*4   laux   (4,naux)
      Real*8      xaux   (naux)           ! Parameters to be parsed

      integer i ,j

      Include    'common3.inc'     ! All MODEST common blocks

*     PERFORM THE OPERATION REQUESTED

      if (direction .eq. 'sel') then             ! select

              call select(ar,nr,ai,ni,xaux,naux,laux,
     &                  ar(pgpar),max(dngpar,1),
     &                  ar(plpar),max(dnlpar,1),
     &                  ar(pxdata),dnx,dmostobs,dnsets,
     &                  ar(pstates0),dnstatea)

      elseif (direction .eq. 'par') then         ! parse

              call parse(ar,nr,ai,ni,xaux,naux,laux,
     &                 ar(pgpar),max(dngpar,1),
     &                 ar(plpar),max(dnlpar,1),
     &                 ar(pxdata),dnx,dmostobs,dnsets,
     &                 ar(pstates0),dnstatea)

      else

              write(*,*) 'something wrong in direction'
              stop

      endif

*     ROUTINE COMPLETED

      return
      end


      Subroutine select(ar,nr,ai,ni,xaux,naux,laux,
     &                  gpar,ngpar,
     &                  lpar,nlpar,
     &                  xdata,nx,mostobs,nsets,
     &                  states0,nstatea)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)


      Integer*4  naux,ngpar,nlpar,nx,mostobs,nsets,nstatea
      Integer*4  laux(4,naux)          ! Pointers
      real*8     xaux(naux)            ! Parameters to be parsed
      real*8     gpar(ngpar)
      real*8     lpar(nlpar,nsets)
      real*8     xdata(nx,mostobs,nsets)
      real*8     states0(nstatea,nsets)

*     DEFINE INTERNAL PARAMETERS

      Integer*4  i                        ! Loop index

      Include    'common3.inc'       ! All MODEST common blocks

      if(debug .gt. 11) write(6,*) ' in select'

          do i=1,naux
              if(laux(1,i).eq.1) then
                  xaux(i) = gpar(laux(2,i))
              else if(laux(1,i).eq.2) then
                  xaux(i) = lpar(laux(2,i),laux(4,i))
              else if(laux(1,i).eq.3) then
                  xaux(i) = xdata(laux(2,i),laux(3,i),laux(4,i))
              else if(laux(1,i).eq.4) then
                  xaux(i) = states0(laux(2,i),laux(4,i))
              else if(laux(1,i).eq.5) then
                     xaux(i)=gpar(laux(2,i))
              else if(laux(1,i).ne.5) then
                  write(6,*) 'something wrong in lest or lopt'
                  stop
              endif
          enddo

      return
      end


      Subroutine parse(ar,nr,ai,ni,xaux,naux,laux,
     &                  gpar,ngpar,
     &                  lpar,nlpar,
     &                  xdata,nx,mostobs,nsets,
     &                  states0,nstatea)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      Integer*4  naux,ngpar,nlpar,nx,mostobs,nsets,nstatea
      Integer*4  laux(4,naux)          ! Pointers
      real*8     xaux(naux)            ! Parameters to be parsed
      real*8     gpar(ngpar)
      real*8     lpar(nlpar,nsets)
      real*8     xdata(nx,mostobs,nsets)
      real*8     states0(nstatea,nsets)

*     DEFINE INTERNAL PARAMETERS

      Integer*4  i                   ! Loop index

      Include    'common3.inc'       ! All MODEST common blocks

          do i=1,naux
              if(laux(1,i).eq.1) then
                  gpar(laux(2,i)) = xaux(i)
              else if(laux(1,i).eq.2) then
                  lpar(laux(2,i),laux(4,i)) = xaux(i)
              else if(laux(1,i).eq.3) then
                  xdata(laux(2,i),laux(3,i),laux(4,i)) = xaux(i)
              else if(laux(1,i).eq.4) then
                  states0(laux(2,i),laux(4,i)) = xaux(i)
              else if (laux(1,i).ne.5) then
                  write(6,*) 'something wrong in lest or lopt'
                  stop
              endif
          enddo

      return
      end



