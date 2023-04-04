
      subroutine direct(ar,nr,ai,ni,
     &                  codir,ncodir,eigdir,neigdir,
     &                  dir,ndir,naux,nest,
     &                  jtj,u,s,v)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4 nest,ncodir,neigdir,ndir,naux
      integer*4 codir(ncodir)
      integer*4 eigdir(neigdir)
      real*8    dir(naux,ndir)
      real*8    jtj(nest,nest)
      real*8    u(nest,nest),s(nest),v(nest,nest)

*     local variables

      integer*4 i,idir,ipath,info,irank,j
      real*8    tol

      include 'common3.inc'

      idir = 0
      if (ncodir.gt.0) then
        do i = 1,ncodir
           dir(codir(i),i) = 1.0d0
        enddo
        idir = ncodir
      endif


      if (neigdir.gt.0) then

*        compute eigenvectors of J'J

         if (nest.eq.0) then
             write(*,*) 'nest is zero in DIRECT'
             stop
         endif
         if (debug .gt. 10 ) write(*,*) 'calling xparser in simulate'
         call xparser(ar,nr,ai,ni,'sel',ar(pxest),nest,ai(plest))
         call hess(ar,nr,ai,ni,
     &             jtj,u,ar(pxest),nest,ai(plest))

         if (debug .gt. 10 ) write(*,*) 'calling dlsvrr in simulate'
         ipath = 11
         tol   = 1d-16
         call dlsvrr(naux, naux, jtj, naux, ipath, info, tol, irank,
     &               s,u,naux,v,naux,
     &               u,ar(pwrk2svd),ar(pwrksvd))

*        add the eigendirections to 'dir'

         do j = 1,neigdir
         do i = 1,naux
            dir(i,idir+j) = v(i,eigdir(j))
         enddo
         enddo
         idir = idir + neigdir

      endif

      if (ndir.ne.idir) then
         write(*,*) 'ndir ne sumdir'
         stop
      endif


      return
      end

