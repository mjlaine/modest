
      subroutine Readbini(ar,nr,ai,ni)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

*     This is the input subroutine for the main program of MODEST.
*     The routine reads the data from the temporary unformatted file modtemp.in.

*     DATA FILES NEEDED

*     - modtemp.in          The data in unformatted form

*     LOCAL PARAMETERS

      integer*4  i
      Integer*4  nstates
      Integer*4  nx

      Integer*4  mostobs
      Integer*4  mostydat
      Integer*4  mostxy
      Integer*4  nsaveset
      Integer*4  nsaux
      Integer*4  nstatea
      Integer*4  nsets
c     Integer*4  nfiles
      Integer*4  ntot
      Integer*4  ngpar
      Integer*4  nlpar
      Integer*4  nsim
      Integer*4  nest
      Integer*4  nexp
      Integer*4  nsen
      Integer*4  nopt
      Integer*4  nnew
      Integer*4  nexpp
      Integer*4  nesopsi
      Integer*4  moststep
      Integer*4  ndir
      Integer*4  ncodir
      Integer*4  neigdir
      Integer*4  ngrid
      Integer*4  nfilepoint
      Integer*4  ndiscpoint
      Integer*4  ncontour
      Integer*4  liwddr
      Integer*4  lrwddr
      Integer*4  liworkod
      Integer*4  lrworkod
      Integer*4  latolod
      Integer*4  lrtolod
      Integer*4  lwanr,liwanr,ljacnr
c      Integer*4  lwadne
      integer*4  NCOAUX
      integer*4  NESOPSI1
      integer*4  NYH
      integer*4  LWM
      integer*4  MAXORD
      character*10 line1


*     PARAMETERS AND COMMON BLOCKS

      include 'common3.inc'

      open(unit=1,file = tempin,status = 'old')

      read(1,'(a)') line1
      if (line1(1:4) .eq. 'stop') then
       write(*,*) 'No proper input file ',tempin
       write(*,*) ' '
       stop
      elseif (line1(1:7) .eq. 'changed') then
       write(*,*) 'Essential changes in modelvar since the last run.'
       write(*,*) 'You probably should recompile the code'
       write(*,*) ' - but I only warn you once! '
       write(*,*) ' '
       stop
      elseif (line1(1:7) .eq. 'missing') then
       write(*,*) 'No previous declarations for modelvar.'
       write(*,*) 'Check if code should be compiled'
       write(*,*) ' - otherwise just run again. '
       write(*,*) ' '
       stop

      elseif (line1(1:5) .eq. 'ready') then
c      nothing to complain ...
      endif


      do i = 1,ndim
         read(1,*) dim(i)
      enddo
      do i = 1,nscali
         read(1,*) scali(i)
      enddo
      do i = 1,nscalr
         read(1,*) scalr(i)
      enddo
      do i = 1,nchar
         read(1,'(a)') char(i)
      enddo
      if (combined.eq.1) then
        read(1,'(a)') datafile(1)
      elseif (combined.eq.0) then
        read(1,'(a)') datafile(1)
c ML 30.4.2002 nsets -> dnsets
        if (dnsets.gt.maxnfile .and. datafile(1).ne.'dummy.dat') then
           write(*,*) 'max', maxnfile, '  separate datafiles allowed'
           write(*,*) 'for more data sets:  combine in one file'
           write(*,*)
           stop
        endif

        if (datafile(1).ne.'dummy.dat') then
          do i = 2,dnsets
            read(1,'(a)') datafile(i)
          enddo
          do i = 1,dnsets
            read(1,'(a)') weightfile(i)
          enddo
        endif
      endif

c      write(*,*) 'debug',debug
c      write(*,*) 'echo',echo
c      write(*,*) 'echodata',echodata
c      read(*,*) echodata

        mostobs     =  dmostobs
        mostydat    =  dmostyda
        mostxy      =  dmostxy
        nsaveset    =  dnsavese
        nsaux       =  dnsaux
        nstatea     =  dnstatea
        nsets       =  dnsets
c       nfiles      =  dnfile
        mostydat    =  dmostyda
        mostobs     =  dmostobs
        ntot        =  dntot
        ngpar       =  dngpar
        nlpar       =  dnlpar
        nsim        =  dnsim
        nest        =  dnest
        nexp        =  dnexp
        nsen        =  dnsen
        nopt        =  dnopt
        nnew        =  dnnew
        nexpp       =  dnexpp
        nesopsi     =  dnesopsi
        moststep    =  dmostste
        ndir        =  dndir
        ncodir      =  dncodir
        neigdir     =  dneigdir
        ngrid       =  dngrid
        nfilepoint  =  dnevapoi
        ndiscpoint  =  ddiscpoi
        ncontour    =  dncontou
        liwddr      =  dliwddr
        lrwddr      =  dlrwddr
        lrworkod    =  dlrwrkod
        liworkod    =  dliwrkod
        lrtolod     =  dlrtolod
        latolod     =  dlatolod
        lwanr       =  dwanr
        liwanr      =  diwanr
        ljacnr      =  djacnr


        nx          = dnx
        nstates     = nstatea - nsaux
        nesopsi1    = nesopsi + 1
        ncoaux      = max(nsim,nsen)


*******************Start allocations ************************

         rdataptr = 1
         idataptr = 1

c     states in ode/imp solvers:

      if (model(1:3) .eq.'ode') then
c 24.10.1999 ML
         if (odesolver(1:4).eq.'odes' 
     &       .or. odesolver(1:4).eq.'lsod') then
            call ralloc(ar,nr,0,pyod,nstatea,nest+1,1,1)
         elseif (odesolver(1:4) .eq. 'ddri' .or.
     &           odesolver(1:4) .eq. 'eule') then
            call ralloc(ar,nr,0,psdr2,nstatea,1,1,1)
         endif
      elseif (model(1:3) .eq. 'imp') then
            call ralloc(ar,nr,0,psimp,nstatea,1,1,1)
      endif

c    states etc:

      call ralloc(ar,nr,0,pstates,nstatea,mostobs,nsaveset,1)
      call ralloc(ar,nr,1,pstates0,nstatea,nsets,1,1)
      call ialloc(ai,ni,1,ps0file, nstatea,nsets,1,1)
      call ralloc(ar,nr,1,pt0,nsets,1,1,1)
      call ralloc(ar,nr,1,ptf,mostobs,nsets,1,1)
      call ralloc(ar,nr,0,pyest,mostydat,mostobs,nsaveset,1)
      call ralloc(ar,nr,0,pxjac,mostydat,nest,mostobs,nsaveset)

c    pointed variables and pointers:

      call ralloc(ar,nr,1,pgpar,ngpar,1,1,1)
      call ralloc(ar,nr,1,plpar,nlpar,nsets,1,1)
      call ialloc(ai,ni,1,plsim,4,nsim,1,1)
      call ialloc(ai,ni,1,plest,4,nest,1,1)
      call ialloc(ai,ni,1,plexp,4,nexp,1,1)
      call ialloc(ai,ni,1,plsen,4,nsen,1,1)
      call ialloc(ai,ni,1,plopt,4,nopt,1,1)
      call ialloc(ai,ni,0,plaux,4,nesopsi,1,1)

c    generally used vectors:

      call ralloc(ar,nr,0,pxsim,nsim,1,1,1)
      call ralloc(ar,nr,0,pxest,nest,1,1,1)
      call ralloc(ar,nr,0,pxexp,nexp,1,1,1)
      call ralloc(ar,nr,0,pxsen,nsen,1,1,1)
      call ralloc(ar,nr,0,pxopt,nopt,1,1,1)
      call ralloc(ar,nr,0,pxaux,nesopsi,1,1,1)
      call ralloc(ar,nr,1,pbounds,2,nesopsi,1,1)
      call ralloc(ar,nr,1,psbounds,2,nesopsi,1,1)
      call ialloc(ai,ni,1,pbtype,nesopsi,1,1,1)
      call ralloc(ar,nr,0,pjtj,nest,nest,1,1)
      call ralloc(ar,nr,0,pjtj1,nest,nest,1,1)
      call ralloc(ar,nr,0,peigval,nest,1,1,1)
      call ralloc(ar,nr,0,peigvec,nest,nest,1,1)

c    data indexes:
      call ralloc(ar,nr,1,prnobs,nsets,1,1,1)
***********
      call ialloc(ai,ni,1,pnobs,nsets,1,1,1)
      call ialloc(ai,ni,1,pnydata,mostobs,nsets,1,1)
      call ialloc(ai,ni,1,pncolxy,nsets,1,1,1)
      call ialloc(ai,ni,1,pindx,nx,nsets,1,1)
      call ialloc(ai,ni,1,pindy,mostydat,nsets,1,1)


      if (task(1:3) .eq. 'exp' .or.
     &   (task(1:3) .eq. 'sen' .and. objfun(1:3).eq.'exp')) then

        call ralloc(ar,nr,0,pexpdu,nest,nest,1,1)
        call ralloc(ar,nr,0,pexpdv,nest,nest,1,1)
        call ialloc(ai,ni,1,pexppar,nexpp,1,1,1)
        call ralloc(ar,nr,0,pjtj0,nest,nest,1,1)
        call ralloc(ar,nr,0,pjtj2,nest,nest,1,1)
        call ralloc(ar,nr,0,psigma,nest,1,1,1)
        call ialloc(ai,ni,1,pnewobs,nnew,2,1,1)
        call ialloc(ai,ni,0,piauxsub,nesopsi,1,1,1)

      endif

        call ialloc(ai,ni,1,pnsimste,ndir,1,1,1)
        call ialloc(ai,ni,1,pcodir,ncodir,1,1,1)
        call ialloc(ai,ni,1,peigdir,neigdir,1,1,1)
        call ralloc(ar,nr,0,pcoline,moststep,ncodir,1,1)
        call ralloc(ar,nr,0,pdir,nesopsi,ndir,1,1)
        call ralloc(ar,nr,0,pxdir,nesopsi,moststep,1,1)
        call ralloc(ar,nr,0,pxdisc,nesopsi,ndiscpoint,1,1)
        call ralloc(ar,nr,0,pvaldisc,ndiscpoint,1,1,1)
        call ralloc(ar,nr,0,poptline,moststep,ndir,1,1)
        call ralloc(ar,nr,0,poptgrid,ncontour,moststep,moststep,1)
        call ialloc(ai,ni,1,pigrid,2,ngrid,1,1)
        call ralloc(ar,nr,0,ptlaux,nesopsi,1,1,1)
        call ralloc(ar,nr,0,ptuaux,nesopsi,1,1,1)


c    arrays used in optimizer subroutines:

      if (task(1:3).eq.'est' .or. task(1:3).eq.'exp' .or.
     &    task(1:3).eq.'opt' .and. optimizer(1:5).eq.'simfl') then


        call ralloc(ar,nr,0,psimplex,nesopsi1,nesopsi1,1,1)
        call ralloc(ar,nr,0,pwrksimf,5*nesopsi+2,1,1,1)
        call ialloc(ai,ni,0,piwrsimf,2*nesopsi+2,1,1,1)

      endif

      if (task(1:3).eq.'est' .and. (optimizer(1:3).eq.'lev' .or.
     &                              optimizer(1:5).eq.'simle')) then
        call ralloc(ar,nr,0,pzlm,nest,3,1,1)
        call ialloc(ai,ni,0,pipivlm,nest,1,1,1)
      endif

c     arrays used in the solver subroutines:

c     Ode solvers

      if (model(1:3) .eq. 'ode' .and.
     &    odesolver(1:4) .eq. 'ddri') then

c     ddriv2
         if (mintdr2 .eq. 1) then
             liwddr = 21
             lrwddr = 16*nstates + 2*nrootdr2 + 204
         elseif (mintdr2 .eq. 2) then
             liwddr = nstates+21
             lrwddr = nstates*nstates+10*nstates+2*nrootdr2+204
         elseif (mintdr2 .eq.3) then
             liwddr = nstates+21
             lrwddr = nstates*nstates+17*nstates+2*nrootdr2+204
         endif

        call ialloc(ai,ni,0,piwrkddr,liwddr,1,1,1)
        call ralloc(ar,nr,0,pworkddr,lrwddr,1,1,1)

      elseif (model(1:3) .eq. 'ode' .and.
     &       ( odesolver(1:4) .eq. 'odes' 
     &     .or. odesolver(1:4).eq.'lsod')) then

c     Odessa

         if (isoptod.eq.0) then
            nyh = nstatea
         elseif (isoptod.eq.1 .and. nest.gt.0) then
            nyh     = nstatea * (nest + 1)
         else
            write(*,*) 'check ISOPT and ESTVAR for ODESSA'
            stop
         endif

         if (method.eq.1) then
            maxord = 12
         elseif (method.eq.2) then
            maxord = 5
         endif

         if (miterod.eq.0) then
            lwm = 0
         elseif (miterod.eq.1 .or. miterod.eq.2) then
            lwm = nstatea*nstatea+2
         elseif (miterod.eq.3) then
            lwm = nstatea+2
         elseif (miterod.eq.4 .or. miterod.eq.5) then
            lwm = (2*mlod+muod+1)*nstatea+2
         else
            write(*,*) 'mjacob incorrect for ODESSA'
         endif

#ifdef USE_LSODE
         if ( odesolver(1:4).eq.'lsod' ) then
ccc for mf = 222 (miten lasketaan lwm???) ML 29.3.2001
c            lwm = 2*nstatea*nstatea + 2*nstatea 
c            lwm = lwm + (nstatea*nstatea+10*nstatea)/2
		lwm = 20 + 2*nstatea*nstatea + 16*nstatea
            lrworkod = 20 + 9*nstatea + lwm 
            liworkod = 30
c
         else                   ! odessa
#endif
            if (isoptod.eq.0) then
               lrworkod = 20 + nyh*(maxord+1) + 3*nstatea + lwm
               if (miterod.eq.0 .or. miterod.eq.3) then
                  liworkod = 20
               else
                  liworkod = 20 + nstatea
               endif
            elseif (isoptod.eq.1) then
               lrworkod = 20 + nyh*(maxord+1) + 2*nyh + lwm + nstatea
               liworkod = 21 + nstatea + nest
            endif

#ifdef USE_LSODE
         endif
#endif

         if (itolod.eq.1 .or. itolod.eq.2) then
            lrtolod = 1
         else
            lrtolod = nstatea*(nest+1)
         endif

         if (itolod.eq.1 .or. itolod.eq.3) then
            latolod = 1
         else
            latolod = nstatea*(nest+1)
         endif


        call ralloc(ar,nr,0,prworkod,lrworkod,1,1,1)
        call ialloc(ai,ni,0,piworkod,liworkod,1,1,1)
        call ralloc(ar,nr,0,prtolod,lrtolod,1,1,1)
        call ralloc(ar,nr,0,patolod,latolod,1,1,1)
        call ralloc(ar,nr,0,pteta0od,nest,1,1,1)

      elseif (model(1:3) .eq. 'ode' .and.
     &        odesolver(1:4) .eq. 'eule') then

        call ralloc(ar,nr,0,pworkeul,2*nstates,1,1,1)

      endif

c    implicit alg solvers

      if (model(1:3) .eq. 'imp') then
       call ralloc(ar,nr,1,psmin,nstates,1,1,1)
       call ralloc(ar,nr,1,psmax,nstates,1,1,1)
       if (equsolver(1:3).eq.'new') then
          lwanr  = 2*nstates*nstates + 10*nstates
          ljacnr = nstates*nstates         !no sparsity here
          liwanr = 2*ljacnr + 9*nstates+2
          call ralloc(ar,nr,0,pwanr,lwanr,1,1,1)
          call ialloc(ai,ni,0,piwanr,liwanr,1,1,1)
          call ialloc(ai,ni,0,piparnr,10,1,1,1)
       else
c       call ralloc(ar,nr,0,pfvecdne,nstates,1,1,1)
c       lwadne  = 3*nstates*(nstates+5)
c       call ralloc(ar,nr,0,pwadne,lwadne,1,1,1)
        write(*,*) 'No such nonlinear equation solver implemented'
        stop
       endif
      endif

      call ralloc(ar,nr,0,pesvd,nest,1,1,1)
      call ralloc(ar,nr,0,pwrksvd,nest,1,1,1)
      call ralloc(ar,nr,0,pwrk2svd,nest,nest,1,1)
******endif

******if ...
      call ralloc(ar,nr,0,pwrkdsco,nesopsi1,nesopsi1,1,1)
******endif

******if ...   !lsqlevmar,  lsqxaux(nest)
c     call ralloc(ar,nr,0,plsqxaux,nest,1,1,1)

c    arrays for data:

      call ralloc(ar,nr,0,pxdata,nx,mostobs,nsets,1)
      call ralloc(ar,nr,0,pydata,mostydat,mostobs,nsets,1)
      call ralloc(ar,nr,0,pxyaux,nx+mostydat,1,1,1)
      call ralloc(ar,nr,0,pweight,mostydat,mostobs,nsets,1)

      if (debug.eq.21) then
        write(*,*) 'real, integer'
        write(*,*)  rdataptr,idataptr
        read(*,*)  i
      endif
************** ALLOCATIONS DONE, update dim  ************************

      dnstates =  nstates
      dnx      =  nx

      dmostobs =  mostobs
      dmostyda =  mostydat
      dmostxy  =  mostxy
      dnsavese =  nsaveset
      dnsaux   =  nsaux
      dnstatea =  nstatea
      dnsets   =  nsets
c     dnfile   =  nfiles
      dmostyda =  mostydat
      dmostobs =  mostobs
      dntot    =  ntot
      dngpar   =  ngpar
      dnlpar   =  nlpar
      dnsim    =  nsim
      dnest    =  nest
      dnexp    =  nexp
      dnsen    =  nsen
      dnopt    =  nopt
      dnnew    =  nnew
      dnexpp   =  nexpp
      dnesopsi =  nesopsi
      dmostste =  moststep
      dndir    =  ndir
      dncodir  =  ncodir
      dneigdir =  neigdir
      dnevapoi =  nfilepoint
      ddiscpoi =  ndiscpoint
      dncontou =  ncontour
      dliwddr  =  liwddr
      dlrwddr  =  lrwddr
      dlrwrkod =  lrworkod
      dliwrkod =  liworkod
      dlrtolod =  lrtolod
      dlatolod =  latolod
      dwanr    =  lwanr
      diwanr   =  liwanr
      djacnr   =  ljacnr

c      dlwadne  =  lwadne

      close(1)
*********************************************************************
      return
      end

      subroutine ralloc(ar,nr,type,ptr,dim1,dim2,dim3,dim4)
c
c     This subroutine allocates space for the integer and real
c     arrays used in modest program and/or puts the available
c     values in their place in the real vector 'ar'.
c           type   0     only compute 'ptr' and update 'rdataptr'.
c                  1     compute 'ptr', update 'rdataptr' and
c                        substitute values
c                  2     only substitute values, no updating
c           ptr          pointer to start of vector in 'ar'
c                        output if type = 0,1; input if type = 2
c           dim's        the actual dimensions (default=1)

      real*8    ar(nr)
      integer*4 type,ptr,dim1,dim2,dim3,dim4,dd

      integer*4 i,j,k,l
      include 'common3.inc'

      dd = dim1*dim2*dim3*dim4
      if (dd.eq.0) then
         ptr = rdataptr
          return
      endif

      if (type .lt. 2) ptr = rdataptr   !get ptr

      if (type.eq.1) then
         do l=0,dim4-1
         do k=0,dim3-1
         do j=0,dim2-1
         do i=0,dim1-1
             read(1,*) ar(ptr+i+j*dim1+k*dim1*dim2+l*dim1*dim2*dim3)
         enddo
         enddo
         enddo
         enddo
      endif
      if (type.lt.2) rdataptr=rdataptr + dd  !update rdataptr

      return
      end

      subroutine ialloc(ai,in,type,ptr,dim1,dim2,dim3,dim4)
c
c     This subroutine allocates space for the integer and real
c     arrays used in modest program and/or puts the available
c     values in their place in the integer vector 'ai',
c           type   0     only compute 'ptr' and update 'rdataptr'.
c                  1     compute 'ptr', update 'rdataptr' and
c                        substitute values
c                  2     only substitute values, no updating
c           ptr          pointer to start of vector in 'ar'
c                        output if type = 0,1; input if type = 2
c           dim's        the actual dimensions (default=1)

      integer*4 ai(in)
      integer*4 type,ptr,dim1,dim2,dim3,dim4,dd

      integer*4 i,j,k,l
      include 'common3.inc'

      dd = dim1*dim2*dim3*dim4
      if (dd.eq.0) then
         ptr = idataptr
         return
      endif

      if (type .lt. 2) ptr = idataptr   !get ptr

      if (type.eq.1) then
         do l=0,dim4-1
         do k=0,dim3-1
         do j=0,dim2-1
         do i=0,dim1-1
             read(1,*) ai(ptr+i+j*dim1+k*dim1*dim2+l*dim1*dim2*dim3)
         enddo
         enddo
         enddo
         enddo
      endif
      if (type.lt.2) idataptr=idataptr+dd  !update rdataptr

      return
      end

