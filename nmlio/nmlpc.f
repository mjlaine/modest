      program nmlio

      implicit none

********Leading dimension declarations for namelist input (only!) *********

      integer*4 maxnstat,maxnobs,maxnsets,maxnstep,maxord,maxndir,
     &          maxny,maxs0
calloc      parameter (maxnstat =  300)           ! max n of equations (& responses)
calloc      parameter (maxnobs  =  2000)          ! max n of observations/data set
calloc      parameter (maxnsets =  100)           ! max n of datasets

      integer*4 maxngpar,maxnlpar
      parameter (maxngpar =  1000)           ! max n of global parameters
      parameter (maxnlpar =  500)           ! max n of local parameters

      integer*4 maxnx,maxnopt,maxnopt2,maxnopt4,maxnnew
      parameter (maxnx     = 500)          ! max n of independent variables
      parameter (maxnopt   = 100)           ! max n of parameters in
                                            ! optimization, simulation, sensitivity
      parameter (maxnopt2  = 2*maxnopt,maxnopt4 = 4*maxnopt)
      parameter (maxndir   = maxnopt2)
calloc      parameter (maxnnew   = maxnopt*maxnsets)
calloc      parameter (maxs0     = maxnsets*maxnstat)
calloc      parameter (maxny     = maxnsets*maxnobs)

***************SMBL variables **********
      integer*4  mxmsent,mxwords ,mxchars
      parameter(mxmsent = 2000,    !max numb. of modelvar-sent.
     *          mxwords = 6000,    !max numb. of words in input-strs
     *          mxchars = 400000)  !max length of strs modelvar,...
c       mxchars = max(len(modelvar,len(target),len(estvar))
       integer*4 error_number,ngrid

calloc      integer*4 ntf(maxnsets)
      integer*4 ntf[ALLOCATABLE](:)

       character*7 type(mxmsent)
       integer charcl(0:mxchars+1),nglobal(mxmsent),
     *          nlocal(mxmsent),ncontrol(mxmsent),
     *          ninitval(mxmsent)
       character*80 word(mxwords),sentence(mxwords)
       character*80 name(mxmsent)

cmatrix-input:
       integer vardim(mxmsent,2)

*******************************
       integer mustcmp
       character*3 ftype
       integer inttype(mxmsent),iperiod
*******************************
       
      real*8  uninit
      integer*4 unini
      parameter (uninit = -9999.0d0)
      parameter (unini  = -9999)

************** LOCAL PARAMETERS ***********************************

      character*80  nmlfile                 ! name of the *.nml file
      character*80  auxfile                 ! aux name for datafile
      character*100 line                    ! a line in datafile
      character*7000 modelvar
      character*1990 target,estvar          ! symbolic names in nml
      character*500 codirvar,gridvar        !   "

      integer*4  error
      real*8     aux
      real*8     xaux,xest,xopt,optgrid
c***      real*8     time[ALLOCATABLE](:,:)
      integer*4  iaux,ntot
      integer*4  i,j,k                      ! counters
      character*80 incfile                  ! project.inc
      character*100 blankol                 ! a blank line
      integer*4  nnew,nexpp,ncoaux
      logical    objunini

********** Namelist variables, scalars in 'common3.inc' ***************
********** tables to be saved defined here  ***********************************

       integer*4    lsim(4,maxnopt),lest(4,maxnopt),lexp(4,maxnopt),
     &              lsen(4,maxnopt),lopt(4,maxnopt),laux(4,maxnopt)
       integer*4    nsim,nest,nexp,nsen,nopt
       integer*4    ndir,ncodir,neigdir,ncontour,nbounds
       integer*4    nfilepoint,ndiscpoint
       integer*4    igrid(2,maxnopt)
       integer*4    codir(maxnopt)
       integer*4    eigdir(maxnopt)

       real*8       bounds(2,maxnopt), sbounds(2,maxnopt)
       integer*4    nsbound

       integer*4    btype(maxnopt)
       integer*4    nevastep(maxndir)
       integer*4    newobs[ALLOCATABLE](:,:)
calloc       integer*4    newobs(maxnnew,2)
*     Model parameters
       integer*4    ngpar,nlpar,nstates,nx,nsaux
       real*8       gpar(maxngpar)
       real*8       lpar[ALLOCATABLE](:,:)
calloc       real*8       lpar(maxnlpar,maxnsets)
       integer*4    ifile,nstates0
       integer*4    s0file[ALLOCATABLE](:,:)
calloc       integer*4    s0file(maxnstat,maxnsets)

*     Parameters of data
       integer*4    nsets,nsaveset

*     Data set parameters
       real*8       rnobs[ALLOCATABLE](:)  !******
       integer*4    nobs[ALLOCATABLE](:),ncolxy[ALLOCATABLE](:)
       integer*4    indx[ALLOCATABLE](:,:)
       integer*4    indy[ALLOCATABLE](:,:)
       integer*4    nydata[ALLOCATABLE](:),nnydata[ALLOCATABLE](:,:)
       real*8       t0[ALLOCATABLE](:),tf[ALLOCATABLE](:,:),
     &              states0[ALLOCATABLE](:,:)

calloc       real*8       rnobs(maxnsets)  !******
calloc       integer*4    nobs(maxnsets),ncolxy(maxnsets)
calloc       integer*4    indx(maxnx,maxnsets)
calloc       integer*4    indy(maxnstat,maxnsets)
calloc       integer*4    nydata(maxnsets),nnydata(maxnobs,maxnsets)
calloc       real*8       t0(maxnsets),tf(maxnobs,maxnsets),
calloc     &              states0(maxnstat,maxnsets)

*     Exp design parameters
       integer*4    exppar(maxnopt)
       integer*4    indaux(maxnopt)


*     Broyden impl alg solver
       real*8       smin[ALLOCATABLE](:),smax[ALLOCATABLE](:)
calloc       real*8       smin(maxnstat),smax(maxnstat)

*     Local dimensions for misc tables
      integer*4 nesopsi,nesopsi1,mostobs,mostydat,mostxy,
     &          moststep,nstatea

*     Local dimensions etc for misc solvers
      integer*4 liwddr,lrwddr,nyh,lwm,lrworkod,liwrkod,
     &          meth,mjacob,itol,itask,iopt,mxstep,
     &          isopt,idf,ml,mu,latolod,lrtolod,
     &          lwanr,liwanr,ljacnr

      real*8    srtol,satol,tcrit

      include 'common3.inc'

************** Namelist declarations & default values  ***********************

      namelist /project/  projectname,user,date

      data projectname   /'dummy.txt'/
      data user          /'dummy.txt'/
      data date          /'dummy.txt'/

      namelist /files/    nsets,reportfile,resultfile,optfile,
     &                    datafile,weightfile,modelfile,
     &                    restartfile,dumpfile,ndumpp,evalfile,discfile,
     &                    formats,formatp,formatm,outputf

       data datafile     /maxnfile*'dummy.dat'/
       data weightfile   /maxnfile*'dummy.dat'/
       data modelfile    /'dummy.dat'/
       data reportfile   /'dummy.dat'/
       data resultfile   /'dummy.dat'/
       data evalfile     /'dummy.dat'/
       data discfile     /'dummy.dat'/
       data dumpfile     /'dummy.dat'/
       data ndumpp       /20/
       data optfile      /'dummy.dat'/
       data restartfile  /'dummy.dat'/
       data formats      /'(1x,a,e11.4)'/
       data formatp      /'(1x,e10.3,2x,e10.3,3x,f8.1,5x,f8.1)'/
       data formatm      /'(1x,7e11.3)'/
       data outputf      /'matlab'/

*     The problem specification:

      namelist /problem/  task,model,optimizer,objfun,
     &                    equsolver,odesolver

      data optimizer /'simflex'/
      data objfun    /'dummy'/
      data equsolver /'newton'/
      data odesolver /'odessa'/

*     Model parameters

      namelist /modelpar/nx,nstates,nsaux,
     &                   modelvar,target,estvar,
     &                   ncodir,codirvar,nsets,
     &                   neigdir,eigdir,
     &                   nevastep,ncontour,gridvar,dstep,
     &                   model

      data nsaux       /0/
c      data ncodir     /0/
      data neigdir     /0/
      data nfilepoint  /0/
      data ndiscpoint  /0/
      data nevastep  /maxndir*12/
      data ncontour  /0/
c      data igrid     /maxnopt2*0/
      data dstep     /1.0d-04/
calloc       data smin   /maxnstat*0/
calloc       data smax   /maxnstat*0/

*     The parameters of the data

      namelist /filepar/combined,usew,obss0,nsaveset,
     &                  nobs, ncolxy, indy, indx, nydata !,t0, tf
      data combined  /0/
      data usew      /0/
      data obss0     /0/
      data nsaveset  /1/
calloc       data nobs       /maxnsets*0/
calloc       data ncolxy     /maxnsets*0/
calloc       data t0         /maxnsets*0/
*      data tf         /maxny*uninit/

*     The printing & dataset options:

      namelist /print/echo,optmonit,stats,jacout,sout,echodata,debug
      data echo      /0/
      data echodata  /0/
      data optmonit  /0/
      data debug     /0/
      data stats     /0/
      data jacout    /0/
      data sout      /1/

*     Exp design parameters

      namelist       /design/optcrit, exppar, nexpp
      data optcrit   /'d'/
      data exppar    /maxnopt*0/
      data nexpp     /1/

*     Simflex optimizer

      namelist /simflex/  abstols, reltols, sizes, itmaxs
      data abstols        /1.0d-6/
      data reltols        /1.0d-6/
      data sizes          /0.1/
      data itmaxs         /50/


*     Lev-Mar optimizer

      namelist /levmar/  reltollm,alfalm,itmaxlm,iterdflm
     &
      data reltollm  /1.0d-3/
      data alfalm    /0.0d0/
      data itmaxlm   /0/
      data iterdflm  /100/

c*     Neqn impl alg solver
c
c      namelist /neqnf/    itmaxn, erreln

*     Newton-Raphson impl alg solver

      namelist    /newton/  tolnr,itmaxnr,iprnr,smin,smax
      data tolnr  /1d-14/
      data itmaxnr /200/
      data iprnr   /0/

c*     Ddriv2 ode solver
c
c      namelist /driv2/interpdr2,nrootdr2,epsdr2,ewtdr2,mintdr2
c      data interpdr2 /1/
c      data nrootdr2  /0/
c      data epsdr2    /1.0d-6/
c      data ewtdr2    /0/
c      data mintdr2   /3/

*     Odessa ode solver

      namelist /odessa/  meth, mjacob, mxstep, itol, srtol, satol,
     &                   itask, tcrit,iopt,isopt,idf,ml,mu
      data meth      /2/
      data mjacob    /2/
      data mxstep    /1000/
      data itol      /1/
      data srtol     /1.0d-6/
      data satol     /1.0d-6/
      data itask     /1/
      data tcrit     /0/
      data iopt      /0/
      data isopt     /0/
      data idf       /0/
      data ml        /1/
      data mu        /1/

      data itolod    /2/
      data itaskod   /1/
      data ioptod    /0/

*     Euler ode solver

      namelist /euler/   stpeul
      data stpeul     /1.0d-3/

*     defaults for values computed by 'smbl' from modelvar, target
*     & estvar:

      data btype     /maxnopt*2/
      data bounds    /maxnopt2*0/
      data nsim      /0/
      data lsim      /maxnopt4*0/
      data nest      /0/
      data lest      /maxnopt4*0/
      data nexp      /0/
      data lexp      /maxnopt4*0/
      data nsen      /0/
      data lsen      /maxnopt4*0/
      data nopt      /0/
      data lopt      /maxnopt4*0/
calloc      data s0file    /maxs0*0/

***** defaults for tables not in namelists

calloc      data states0   /maxs0*uninit/

      data target     /' '/
      data estvar     /' '/
      data codirvar   /' '/
      data gridvar    /' '/

      open(unit=2, file = tempin,status='unknown')
      write(2,'(a)') 'stop '
      close(2)

c16-06-94 ****************************************
      do i=1,100
        modelvar(i:i) = ' '
      enddo
c16-06-94 ****************************************

calloc      do k = 1,maxnsets
calloc         nobs(k)   = unini
calloc         ncolxy(k) = unini
calloc         nydata(k) = unini
calloc         t0(k)     = 0.0d0
calloc         do i = 1,maxnx
calloc           indx(i,k) = unini
calloc         enddo
calloc         do i = 1,maxnstat
calloc           indy(i,k)    = unini
calloc           s0file(i,k)  = 0
calloc           states0(i,k) = uninit
calloc         enddo
calloc         do i = 1,maxnobs
calloc           tf(i,k)      = 0.0d0
calloc         enddo
calloc      enddo

************** Initializations ***********************
      idataptr = 1
      rdataptr = 1
      jcount   = 0                        ! to initialize jacob
      sencount = 0                        ! to initialize sensit
************** Read Namelists  ***********************

      call CGETARG(1,nmlfile,ftype)

c     
     
     
     

      read(1,nml=project)
      read(1,nml=files)
      read(1,nml=problem)
      read(1,nml=modelpar)

      if (nsets.lt.1) then
          write(*,*) 'no n of data sets given'
          stop
      endif

      nstatea  = nstates + nsaux
      maxnsets = max(nsets,100) !just to allow more data set definitions
      maxnstat = nstatea

      allocate (ntf(maxnsets),stat=error)
      if (error.ne.0) stop 'Not enough storage for ntf '

      allocate (lpar(maxnlpar,maxnsets),stat=error)
      if (error.ne.0) stop 'Not enough storage for locals'

      allocate (s0file(nstatea,maxnsets),stat=error)
      if (error.ne.0) stop 'Not enough storage for s0 '

      allocate (rnobs(maxnsets),stat=error)
      if (error.ne.0) stop 'Not enough storage for nobs'

      allocate (nobs(maxnsets),stat=error)
      if (error.ne.0) stop 'Not enough storage for nobs'

      allocate (ncolxy(maxnsets),stat=error)
      if (error.ne.0) stop 'Not enough storage for ncolxy'

      allocate (indx(maxnx,maxnsets),stat=error)
      if (error.ne.0) stop 'Not enough storage for indx'

      allocate (indy(maxnstat,maxnsets),stat=error)
      if (error.ne.0) stop 'Not enough storage for indy'

      allocate (nydata(maxnsets),stat=error)
      if (error.ne.0) stop 'Not enough storage for nydata'

      allocate (t0(maxnsets),stat=error)
      if (error.ne.0) stop 'Not enough storage for t0'

      allocate (states0(maxnstat,maxnsets),stat=error)
      if (error.ne.0) stop 'Not enough storage for states0'

      allocate (smin(maxnstat),stat=error)
      if (error.ne.0) stop 'Not enough storage for smin'

      allocate (smax(maxnstat),stat=error)
      if (error.ne.0) stop 'Not enough storage for smax'

      do k = 1,maxnsets
         nobs(k)   = unini
         nydata(k) = unini
         ncolxy(k) = unini
         t0(k)     = 0.0d0
         do i = 1,maxnx
           indx(i,k) = unini
         enddo
         do i = 1,maxnstat
           indy(i,k)    = unini
           s0file(i,k)  = 0
           states0(i,k) = uninit
         enddo
      enddo
      do i = 1,maxnstat
         smin(i) = 0.0d0
         smax(i) = 0.0d0
      enddo

      read(1,nml=filepar)

      if (nobs(1).lt.1) then
          write(*,*) 'no n of observations given for 1. data set'
          stop
      endif

      if (nsets.gt.1) then
       do k = 2,nsets
         if (nobs(k).eq.unini)   nobs(k)   = nobs(1)
         if (nydata(k).eq.unini) nydata(k) = nydata(1)
       enddo
      endif

       do k=1,nsets
         maxnobs = max(maxnobs,nobs(k))
       enddo
       maxnobs   = maxnobs+1
       maxny     = maxnobs*maxnsets
       maxnnew   = maxny
       allocate (nnydata(maxnobs,maxnsets),stat=error)
       if (error.ne.0) stop 'Not enough storage for nnydata'

       allocate (tf(maxnobs,maxnsets),stat=error)
       if (error.ne.0) stop 'Not enough storage for tf'

       allocate (newobs(maxnnew,2),stat=error)
       if (error.ne.0) stop 'Not enough storage for newobs '

      do k = 1,nsets
      do i = 1,maxnobs
        nnydata(i,k)  = unini
        tf(i,k)      = 0.0d0
      enddo
      enddo

      read(1,nml=print)
      read(1,nml=design)
      read(1,nml=simflex)
      read(1,nml=levmar)
c      read(1,nml=neqnf)
      read(1,nml=newton)
c      read(1,nml=driv2)
      read(1,nml=odessa)
      read(1,nml=euler)

      close(1)


      if (odesolver(1:4) .eq. 'odes' .or.
     &     odesolver(1:4) .eq. 'lsod') then
                                        
         if (meth.eq.1 .or. meth.eq.2 .and.
     &        mjacob.ge.0 .and. mjacob.le.5) then
              mfod = 10*meth + mjacob
         endif
         if (isopt.eq.1 .and. mjacob.eq.0 .or. mjacob.eq.3) then
            write(*,*) 'Odessa variables isopt & mjacob incompatible'
            stop
         endif

          method    = meth
          miterod   = mjacob
          itolod    = itol
          itaskod   = itask
          ioptod    = iopt
          isoptod   = isopt
          idfod     = idf
          mlod      = ml
          muod      = mu
          mxstepod  = mxstep
          srtolod   = srtol
          satolod   = satol
          tcritod   = tcrit

c         do i = 1,nstates*(nest+1)
c            rtolod(i) = srtol  ! muutosta!!!!
c            atolod(i) = satol
c         enddo

      endif

        call smbl(task,model,modelvar,target,estvar,nsets,
     *            gpar,lpar,states0,s0file,lsim,lest,lexp,lsen,
     *            lopt,bounds,nbounds,ngpar,nlpar,nstates0,nsim,nest,
     *            nexp,nsen,nopt,nmlfile,error_number,maxnstat,
     *            maxnlpar,mxchars,mxwords,mxmsent,maxnopt,type,
     *            charcl,nglobal,nlocal,ncontrol,ninitval,laux,
     *            word,sentence,name,ncodir,codir,
     *            evalfile,discfile,nfilepoint,ndiscpoint,
     *            ngrid,igrid,codirvar,gridvar,maxnobs,
     *            t0,tf,ntf,vardim,nx,nsbound,sbounds,optcrit,
     *            objfun,projectname,mustcmp,ftype,inttype,iperiod)

      if(error_number.eq.8) then
         write(*,*) 'Something wrong in namelist modelpar'
         stop
      endif

      ncontour = ngrid


************** defaults for data parameters incompletely given in nml****

************** output file names *********************************

c     result/report/optfile according to projectnme
c     dumpfile              according to resultfile

      if (resultfile .eq.   'dummy.dat' ) then
          resultfile(1:9) = '         '
          if (task(1:3).eq.'sim') then
             resultfile = projectname(1:iperiod) //'.sim'
          elseif (task(1:3).eq.'est') then
             resultfile = projectname(1:iperiod) //'.est'
          elseif (task(1:3).eq.'exp') then
             resultfile = projectname(1:iperiod) //'.exp'
          elseif (task(1:3).eq.'sen') then
             resultfile = projectname(1:iperiod) //'.sen'
          elseif (task(1:3).eq.'opt') then
             resultfile = projectname(1:iperiod) //'.opt'
          endif
      endif

      if (reportfile(1:9) .eq.   'dummy.dat'  .and.
     &    task(1:3)  .eq.   'est' .and.
     &    stats      .eq.    1) then
          reportfile(1:9) = '         '
          reportfile = projectname(1:iperiod) //'.sta'
      endif

      if (optfile(1:9) .eq.   'dummy.dat' .and.
     &     (task(1:3).eq.'est' .or.
     &     task(1:3).eq.'exp' .or.
     &     task(1:3).eq.'opt') ) then
           optfile(1:9) = '         '
           optfile = projectname(1:iperiod) //'.arg'
      endif

      if(dumpfile(1:9) .eq. 'dummy.dat') then
         j = index(resultfile,'.')
         if(j .gt. 0) then
             dumpfile(1:j-1)   = resultfile(1:j-1)
             dumpfile(j:j+3) = '.dmp'
             do k = j+4,len(dumpfile)
                 dumpfile(k:k) = ' '
             enddo
         else
             j                 = index(resultfile,' ')
             dumpfile(1:j-1)   = resultfile(1:j-1)
             dumpfile(j:j+3)   = '.dmp'
             do k = j+4,len(dumpfile)
                 dumpfile(k:k) = ' '
             enddo
         endif
      endif

************** weight files  ************************************

      if (usew .eq.1 .and. combined.eq.0) then

         if (weightfile(1) .eq. 'dummy.dat') then
            write(*,*) 'No weight file for data set 1'
            write(*,*) ' '
            stop
         endif
         do k=2,nsets
            if (weightfile(k) .eq. 'dummy.dat') then
                weightfile(k) = weightfile(1)
            endif
         enddo

      endif

************** nnydata,nobs,ncolxy,indx,indy,s0,tf ***************************
      do k=1,nsets
         nnydata(1,k) = nydata(k)
      enddo
******change: remove obsolete nnydata!! ***************

      k = 1
      if (nnydata(1,1).gt.0) then
         do j = 2,nobs(k)
            if (nnydata(j,k).eq.unini) nnydata(j,k) = nnydata(1,k)
         enddo
      endif
      do i = 1,nstates
         if (states0(i,k).eq.uninit) states0(i,k) = 0.0d0
      enddo

      if (nsets.gt.1) then
       do k = 2,nsets

         if (nobs(k).eq.unini)   nobs(k)   = nobs(1)
         if (ncolxy(k).eq.unini) ncolxy(k) = ncolxy(1)

         do i = 1,nstates
            if (states0(i,k).eq.uninit) states0(i,k) = states0(i,1)
         enddo

         if (nnydata(1,k).eq.unini) nnydata(1,k) = nnydata(1,1)
         if (tf(1,k).eq.uninit) tf(1,k)     = tf(1,1)
         do j = 2,nobs(k)
            if (nnydata(j,k).eq.unini) nnydata(j,k) = nnydata(1,k)
            if (tf(j,k).eq.uninit) tf(j,k)     = tf(j,1)
            do i = 1,nx
               if (indx(i,k).eq.unini) indx(i,k)  = indx(i,1)
            enddo
            do i = 1,nnydata(j,k)
               if (indy(i,k).eq.unini) indy(i,k)  = indy(i,1)
            enddo
         enddo

       enddo
      endif

      ntot = 0
      do k = 1,nsets
      do j = 1,nobs(k)
         ntot = ntot + nnydata(j,k)
      enddo
      enddo

****** the objective function ****************************************
      objunini = (objfun(1:5) .eq. 'dummy')
      if (task(1:3).eq.'est' .and. objunini) then
         objfun(1:3) = 'lsq'
      elseif (task(1:3).eq.'exp' .and. objunini) then
         objfun(1:3) = 'exp'
      elseif (task(1:3).eq.'exp' .and. optcrit(1:6).eq.'glolsq'
     &                           .and. objunini) then
         objfun(1:3) = 'exp'
      elseif (task(1:3).eq.'opt' .and. objunini) then
         objfun(1:3) = 'opt'
      elseif (task(1:3).eq.'sen' .and. objunini) then
         objfun(1:3) = 'lsq'
      endif

***** dimensions and defaults for simulations & sensitivity ******************

      neigdir = 0
      do i = 1,maxnopt
         if (eigdir(i).ne.0) neigdir = neigdir+1
      enddo
      
      if (task(1:3).eq.'exp') then
       if (optcrit(1:6) .ne. 'glolsq')  jacout = 1
         if (optcrit(1:3).eq.'glo' .or.
     &       optcrit(1:3).eq.'GLO') then
             nsen = nest
         endif
      endif

      nesopsi  = max(nsim,nest,nexp,nsen,nopt)    ! dimension for various arrays
      nesopsi1 = nesopsi+1                   ! ??
      ncoaux   = max(nsim,nsen)
      ndir     = ncodir + neigdir            ! colines, eiglines, file below

      if (nfilepoint.gt.0) then
         ndir = ndir + 1                 ! file'line'
         nevastep(ndir) = nfilepoint-1   ! n of 'steps' = n of points in file
      endif

******nsim or nsen > 0, no dir/filepoints given: default 'all 1 by 1'
      if (ncoaux.gt.0 .and. ndir.eq.0) then
          ncodir = ncoaux
          do i = 1,ncodir
             codir(i) = i
          enddo
          ndir   = ncodir
      endif

      if (neigdir.gt.0) then ! nest,lest needed to calc eig.vectors
          if (task(1:3).eq.'sim') then
             nest = nsim
             do i = 1,4
             do j = 1,nsim
                lest(i,j)=lsim(i,j)
             enddo
             enddo
          elseif (task(1:3).eq.'sen') then
             nest = nsen
             do i = 1,4
             do j = 1,nsen
                lest(i,j)=lsen(i,j)
             enddo
             enddo
          endif
      endif

      if (task(1:3).eq.'sen' .and. objfun(1:3).eq.'exp') then
          nexp = nsen
          do i = 1,4
          do j = 1,nsen
             lexp(i,j)=lsen(i,j)
          enddo
          enddo
      endif

      if (task(1:3).eq.'sen' .and. objfun(1:3).eq.'opt') then
          nopt = nsen
          do i = 1,4
          do j = 1,nopt
             lopt(i,j)=lsen(i,j)
          enddo
          enddo
      endif

      moststep = 0
      if (ndir .gt. 0) then
         do i=1,ndir
            moststep=max(moststep, nevastep(i))
         end do
      end if
      moststep = moststep+1               ! one extra place for 'x0' point
      moststep = max(moststep,nfilepoint) ! note n of points in file'line'

***** dimensions for simulations & sensitivity done ******************

*     substitute in the matrix 'newobs' the indexes of the new
*     observations to be optimized (in case of task 'exp' or 'opt')

      if (task(1:3).eq.'exp' .or.
     &  (task(1:3).eq.'sen' .and. objfun(1:3).eq.'exp')) then

        call newexp(lexp,nexp,nobs,nsets,newobs,maxnnew,nnew)
        if (optcrit(1:1).eq.'s') call subset(exppar,nest,nexpp,indaux)

      elseif (task(1:3).eq.'opt' .or.
     &  (task(1:3).eq.'sen' .and. objfun(1:3).eq.'opt')) then

        call newexp(lopt,nopt,nobs,nsets,newobs,maxnnew,nnew)

      endif

************** Actual max dimensions of the large arrays ***************

c     some computations

      nstatea = nstates + nsaux             ! actual and auxiliary states

      mostobs=0
      do i=1,nsets
         mostobs=max(mostobs, nobs(i))      ! max n of obss in a set
      end do
      ifile = 0
      do i = 1,nstatea
      do k = 1,nsets
         ifile = ifile + s0file(i,k)
      enddo
      enddo
c      if (model(1:3) .eq. 'ode' .and. ifile.eq.0) then
         mostobs = mostobs+1                   ! make room for initial values
c      endif

      mostydat=1
      do k=1,nsets
         do j=1,nobs(k)
           mostydat=max(mostydat, nnydata(j,k)) ! max n of ydata in an obs
         end do
      end do

      mostxy=1
      do k=1,nsets
           mostxy=max(mostxy, ncolxy(k)) ! max n of col
      end do
********************just for debug *******

c      include 'debug.hlp'


************** Set actual values for dimension names ***************
************** unspecified: 0's

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
      dngrid   =  ngrid
      dnevapoi =  nfilepoint
      ddiscpoi =  ndiscpoint
      dncontou =  ncontour
      dliwddr  =  liwddr       ! values in modio
      dlrwddr  =  lrwddr
      dlrwrkod =  lrworkod
      dliwrkod =  liwrkod
      dlrtolod =  lrtolod
      dlatolod =  latolod

      dwanr    =  lwanr
      diwanr   =  liwanr
      djacnr   =  ljacnr
c      dlwadne  =  lwadne

*********** write the dimensions, scalars and characters *******
*********** to temporary io file ******
      open(unit=1, file   = tempin,status = 'unknown')

      if (mustcmp.eq.0) then
         write(1,'(a)') 'ready'
      elseif (mustcmp.eq.1) then
         write(1,'(a)') 'changed'
      elseif (mustcmp.eq.2) then
         write(1,'(a)') 'missing'
      endif

      do i = 1,ndim
         write(1,*) dim(i)
      enddo
      do i = 1,nscali
         write(1,*) scali(i)
      enddo
      do i = 1,nscalr
         write(1,*) scalr(i)
      enddo
      do i = 1,nchar
        write(1,'(a)') char(i)
      enddo
      if (combined.eq.1) then
        write(1,'(a)') datafile(1)
      elseif (combined.eq.0) then
        write(1,'(a)') datafile(1)
        if (nsets.gt.maxnfile .and. datafile(1).ne.'dummy.dat') then
           write(*,*) 'max ', maxnfile, '  separate datafiles allowed'
           write(*,*) 'for more data sets:  combine in one file'
           write(*,*) ' '
           stop
        endif

        if (datafile(1).ne.'dummy.dat') then
          do i = 2,nsets
            write(1,'(a)') datafile(i)
          enddo
          do i = 1,nsets
            write(1,'(a)') weightfile(i)
          enddo
        endif
      endif


************ write the tables to temporary io file *****************************

c    states etc:
      call wralloc(nstatea,nsets,1,1,
     &                     states0,maxnstat,maxnsets,1,1)
      call wialloc(nstatea,nsets,1,1,
     &                     s0file,maxnstat,maxnsets,1,1)
      call wralloc(nsets,1,1,1,t0,maxnsets,1,1,1)
      call wralloc(mostobs,nsets,1,1,
     &                     tf,maxnobs,maxnsets,1,1)

c    pointed variables and pointers:
      call wralloc(ngpar,1,1,1,gpar,maxngpar,1,1,1)
      call wralloc(nlpar,nsets,1,1,
     &                    lpar,maxnlpar,maxnsets,1,1)
      call wialloc(4,nsim,1,1,lsim,4,maxnopt,1,1)
      call wialloc(4,nest,1,1,lest,4,maxnopt,1,1)
      call wialloc(4,nexp,1,1,lexp,4,maxnopt,1,1)
      call wialloc(4,nsen,1,1,lsen,4,maxnopt,1,1)
      call wialloc(4,nopt,1,1,lopt,4,maxnopt,1,1)

c    generally used vectors:
      call wralloc(2,nesopsi,1,1,bounds,2,maxnopt,1,1)
      call wralloc(2,nesopsi,1,1,sbounds,2,maxnopt,1,1)
      call wialloc(nesopsi,1,1,1,btype,maxnopt,1,1,1)

c    data indexes:
*********************
      do i = 1,nsets
         rnobs(i) = dble(nobs(i))
      enddo
      call wralloc(nsets,1,1,1,rnobs,maxnsets,1,1,1)
*********************
      call wialloc(nsets,1,1,1,nobs,maxnsets,1,1,1)
      call wialloc(mostobs,nsets,1,1,
     &                     nnydata,maxnobs,maxnsets,1,1)
      call wialloc(nsets,1,1,1,ncolxy,maxnsets,1,1,1)
      call wialloc(nx,nsets,1,1,indx,maxnx,maxnsets,1,1)
      call wialloc(mostydat,nsets,1,1,
     &                     indy,maxnstat,maxnsets,1,1)

      if (task(1:3) .eq. 'exp' .or.
     &   (task(1:3) .eq. 'sen' .and. objfun(1:3).eq.'exp')) then

        call wialloc(nexpp,1,1,1,exppar,maxnopt,1,1,1)
        call wialloc(nnew,2,1,1,newobs,maxnnew,2,1,1)

      endif

      call wialloc(ndir,1,1,1,nevastep,maxndir,1,1,1)
      call wialloc(ncodir,1,1,1,codir,maxnopt,1,1,1)
      call wialloc(neigdir,1,1,1,eigdir,maxnopt,1,1,1)
      call wialloc(2,ngrid,1,1,igrid,2,maxnopt,1,1)

c      if (odesolver(1:4).eq.'odes') then
c          rtol,atol
c      endif

      if (model(1:3).eq.'imp' .and. equsolver(1:3).eq.'new') then
        call wralloc(nstates,1,1,1,smin,maxnstat,1,1,1)
        call wralloc(nstates,1,1,1,smax,maxnstat,1,1,1)
      endif

       close(1)

************** io savings DONE ************************
      end

      subroutine wralloc(dim1,dim2,dim3,dim4,
     &                   table,ld1,ld2,ld3,ld4)
c     This subroutine saves the real arrays used in modest program
c           dim's        the actual dimensions (default=1)
c           table        values saved
c           ld's         leading dimensions of 'table', as declared in
c                        the calling subprogram (default=1)

      integer*4 dim1,dim2,dim3,dim4,dd
      integer*4 ld1,ld2,ld3,ld4
      real*8    table(ld1,ld2,ld3,ld4)

      integer*4 i,j,k,l
      include 'common3.inc'

      dd = dim1*dim2*dim3*dim4
      if (dd.eq.0) then
          return
      endif

         do l=0,dim4-1
         do k=0,dim3-1
         do j=0,dim2-1
         do i=0,dim1-1
         write(1,*) table(i+1,j+1,k+1,l+1)
         enddo
         enddo
         enddo
         enddo

      return
      end

      subroutine wialloc(dim1,dim2,dim3,dim4,
     &                   table,ld1,ld2,ld3,ld4)
c     This subroutine saves the integer arrays used in modest program
c           dim's        the actual dimensions (default=1)
c           table        values saved
c           ld's         leading dimensions of 'table', as declared in
c                        the calling subprogram (default=1)

      integer*4 dim1,dim2,dim3,dim4,dd
      integer*4 ld1,ld2,ld3,ld4
      integer*4 table(ld1,ld2,ld3,ld4)

      integer*4 i,j,k,l
      include 'common3.inc'

      dd = dim1*dim2*dim3*dim4
      if (dd.eq.0) then
         return
      endif

         do l=0,dim4-1
         do k=0,dim3-1
         do j=0,dim2-1
         do i=0,dim1-1
             write(1,*) table(i+1,j+1,k+1,l+1)
         enddo
         enddo
         enddo
         enddo

      return
      end




      subroutine newexp(laux,naux,
     &                  nobs,nsets,
     &                  newobs,maxnnew,nnew)

*     substitutes in the matrix 'newobs' the indexes of the new
*     observations to be optimized (in case of task 'exp' or 'opt')

      implicit none

      integer*4 naux,nsets,maxnnew,nnew
      integer*4 laux(4,naux)
      integer*4 nobs(nsets),newobs(maxnnew,2)

      include 'common3.inc'

*     LOCAL VARIABLES

      integer*4      count,all,type,prevk
      integer*4      i,j,k,ii,jj,jj0      ! loop indices

*     THE ALGORITHM

      count = 0
      prevk = -1
      i     = 1

      do while (i .le. naux)

          type = laux(1,i)

          if (type.eq.1) then
              do k = 1,nsets
                 do j = 1,nobs(k)
                    ii = count + j
                    newobs(ii,1) = j
                    newobs(ii,2) = k
                 enddo
                 count = count + nobs(k)
              enddo
              return
          endif

          k   = laux(4,i)

          if (model.eq.'ode') then

             if (k.ne.prevk) then
                all   = 1
                i     = i + 1
                prevk = k
             else
                i   = i + 1
                all = 0
             endif

          elseif (model.eq.'alg' .or. model.eq.'imp') then

             all = 0
             jj  = 0
             jj0 = -1
             do while ((i.le.naux) .and. (k.eq.laux(4,i)))
                 type = laux(1,i)
                 if (type.eq.2 .or. type.eq.4) then
                     all = 1
                 elseif (type.eq.3) then
                     jj = laux(3,i)
                     if (jj.ne.jj0) then
                        jj0 = jj
                        count = count + 1
                        newobs(count,1) = jj
                        newobs(count,2) = k
                     endif
                 else
                     write(*,*) 'the type in lopt must be 1/2/3/4'
                     stop
                 endif
                 i = i + 1
             enddo

          endif

          if (all.eq.1) then
              do j = 1,nobs(k)
                  ii = count + j
                  newobs(ii,1) = j
                  newobs(ii,2) = k
              enddo
              count = count + nobs(k)
          endif

      enddo

      nnew = count

      return
      end


      subroutine subset(exppar,nest,nexpp,iaux)

*     initially, the vector 'exppar' contains the parameters of
*     interest for subset exp.design. In the 'subset' case 'expdes'
*     uses the complementary parameters. The content of 'exppar' is
*     switched here. 'nexpp' gives the n. of nuissance parameters.

      implicit none
*     arguments
      integer*4  nest,nexpp,iaux(nest),exppar(nest)

*     local variables
      integer*4  i

      i = 1
      do while (exppar(i) .gt.0)
          iaux(exppar(i)) = exppar(i)
          i = i+1
      enddo

      nexpp = 0
      do i = 1,nest
          if (i .ne. iaux(i)) then
              nexpp = nexpp + 1
              exppar(nexpp) = i
          endif
      enddo

      return
      end


       subroutine cgetarg(argno, argbuf,ftype)
       integer*4     argno
       character*80  argbuf
       character*3   ftype

c     CONVEX,SUN,STARDENT:
c       call getarg(argno,argbuf)
c       ftype = 'f  '
c     HP:
c       call igetarg(argno, argbuf, 80)
c       ftype = 'f  '
*     PC:
      integer*2     argstat
      ftype = 'for'
      call getarg(argno, argbuf, argstat)

*     VAX:
*      integer       cli$get_value, status, length
*      character*(*) arglist
*      parameter(arglist='P1P2P3P4P5P6P7P8P9')
*
*      ftype = 'for'
c      status = cli$get_value(arglist(argno*2-1:argno*2),argbuf,length)
c      if (.not. status) call lib$signal(%val(status))

       return
       end




