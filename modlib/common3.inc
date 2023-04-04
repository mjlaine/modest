
      integer*4 arrays,rdataptr,idataptr,nscali,nscalr,nchar,
     &          ndim,maxnfile
      parameter (arrays   = 100) ! # of arrays
      parameter (nscali   = 50,nscalr=50,ndim=50)
      parameter (maxnfile = 100)
      parameter (nchar    = 22)

      real*8    scalr(nscalr)
      integer*4 scali(nscali)
      integer*4 arr(arrays)
      integer*4 dim(ndim)
      character*80 char(nchar)
      common   /rdata/rdataptr
      common   /idata/idataptr

****** All base pointers ****************************
      integer*4 pyod,psdr2,psimp,pstates,pstates0,ps0file,
     &          pstatesa,pt0,ptf,pyest,pxjac,pgpar,plpar,
     &          plsim,plest,plexp,plsen,plopt,
     &          pxsim,pxest,pxexp,pxsen,pxopt,pxaux,pbounds,pbtype,pjtj,
     &          psbounds,prnobs,
     &          pnobs,pnydata,pncolxy,pindx,pindy,pymin,pymax,pfvecdne,
     &          pwabro,pwadne,pexpdu,pexpdv,pexppar,pjtj0,pjtj1,pjtj2,
     &          psigma,ppenpar,pnewobs,pjacxorg,pcoline,pigrid,peigdir,
     &          pcodir,pnsimste,peigval,peigvec,
     &          poptline,poptgrid,plsqxaux,pvaldisc,
     &          psimxaux,ptlaux,ptuaux,pdir,pxdir,pxdisc,
     &          psimplex,pwrksimf,piwrsimf,
     &          pzlm,pipivlm,
     &          plaux,piauxsub,piwrkddr,pworkddr,pesvd,pwrksvd,
     &          pwrk2svd,pwrkdsco,
     &          prtolod,patolod,piworkod,prworkod,pteta0od,
     &          pworkeul,
     &          psmin,psmax,pwanr,piparnr,piwanr,
     &          pxdata,pydata,pweight,pxyaux

****** Dimension names in common 'dim' ******************************

      integer*4 dnstates,dnx,dmostobs,dnsavese,dnsaux,dnstatea,
     &          dnsets,dmostyda,dmostxy,dngpar,dnlpar,
     &          dnsim,dnest,dnexp,dnsen,dnopt,
     &          dntot,dnnew,dnexpp,dnesopsi,dmostste,
     &          dndir,dncodir,dneigdir,dngrid,dnevapoi,dncontou,
     &          ddiscpoi,
     &          dliwddr,dlrwddr,
     &          dlatolod,dliwrkod,dlrtolod,dlrwrkod,
     &          dwanr,diwanr,djacnr

****** temporary files ***********************************************

      character*80 tempin, tempsav, tempout, tempaux

      parameter (tempin   = 'modtemp.in      ')
      parameter (tempsav  = 'modtemp.sav     ')
      parameter (tempout  = 'modtemp.out     ')
      parameter (tempaux  = 'modtemp.aux     ')

********** Namelist variables, scalars here in commons, ***************
********** allocatable tables in 'setdim', dimensions in 'arr' ********

*     Project,file names, formats, output form
       character*80  projectname,user,date
       character*80  datafile(maxnfile),weightfile(maxnfile)
       character*80  reportfile,resultfile,optfile,modelfile
       character*80  restartfile,dumpfile,evalfile,discfile
       character*80  formats,formatp,formatm
       character*80  outputf

*     The problem specification
       character*80 task,model,optimizer,objfun,equsolver,odesolver
       integer*4    exjtj
       integer*4    ndumpp,dumpcall

*     Parameters of data
       integer*4    combined,usew,obss0 !,nsaveset

*     The printing & dataset options:
       integer*4 echo,optmonit,stats,jacout,sout,echodata,debug

*     Exp design parameters
       integer*4    jcount,sencount
       character*80 optcrit

*     Jacobian
       real*8       dstep

*     Simflex optimizer
       real*8       abstols,reltols,sizes
       integer*4    itmaxs

*     Lmdiff optimizer
       real*8       reltollm, alfalm
       integer*4    itmaxlm,iterdflm

*     Neqn impl alg solver
       real*8       erreln
       integer*4    itmaxn

*     Broyden impl alg solver

       real*8       tolnr,resnr
       integer*4    itmaxnr,iprnr,itcntnr

*     Ddriv2 ode solver

       integer*4    interpdr2,nrootdr2,mintdr2
       real*8       epsdr2,ewtdr2

*     Odessa ode solver

       integer*4    itolod,itaskod,ioptod,isoptod,idfod
       integer*4    mlod,muod,mxstepod
       real*8       srtolod, satolod, tcritod
       integer*4    method,miterod,mfod

*     Euler ode solver

       real*8       stpeul

***** Miscellaneous var *******************************

      integer*4 iobs,iset

******************************************************

      equivalence(arr,pyod)
      equivalence(dim,dnstates)
      equivalence(scali,exjtj)
      equivalence(scalr,dstep)
      equivalence(char,projectname)

****** All commons ************************************

      common /arr/
     &       pyod,psdr2,psimp,pstates,pstates0,ps0file,
     &       pstatesa,pt0,ptf,pyest,pxjac,pgpar,plpar,
     &       plsim,plest,plexp,plsen,plopt,
     &       pxsim,pxest,pxexp,pxsen,pxopt,pxaux,pbounds,pbtype,pjtj,
     &       psbounds,prnobs,
     &       pnobs,pnydata,pncolxy,pindx,pindy,pymin,pymax,pfvecdne,
     &       pwabro,pwadne,pexpdu,pexpdv,pexppar,pjtj0,pjtj1,pjtj2,
     &       psigma,ppenpar,pnewobs,pjacxorg,pcoline,pigrid,peigdir,
     &       pcodir,pnsimste,peigval,peigvec,
     &       poptline,poptgrid,plsqxaux,pvaldisc,
     &       psimxaux,ptlaux,ptuaux,pdir,pxdir,pxdisc,
     &       psimplex,pwrksimf,piwrsimf,
     &       pzlm,pipivlm,
     &       plaux,piauxsub,piwrkddr,pworkddr,pesvd,pwrksvd,
     &       pwrk2svd,pwrkdsco,
     &       prtolod,patolod,piworkod,prworkod,pteta0od,
     &       pworkeul,
     &       psmin,psmax,pwanr,piparnr,piwanr,
     &       pxdata,pydata,pweight,pxyaux

      common /dim/
     &       dnstates,dnx,dmostobs,dnsavese,dnsaux,dnstatea,
     &       dnsets,dmostyda,dmostxy,dngpar,dnlpar,
     &       dnsim,dnest,dnexp,dnsen,dnopt,
     &       dntot,dnnew,dnexpp,dnesopsi,dmostste,
     &       dndir,dncodir,dneigdir,dngrid,dnevapoi,dncontou,
     &       ddiscpoi,
     &       dliwddr,dlrwddr,
     &       dlatolod,dliwrkod,dlrtolod,dlrwrkod,
     &       dwanr,diwanr,djacnr

      common /char/  projectname,user,date,
     &               restartfile,dumpfile,modelfile,resultfile,
     &               evalfile,discfile,reportfile,optfile,
     &               formats,formatp,formatm,outputf,
     &               task,model,optimizer,odesolver,equsolver,optcrit,
     &               objfun

      common /dfile/ datafile,weightfile

      common /scali/ exjtj,jcount,sencount,ndumpp,
     &               combined,usew,obss0,dumpcall,
     &               echo,optmonit,stats,jacout,sout,echodata,debug,
     &               itmaxs, itmaxlm,iterdflm,
     &               itmaxn,
     &               itmaxnr,iprnr,itcntnr,
     &               interpdr2,nrootdr2,mintdr2,
     &               itolod,itaskod,ioptod,isoptod,idfod,
     &               mlod,muod,mxstepod,method,miterod,mfod

      common /scalr/ dstep,
     &               abstols,reltols,sizes,
     &               reltollm, alfalm,
     &               erreln,
     &               tolnr,resnr,
     &               epsdr2,ewtdr2,
     &               srtolod, satolod, tcritod,
     &               stpeul

      common /misc/  iobs,iset

