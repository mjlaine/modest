C $Id: mdstlsq2.f,v 1.6 2010/01/18 15:38:48 mjlaine Exp $
C ------------------------------------------------------------------------
C mdstmcmc library 
C File: mdstlsq2.f
C Purpose: new lsq target for modest that uses mcmclib -2*log(likelihood)
C          and prior, parameters are read from mcmc namelist
C
C Marko Laine <Marko.Laine@Helsinki.fi>
C ------------------------------------------------------------------------
C $Author: mjlaine $
C $Revision: 1.6 $
C $Date: 2010/01/18 15:38:48 $
C ------------------------------------------------------------------------
C $Log: mdstlsq2.f,v $
C Revision 1.6  2010/01/18 15:38:48  mjlaine
C mcmccovn.dat file processing and extra file io routines
C
C Revision 1.5  2010/01/11 14:44:31  mjlaine
C utils: file locking, mcmc: cov files parametrized
C
C Revision 1.4  2009/02/10 11:38:11  mjlaine
C precision checks on constants
C
C Revision 1.3  2009/01/20 12:43:15  mjlaine
C mdstlsq2.f is used (I forgot this)
C
C Revision 1.2  2007/09/19 07:00:22  mjlaine
C gfortran compile fix
C
C Revision 1.1  2007/09/18 18:19:28  mjlaine
C added function to be used as modest optimization target
C
C ------------------------------------------------------------------------
C
C F77 version
      function mdstlsq2(ar,nr,ai,ni,xaux,naux)

c no do, mcmc has not been inited yet
c      use mcmc, only: sstype, sstrans 

      implicit none
      real*8 mdstlsq2
      integer*4      nr,ni
      real*8         ar(nr)
      integer*4      ai(ni)
      integer*4      naux           
      real*8         xaux(naux+1) 
	
      integer*4 ny
      integer*4, save :: sstype
      real*8, save :: sstrans
      real*8, save, allocatable :: mu(:), sig(:)
      logical, save :: inited = .false.

      real*8 prior2lfun
      external prior2lfun

      interface
         function mdstlsqs(ar,nr,ai,ni,xaux,naux,ny,sstype,sstrans)
         implicit none
         integer*4 ny
         real*8         mdstlsqs(ny)
         integer*4      nr,ni
         real*8         ar(nr)
         integer*4      ai(ni)
         integer*4      naux
         real*8         xaux(naux+1)
         integer*4      sstype
         real*8         sstrans
         end function mdstlsqs
      end interface

      include 'common3.inc'
            
      ny = ai(pnydata)

c need to get these from the namelist
      if (.not.inited) then
         allocate(mu(naux),sig(naux))
         call getssparams(sstype,sstrans,mu,sig,naux)
         inited = .true.
      end if

      mdstlsq2 = sum(mdstlsqs(ar,nr,ai,ni,xaux,naux,ny,sstype,sstrans))
c add prior "sum of squares"
      mdstlsq2 = mdstlsq2 + prior2lfun(xaux,mu,sig,naux)

      return
      end
c-------------------------------------------------------------------------
c
c we need some parameters from the mcmc namelist
c
      subroutine getssparams(sstype,sstrans,mu,sig,len)
      implicit none
      integer*4 sstype
      real*8 sstrans
      integer*4 len
      real*8 mu(len), sig(len)

      integer, parameter :: dbl=8
c we need to define the whole mcmc namelist here
c this is copied from mcmc.F90

cc mcmc namelist variables
      integer, save :: nsimu    ! length of the chain
      integer, save :: doadapt  ! do we adapt?
      integer, save :: doburnin ! do we have burn in time?
      integer, save :: adaptint ! adaptation interval
      integer, save :: adapthist ! history size
      integer, save :: initcmatn ! n for initial cmat
      integer, save :: burnintime ! time to burn in
      integer, save :: badaptint ! burn-in adapt interval
      integer, save :: adaptend ! end adaptation
      integer, save :: greedy   ! greedy burn in adaptation
      real(kind=dbl), save :: scalelimit ! when to scale during burn in
      real(kind=dbl), save :: scalefactor ! scale this much
      real(kind=dbl), save :: drscale ! scale in dr
      integer, save :: svddim   ! not used, use condmax instead
      real(kind=dbl), save :: condmax ! use svd with condmax
      real(kind=dbl), save :: condmaxini ! modest jtj reqularization
      real(kind=dbl), save :: N0, S02 ! error variance prior
c      integer, save :: sstype   ! type of the likelihood
c      real(kind=dbl), save :: sstrans ! transform in mdstlsqs
      integer, save :: filepars ! read initial values from files
      integer, save :: printint ! how often to give statistics
      integer, save :: updatesigma ! do we update sigma2
      integer, save :: usrfunlen ! how many extra values to calculate?
      integer, save :: dumpint  ! dump interval
      character(len=128), save :: chainfile ! name of the chain file
      character(len=128), save :: s2file ! name of the s2chain file
      character(len=128), save :: ssfile ! name of the sschain file
      character(len=128), save :: priorsfile ! name of the priors file
      character(len=128), save :: cov0file ! name of the initial covmat file
      character(len=128), save :: covffile ! name of the final covmat file
      character(len=256), save :: covnfile ! file to hold the weight of covmat
      integer, save :: verbosity ! how much to print
      character(len=10), save :: method ! the method used (dram, scam)

c      namelist /lsq/ sstype, sstrans

      namelist /mcmc/ nsimu, doadapt, doburnin, adaptint, adapthist, 
     b   scalelimit, scalefactor, drscale, 
     b   badaptint, adaptend, initcmatn, N0, S02, 
     b   filepars, burnintime, greedy, printint, updatesigma, 
     b  usrfunlen, chainfile, s2file, ssfile, svddim, condmax, 
     b  cov0file, covffile, covnfile,
     b  condmaxini, sstype, sstrans, dumpint, priorsfile, verbosity, 
     b  method



      logical, save :: inited=.false.
      integer fstat,i
      logical fexist
      character(len=128) :: nmlfile 


      if (.not. inited) then
         sstrans = 1.0d0
         sstype = 1
         nmlfile = ''
         priorsfile = ''
c read the modest namelist file name from mcmcnml.txt
         inquire(file='mcmcinit.nml',exist=fexist)
         if (fexist) then
            nmlfile = 'mcmcinit.nml'
         else
            inquire(file='mcmcnml.txt',exist=fexist)
            if (fexist) then
               open(unit=10, file='mcmcnml.txt', status='old', 
     b              iostat=fstat)
               read(10,*) nmlfile
               close(10)
               nmlfile=trim(nmlfile)
            end if
         end if
         if (len_trim(nmlfile) == 0) then
            write(*,*) 'namelist file not found:',nmlfile
            inited = .true.
            return
         end if
!!!   read some parameters from the namelist file
         open(unit=10, file=nmlfile, status='old', iostat=fstat)
         if (fstat /= 0) then
            write(*,*) ''
            write(*,*) 'File ',trim(nmlfile),' not found?'
            inited = .true.
            return
         else
            read(10,nml=mcmc,iostat=fstat)
            if (fstat /= 0) then
               write(*,*) 'Could not read mcmc namelist from ',
     b              trim(nmlfile)
               write(*,*) ' status:',fstat
               write(*,*) 'Using defaults'
            end if
            close(unit=10)
         end if

         write(*,*) 'modest uses: sstype = ',sstype,
     b        ' sstrans=',real(sstrans)

c
c         prior mu & sigma from priorsfile
c
         mu = 0.0d0
         sig = 0.0d0
         if (len_trim(priorsfile) > 0) then
            open(333,file=priorsfile,status='old',iostat=fstat)
            if (fstat /= 0) then
               write(*,*) 'could not open file ',trim(priorsfile)
            else
               read(333,*,iostat=fstat) (mu(i), i=1,len)
               if (fstat /= 0) then
                  write(*,*) 'priorsfile should have 2*npar elements'
                  write(*,*) 'npar=',len
                  stop
               end if
               read(333,*,iostat=fstat) (sig(i), i=1,len)
               if (fstat /= 0) then
                  write(*,*) 'priorsfile should have 2*npar elements'
                  write(*,*) 'npar=',len
                  stop
               end if
               close(333)
            end if
         end if

         inited = .true.
      endif

      return
      end
c-------------------------------------------------------------------------
c
c calculate -2*log(p(par), prior "sum-of-squares" for the
c Gaussian case, where the prior parameters are read from a file
c
      function prior2lfun(par,mu,sig,len)
      implicit none
      real*8 prior2lfun
      integer*4 len
      real*8 par(len), mu(len),sig(len)

      real*8 p, e
      integer*4 i, np

      np = count(sig > 0.0d0)

      p = 0.0d0
      if (np > 0) then
         do i=1,len
            if (sig(i)>0.0d0) then
               e = (par(i)-mu(i))/sig(i)
               p = p + e*e
            end if
         end do
      end if
      prior2lfun = p

      return
      end
c-------------------------------------------------------------------------
