!!! $Id: mdstexpanneal.f90,v 1.5 2011/11/04 14:25:21 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mdstmcmc library 
!!! File: mdstexpanneal.f90
!!! Purpose: modest experimental desing annealing, 
!!!          interface from modest to anneal code.
!!!
!!! Marko Laine <Marko.Laine@Helsinki.fi>
!!! ------------------------------------------------------------------------
!!! $Author: mjlaine $
!!! $Revision: 1.5 $
!!! $Date: 2011/11/04 14:25:21 $
!!! ------------------------------------------------------------------------
!!! $Log: mdstexpanneal.f90,v $
!!! Revision 1.5  2011/11/04 14:25:21  mjlaine
!!! precision -> mcmcprec, random -> mcmcrand
!!!
!!! Revision 1.4  2011/03/14 13:05:04  mjlaine
!!! bug fixes, simuind etc
!!!
!!! Revision 1.3  2009/02/10 13:17:42  mjlaine
!!! more cleanups
!!!
!!! Revision 1.2  2009/01/08 15:00:14  mjlaine
!!! mcmc.mod -> mcmcmod.mod, mcmc namelist now public in mcmcmod
!!!
!!! Revision 1.1.1.1  2006/05/21 13:10:04  mjlaine
!!! MCMC f90 code import
!!!
!!! Revision 1.1  2002/05/17 11:14:20  mjlaine
!!! Initial revision
!!!
!!! ------------------------------------------------------------------------
!!!
subroutine mdstexpanneal(ar,nr,ai,ni)

  use matutils
  use mcmcrand
  use matfiles
  use mcmcmod , only : Mcmc_Code_Version
  implicit none
  integer(4)  ::    nr,ni
  real(8)     ::    ar(nr)
  integer(4)  ::    ai(ni)

  real(8) :: ee1, ee2, T, r, rscale
  real(8), allocatable :: xpar(:), xlb(:), xub(:), xparnew(:)
  integer :: nxpar
  character(len=128) :: xchainfile

  real(8), allocatable :: xchain(:,:)
  integer :: nanneal, i, fstat, accept, borderstyle

  namelist /anneal/ T, nanneal, rscale, borderstyle, xchainfile

  include 'mdst_inc.h'   ! interfaces to modest ar and ai
  
  write(*,*) 'starting expdes anneal, code version: ',Mcmc_Code_Version
  write(*,*) '  initializing ...'

  !! anneal.nml control variables initialized here
  nanneal     = 10000
  T           = 0.0d0
  rscale      = 2.0d0
  borderstyle = 1
  xchainfile  = 'anneal.mat'

  accept  = 0 ! count the x values accepted

  open(unit=10, file='anneal.nml', status='old', iostat=fstat)
  if (fstat /= 0) then
     write(*,*) ''
     write(*,*) 'File anneal.nml not found, no annealing!!'
     return
  else
     read(10,nml=anneal,iostat=fstat)
     if (fstat /= 0) then
        write(*,*) 'Error reading file anneal.nml, status:',fstat
        stop
     end if
     close(unit=10)
  end if

  call random_initialize()

  !! fetch and allocate xpar and other modest parameters
  call mdstexpnpar(ar,nr,ai,ni,nxpar)
  allocate(xpar(nxpar), xlb(nxpar), xub(nxpar), xparnew(nxpar))
  allocate(xchain(nanneal,nxpar+1))
  call mdstexppar(ar,nr,ai,ni,xpar,nxpar,xlb,xub)

  write(*,*) 'nxpar = ',nxpar
  write(*,*) 'xpar = ',xpar
  write(*,*) ' xlowerb   xupperb:'
  do i=1,nxpar
     write(*,*) xlb(i), xub(i)
  end do

  ee1 = mdstexptarget(ar,nr,ai,ni,xpar,nxpar)
  if ( T == 0.0d0 ) then
     T = ee1
     write(*,*) 'Using T = ', T
  else
     write(*,*) '       T = ',T
     write(*,*) 'Target 1 = ',ee1
  end if

  xchain(1,1:nxpar) = xpar
  xchain(1,nxpar+1) = ee1

  do i=2,nanneal
     xparnew = xpar + (xub-xlb)/rscale*random_normal(nxpar)
     if (borderstyle == 0) then
        !! reject outbound
        if (any(xparnew < xlb) .or. any(xparnew > xub)) then
           xchain(i,1:nxpar) = xpar
           xchain(i,nxpar+1) = ee1
           cycle
        end if
     else
        !! push outbound to the boundary
        xparnew = max(xlb,xparnew)
        xparnew = min(xub,xparnew)
     end if ! borderstyle
     ee2 = mdstexptarget(ar,nr,ai,ni,xparnew,nxpar)
     call random_number(r)
     if ( (ee2>ee1) .or. (exp(-(ee1-ee2)/T) > r) ) then
        ee1    = ee2
        xpar   = xparnew
        accept = accept+1
     end if
     xchain(i,1:nxpar) = xpar
     xchain(i,nxpar+1) = ee1
  end do

  call writemat(xchainfile,xchain,'xchain',fstat)
  write(*,*) 'Wrote ',trim(xchainfile)

  deallocate(xpar, xlb, xub, xparnew)
  deallocate(xchain)
  call random_eoj()             ! random number generation end-of-job

!  write(*,*) 'accept :',accept, '(',real(accept)/real(nanneal),')'
  write(*,'(A,F6.2,A)') ' Accept: ',100.0*real(accept)/real(nanneal),'%'

  return
end subroutine mdstexpanneal
