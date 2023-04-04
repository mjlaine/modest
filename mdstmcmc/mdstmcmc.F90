!!! $Id: mdstmcmc.F90,v 1.14 2012/07/03 10:57:34 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mdstmcmc library
!!! File: mdstmcmc.F90
!!! Purpose: interface from modest to mcmc
!!!
!!! Marko Laine <Marko.Laine@Helsinki.fi>
!!! ------------------------------------------------------------------------
!!! $Author: mjlaine $
!!! $Revision: 1.14 $
!!! $Date: 2012/07/03 10:57:34 $
!!! ------------------------------------------------------------------------
!!!
!!! Marko Laine 2001  <marko.laine@helsinki.fi>
!!!
subroutine mdstmcmc(ar,nr,ai,ni)

  use matutils
  use mcmcrand
  use mcmcmod
  use mdst

  implicit none
  integer(4)      nr,ni
  real(8)         ar(nr)
  integer(4)      ai(ni)

  integer(4) fstatus

  write(*,*) 'MCMC code version: ',Mcmc_Code_Version

  call mdst_init_arai(ar,nr,ai,ni)

  call MCMC_init()

  !! Some sstype stuff relevant only to 'modest' code
  if (sstype .eq. 3) then
     sstype = 3
     sstrans = -999.e0_dbl
     sigma2 = 1.0e0_dbl
     updatesigma = 0
     if (verbosity>0) write(*,*) 'note: using Poisson likelihood'
  elseif (sstype .eq. 4) then
     sstype = 3
     sstrans = -999.e0_dbl
     sigma2 = 1.0e0_dbl
     updatesigma = 0
     if (verbosity>0)  write(*,*) 'note: using Binomial likelihood'
     if (verbosity>0)  write(*,*) 'note: DO NOT USE sstype = 4 yet'
  elseif (sstype /= 1 .or. sstrans /= 1.0_dbl) then
     if (verbosity>0) write(*,*) 'note: sstype =',sstype,' sstrans =',real(sstrans)
  endif


  !! some checks
  if (nsimu <= 0) then
     return
  end if
  if (dumpint>0 .and. usrfunlen>0) then
     write(*,*) 'You can not have both dumpint>0 and usrfunlen>0'
     write(*,*) 'Assumig dump'
     usrfunlen = 0
  end if

  call MCMC_run()

  call MCMC_writechains(fstatus)

  call MCMC_cleanup(fstatus)

  write(*,*) 'MCMC done'

  return
end subroutine mdstmcmc
