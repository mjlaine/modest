!!! $Id: mdstinitialize.F90,v 1.2 2015/01/29 11:31:37 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: mdstinitialize.F90
!!! Purpose: initialization for modest
!!!
!!! Marko Laine

subroutine initialize(par0,npar,cmat0,initcmatn,sigma2,nobs,nycol)

  use mdst, only :  mdst_nr, mdst_ni, mdst_ar, mdst_ai, &
       updatesigma, condmaxini, sstype, sstrans, S02, dbl, & 
       loaddata, writedata


  implicit none
  integer, intent(inout) :: npar, initcmatn, nycol
  integer, intent(inout), allocatable :: nobs(:)
  real(kind=dbl), intent(inout), allocatable :: par0(:), cmat0(:,:)
  real(kind=dbl), intent(inout), allocatable :: sigma2(:)

  real(kind=dbl), pointer :: cmat(:,:)
  integer :: allocstat, fstat

  include 'mdst_inc.h'   ! interfaces to modest ar and ai

  write(*,*) 'note: in mdst initialize'

  call mdstnpar(mdst_ar,mdst_nr,mdst_ai,mdst_ni,npar,nycol) ! fetch npar and nycol

  allocate(sigma2(nycol), stat=allocstat)
  if (allocstat /= 0) stop 'ERROR: allocate sigma2'
  allocate(nobs(nycol), stat=allocstat)
  if (allocstat /= 0) stop 'ERROR: allocate nobs'

  allocate(par0(npar), stat=allocstat)
  if (allocstat /= 0) stop 'allocate'
  
  allocate(cmat0(npar,npar), stat=allocstat)
  if (allocstat /= 0) stop 'allocate'

  allocate(cmat(npar,npar), stat=allocstat)
  if (allocstat /= 0) stop 'allocate'

  call mdstpar(mdst_ar,mdst_nr,mdst_ai,mdst_ni,par0,npar)   ! fetch par, 
  ! the current parameter value
  ! fetch cmat, sigma2, nobs
  call mdstcov(mdst_ar,mdst_nr,mdst_ai,mdst_ni,cmat,npar,sigma2,nobs,nycol, &
       condmaxini,sstype,sstrans)

  !! if getting cmat failed
  if (cmat(1,1) == -1.0_dbl) then
     deallocate(cmat)
     write(*,*) 'No cmat from modest trying mcmccov.dat' 
     call loaddata('mcmccov.dat',cmat,fstat)
     if (fstat /= 0 .or. size(cmat,1) /= npar .or. &
          size(cmat,2) /= npar ) then
        write(*,*) 'ERROR: Error reading file mcmccov.dat, status:',fstat
        write(*,*) '       No proposal covariance. I will stop now'
        stop
     end if
  else
     call writedata('mcmccov0.dat',cmat)
  end if
  call writedata('mcmcpar0.dat',par0)
  call writedata('mcmcsigma20.dat', &
       transpose(reshape( (/sigma2, dble(nobs)/), (/nycol,2/) )))
  
  cmat0 = cmat
  deallocate(cmat)

  if ( updatesigma == 0 .and. S02 > 0.0_dbl ) then
     sigma2 = S02
  end if

end subroutine initialize
