!!!
!!! dump version of usrfun
!!!
!!!
function usrfun(len,s,ns,yest,ny, &
     xdata,nx,nobs, &
     nsaux,nstatea, &
     ydata,weight,mostydat, &
     states0, &
     gpar,ngpar, &
     lpar,nlpar, &
     iset)

  implicit none
  integer*4 len
  real*8 usrfun(len)
  integer*4 ns,nx,nobs,nsaux,nstatea,mostydat,ngpar,nlpar,iset
  integer*4 ny(nobs)
  real*8    s(nstatea,nobs)
  real*8    yest(mostydat,nobs)
  real*8    xdata(nx,nobs)
  real*8    ydata(mostydat,nobs)
  real*8    weight(mostydat,nobs)
  real*8    states0(nstatea)
  real*8    gpar(ngpar)
  real*8    lpar(nlpar)

  integer, save :: first_time=1
  integer i, j, rlen
  integer, save :: nextrec

  integer(kind=4),save :: type  = 0
  integer(kind=4),save :: mrows = 0
  integer(kind=4),save :: ncols = 0
  integer(kind=4),save :: imagf = 0
  integer(kind=4),save :: namelen = 10 ! length of name + 1
  character(len=9) :: name = 'statedump'

  if (first_time == 1) then
     inquire(iolength=rlen) type

#ifdef __GFORTRAN__
     open(unit=111,file='dumppi.mat',form='unformatted',status='replace', &
          access='direct',recl=rlen)
#else
     open(unit=111,file='dumppi.mat',form='binary',status='replace', &
          access='direct',recl=rlen)
#endif

!    ncols = nobs*nstatea
     mrows = nobs*nstatea

     write(111, rec=0*4+1) type
     write(111, rec=1*4+1) mrows
     write(111, rec=2*4+1) ncols
     write(111, rec=3*4+1) imagf
     write(111, rec=4*4+1) namelen
     write(111, rec=5*4+1) trim(name)//achar(0)
     inquire(111,nextrec=nextrec)
     first_time = 0
  end if

! mrows = mrows+1;
! write(111,rec=5) mrows
  ncols = ncols+1
  write(111,rec=9) ncols

  write(111,rec=nextrec) transpose(s(:,:))
  inquire(111,nextrec=nextrec)

  usrfun = 0.0d0

  return
end function usrfun
