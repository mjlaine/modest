!!!
!!! dummy version of usrfun
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

  usrfun = 0.0d0
  write(*,*) 'dummy usrfun'
  stop

  return
end function usrfun
