!!! ssfunction for modest

function ssfunction(theta,npar,ny) result(ss)

  use mdst, only :  mdst_nr, mdst_ni, mdst_ar, mdst_ai, sstype, sstrans
  
  implicit none
  integer*4 npar, ny
  real*8 theta(npar)
  real*8 ss(ny)

  include 'mdst_inc.h'

  ss = mdstlsqs(mdst_ar,mdst_nr,mdst_ai,mdst_ni,theta,npar,ny,sstype,sstrans)

end function ssfunction
