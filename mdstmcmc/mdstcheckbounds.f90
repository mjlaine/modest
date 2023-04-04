!!! checkbounds for modest
function checkbounds(theta)

  use mdst, only :  mdst_nr, mdst_ni, mdst_ar, mdst_ai, npar

  implicit none
  real*8 theta(:)
  logical checkbounds

  include 'mdst_inc.h'

  checkbounds = mdstbnd(mdst_ar,mdst_nr,mdst_ai,mdst_ni,theta,npar)
  return
  
end function checkbounds

