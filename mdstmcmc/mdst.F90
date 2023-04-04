!!!
!!! modest module to store links to modest data
!!!
module mdst

  use mcmcprec
  use mcmcinit, only : initcmatn, condmaxini, sstype, sstrans, updatesigma, S02
  use matutils
  use mcmcmod, only : npar

  implicit none

  public

  integer(4), save :: mdst_nr, mdst_ni
  real(8), pointer, save :: mdst_ar(:)
  integer(4), pointer, save :: mdst_ai(:)

contains

  subroutine mdst_init_arai(ar,nr,ai,ni)
    implicit none
    integer(4) nr, ni
    real(8), target :: ar(:)
    integer(4),target :: ai(:)

    mdst_nr = nr
    mdst_ni = ni
    mdst_ar => ar
    mdst_ai => ai

  end subroutine mdst_init_arai

end module mdst
