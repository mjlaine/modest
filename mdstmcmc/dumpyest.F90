!!!
!!! dumpyest -- dump yest to mat file
!!!
subroutine dumpyest(yest,mostydat,nobs, iset, nobsall)

  implicit none
  integer*4 nobs,mostydat,nobsall,iset
  real*8    yest(mostydat,nobs)

  integer, save :: first_time=1
  integer rlen
  integer, save :: nextrec

  integer(kind=4),save :: type  = 0
  integer(kind=4),save :: mrows = 0
  integer(kind=4),save :: ncols = 0
  integer(kind=4),save :: imagf = 0
  integer(kind=4),save :: namelen ! length of name + 1
  character(len=8) :: name = "yestdump"

  if (first_time == 1) then
     ! get record length of integer*4, which sholud be 1 (?)
     inquire(iolength=rlen) type
     ! open dump file, it is closed in mcmccleanup
#ifdef __GFORTRAN__
     open(unit=111,file='yestdump.mat',form='unformatted',status='replace', &
          access='direct',recl=rlen)
#else
     open(unit=111,file='yestdump.mat',form='binary',status='replace', &
          access='direct',recl=rlen)
#endif

     mrows = nobsall*mostydat

     write(111, rec=0*4+1) type
     write(111, rec=1*4+1) mrows
     write(111, rec=2*4+1) ncols
     write(111, rec=3*4+1) imagf
     namelen = len_trim(name) + 1
     write(111, rec=4*4+1) namelen
     write(111, rec=5*4+1) trim(name)//achar(0)
     inquire(111,nextrec=nextrec)
     first_time = 0
  end if

  if (iset == 1) then
     ncols = ncols+1
     write(111,rec=9) ncols
  end if

  write(111,rec=nextrec) transpose(yest(:,:))
  inquire(111,nextrec=nextrec) ! save location

  return
end subroutine dumpyest
