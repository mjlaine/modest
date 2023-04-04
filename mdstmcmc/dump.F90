!!! $Id: dump.F90,v 1.2 2012/07/02 09:22:55 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mdstmcmc library
!!! File: dump.F90
!!! Purpose: handle the Modest dump files
!!!
!!! Marko Laine 2008 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
!!! contains: dump_init, dump_end, dump
!!!


!!! now dummy until I fix the code below
! subroutine dump_init()
!   implicit none
! end subroutine dump_init

! subroutine dump_end()
!   implicit none
! end subroutine dump_end

! subroutine dump(oldpar)
!   use mcmcprec
!   implicit none
!   real(kind=dbl), intent(in) :: oldpar(:)
! end subroutine dump

!! maybe now it works


subroutine dump_init()
  !! dump states at every dumpint interval
  !! channel 111 is written in mdstdump
  !! channel 112 is s2chain.dmp
  !! channel 113 is chain.eva

  use mdst
  use mcmcmod
  use mcmcinit

  implicit none
  integer :: fstat 
  !! call dump (this should be done in lsqfun)
  if (dumpint > 0) then
     ! no need to do it in 1st run, allready done in modest
     !       call mdstdump(ar,nr,ai,ni,oldpar,npar)
     ! open dump file for s2chain
     if (updatesigma /= 0) then
        open(unit=112,file=s2dmpfile,status='replace',iostat=fstat)
        if (fstat /= 0) then
           write(*,*) 'ERROR: File ',trim(s2dmpfile),' could not be opened?'
           write(*,*) ' I will stop now, sorry'
           stop
        else
           write(*,*) 'Opened ',trim(s2dmpfile),' for s2 dump'
           call writevec(112,s2chain(1,1:nycol))
        end if
     end if
     open(unit=113,file=chdmpfile,status='replace',iostat=fstat)
     if (fstat /= 0) then
        write(*,*) 'ERROR: File ',trim(chdmpfile),' could not be opened?'
        write(*,*) ' I will stop now, sorry :-('
        stop
     else
        write(*,*) 'Opened ',trim(chdmpfile),' for chain evafile'
        call writevec(113,chain(1,1:npar))
     end if
  end if
end subroutine dump_init

!!! close the dump files
subroutine dump_end()

  implicit none
  logical :: fo

!!! unit 111 is the dump file opened in mdstdump.f
  inquire(unit=111,opened=fo)
  if (fo) then
     close(111)
  end if
!!! unit 112 is the dump file for s2chain.dmp ascii dump
  inquire(unit=112,opened=fo)
  if (fo) then
     close(112)
  end if
!!! unit 113 is the dump file for chain.eva ascii dump
  inquire(unit=113,opened=fo)
  if (fo) then
     close(113)
  end if

end subroutine dump_end

!!! do dump of states etc with Modest
subroutine dump(oldpar)

  use mdst
  use mcmcmod
  use mcmcinit

  implicit none
  real(kind=dbl), intent(in) :: oldpar(:)

  include 'mdst_inc.h'

  if (dumpint > 0 .and. (mod(simuind,dumpint) == 0) ) then
     call mdstdump(mdst_ar,mdst_nr,mdst_ai,mdst_ni,oldpar,npar)
     if (updatesigma /= 0) then
        call writevec(112,sigma2)
     end if
     call writevec(113,chain(chainind,1:npar))
  end if

end subroutine dump

