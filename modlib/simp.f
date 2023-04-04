      subroutine SIMP(ar,nr,ai,ni,
     &                s, s0, states, nstatea,
     &                yest,mostydat,mostobs,nydata,
     &                nsaveset, xdata, nx, nsets, nobs,
     &                gpar, ngpar, lpar, nlpar)

*     CALLED BY

*     - cstates

*     SUBROUTINES NEEDED

*     - Implicit solvers

*     PURPOSE

*     This subroutine calculates the states in the user given points (in x)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

      integer*4 nstatea,mostydat,mostobs,nsaveset,nx,nsets
      integer*4 ngpar,nlpar
      real*8    s(nstatea)                       ! s is actually the whole ar
      real*8    s0(nstatea,nsets)
      real*8    states(nstatea,mostobs,nsaveset)
      real*8    yest(mostydat,mostobs,nsaveset)
      integer*4 nydata(mostobs,nsets)
      real*8    xdata(nx,mostobs,nsets)
      integer*4 nobs(nsets)
      real*8    gpar(ngpar)
      real*8    lpar(nlpar,nsets)

*     LOCAL VARIABLES

      integer*4 j,k,k1,info
      integer*4 nstates,ncount
      external  newfimp

*     COMMON BLOCKS

      include 'common3.inc'

      if(debug .gt. 11) write(*, *) ' in simp'

      nstates = nstatea - dnsaux                 ! dnsaux comes from common3.inc

*     Fimp gives the model as an implicit algebraic system of equations.
*     This subroutine calculates the states in the user given points (x)
*     by solving the system for 'y' in the explicit form 'y = F(x)'.

      ncount = ncount+1
      k  = iset
      k1 = mod(k,nsaveset) + 1

      call dcopy(nstates,s0(1,k),1,s(1),1)  !subst the initial guess

      if (iobs .le. 0) then     ! all observations in the data set

          do j = 1,nobs(k)

             iobs = j                                !index of observation

             if (equsolver(1:3) .eq. 'new') then
                 
                if(debug .gt. 10) write(*, *) ' calling newton'

                call nrjr(newfimp,nstates,s,ar(psmin),ar(psmax),
     &                    tolnr,resnr,itmaxnr,ar(pwanr),dwanr,
     &                    iprnr,itcntnr,ai(piparnr),ai(piwanr),
     &                    diwanr,dstep,djacnr)
c                call ebroyd(newfimp,nstates,s,ar(psmin),ar(psmax),
c     &                    tolnr,resnr,itmaxnr,ar(pwanr),dwanr,
c     &                    iprnr,itcntnr,info) !ai(piparnr),ai(piwanr),
c     &                    diwanr,dstep,djacnr)
c                call ebroyd(newfimp,nstates,s,ar(psmin),ar(psmax),
c     &                    tolnr,resnr,itmaxnr,ar(pwanr),6,
c     &                    iprnr,itcntnr) !ai(piparnr),ai(piwanr),
c     &                    diwanr,dstep,djacnr)

         
*         Auxiliary states:
              if (dnsaux.gt.0) then
                 call newfimp(nstates, s, ar(pwanr), info)
              endif
              call dcopy(nstatea,s(1),1,states(1,j,k1),1)
c              if (optmonit.eq.4) write(*, *) 'Res in Simp',resnr
            else
              write(*, *) 'the solver not implemented'
              stop
            endif

          enddo  ! j
          iobs = -1

      else       ! the  observation 'iobs' only

          if (equsolver(1:3) .eq. 'new') then

             if(debug .gt. 10) write(*, *) ' calling newton'

             call nrjr(newfimp,nstates,s,ar(psmin),ar(psmax),
     &                 tolnr,resnr,itmaxnr,ar(pwanr),dwanr,
     &                 iprnr,itcntnr,ai(piparnr),ai(piwanr),
     &                 diwanr,dstep,djacnr)
c                call ebroyd(newfimp,nstates,s,ar(psmin),ar(psmax),
c     &                    tolnr,resnr,itmaxnr,ar(pwanr),dwanr,
c     &                    iprnr,itcntnr,info) !ai(piparnr),ai(piwanr),
c     &                    diwanr,dstep,djacnr)
c                call ebroyd(newfimp,nstates,s,ar(psmin),ar(psmax),
c     &                    tolnr,resnr,itmaxnr,ar(pwanr),6,
c     &                    iprnr,itcntnr) !ai(piparnr),ai(piwanr),
c     &                    diwanr,dstep,djacnr)

*         Auxiliary states:
              if (dnsaux.gt.0) then
                 call newfimp(nstates, s, ar(pwanr), info)
              endif
              call dcopy(nstatea,s(1),1,states(1,iobs,k),1)
c              if (optmonit.eq.4) write(*, *) 'Res in Simp',resnr
          else
            write(*, *) 'the solver not implemented'
            stop
          endif
      
      endif ! iobs
      return
      end

