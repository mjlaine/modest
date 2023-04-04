
      subroutine writetas(ar,nr,ai,ni)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

*     LOCAL PARAMETERS

      integer*4  nsim,nest,nexp,nsen,nopt,nesopsi

      include 'common3.inc'

      if (debug .gt. 11) write(6,*) ' in writetas'

      if (outputf(1:4).eq.'matl') then

*        write out for Matlab

         call writopml(ar,nr,ai,ni,
     & ar(pxest),max(dnest,1),ai(plest),
     & ar(pxexp),max(dnexp,1),ai(plexp),
     & ar(pxopt),max(dnopt,1),ai(plopt),
     & ar(pxaux),max(dnesopsi,1),ai(plaux),
     & ai(pnobs),dnsets,
     & ar(pstates),ar(pstates0),dnstatea,dmostobs,dnsavese,
     & ar(pxdata),dnx,
     & ar(pydata),ar(pweight),
     & ar(pyest),ar(pxjac),ar(pjtj),ar(pjtj1),
     & dmostyda,ai(pnydata),
     & ar(pgpar),max(dngpar,1),
     & ar(plpar),max(dnlpar,1))

      elseif (outputf(1:4).eq.'asci') then

*        write out in ascii files

c         call writeasc(ar,nr,ai,ni,

      elseif (outputf(1:4).eq.'exce') then

*         write out for Excel

c         call writeexc(ar,nr,ai,ni,

*        write out for Mathematica

c         call writemat(ar,nr,ai,ni,

      else

        write(6,*) 'no such output format'
        stop

      endif

      return
      end


