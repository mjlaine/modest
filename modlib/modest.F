      subroutine Modest(ar,nr,ai,ni)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)

*     The Main Program -- just technically a subroutine
*
*     PURPOSE
*
*     a program which
*         - estimates unknown parameters in models
*         - analyses identifiability
*         - computes optimal designs
*         - simulates models
*         - computes optima for different targets
*
*     Subroutines called by Main:
*
*         - Readbini:     Reads the input from an unformatted file
*         - Xparser:      Selects the variables to be optimized
*         - Optimizers:   e.g. Simflex, performs the optimization
*         - Wrufout:      Writes the results to unformatted output files
*
*     User written subroutines needed:
*
*         - Falg:         The model subr. if the model algebraic
*         - Fimp:              - " -      if the model implicit algebraic
*         - Fode:              - " -      if the model described by ODE's.
*         - Fpde:              - " -      if the model described by PDE's
*         - Fdae:              - " -      if the model described by DAE system
*         - Penalty:      Describes the non-standard constraints. Dummy, if
*                         not needed.
*

*     PARAMETERS AND COMMON BLOCKS

      include 'common3.inc'

*     Read model & problem definitions and data.
*     The information is passed via common blocks.

      call readbini(ar,nr,ai,ni)   ! read the var from temp file

      call readdata(ar,nr,ai,ni,   ! read data from data files
     &              ar(pxdata),dnx,dmostobs,dnsets,
     &              ar(pydata),dmostyda,
     &              ar(pstates0),dnstatea,
     &              ar(pt0),
     &              ar(ptf),
     &              ai(pnobs),
     &              ai(pnydata),
     &              ar(pweight))

*     Select the operation according to the task


      if    (task(1:3).eq.'est') then
         if (dnest.lt.1) then
            write(*,*) 'Nothing to estimate!'
            stop
         endif
         call optimization(ar,nr,ai,ni)
c mcmc stuff
#ifdef USE_MCMC
         call mdstmcmc(ar,nr,ai,ni)
#endif

      elseif (task(1:3).eq.'exp') then
         if (dnexp.lt.1) then
            write(*,*) 'Nothing to optimize!'
            stop
         endif
         call optimization(ar,nr,ai,ni)
#ifdef USE_MCMC
         call mdstexpanneal(ar,nr,ai,ni)
#endif

      elseif (task(1:3).eq.'opt') then
         if (dnopt.lt.1) then
            write(*,*) 'Nothing to optimize!'
            stop
         endif
         call optimization(ar,nr,ai,ni)

      elseif (task(1:3).eq.'sim') then
         call simulation(ar,nr,ai,ni)

      else if (task(1:3).eq.'sen') then
         call sensitiv(ar,nr,ai,ni)


      else

         write(*,*) 'The task not specified!'

      endif

      open(unit=3,file = 'mod.end',status = 'unknown')
      close(3)

      stop
      end

       subroutine optimization (ar,nr,ai,ni)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)


*     LOCAL PARAMETERS

*    COMMON BLOCKS

      include 'common3.inc'

      if (restartfile(1:5) .ne. 'dummy') then
          call restart(ar,nr,ai,ni,ar(pxaux))
      endif

      if (optimizer(1:5).eq.'simfl') then

         call csimflex(ar,nr,ai,ni)

      elseif (optimizer(1:3).eq.'lev') then

         call clmdiff(ar,nr,ai,ni)

      elseif (optimizer(1:5).eq.'simle') then

       call csimflex(ar,nr,ai,ni)
       call clmdiff(ar,nr,ai,ni)

      else

         write(*,*) 'No such optimizer'
         stop

      endif

*     write the result of the optimization and the state of the
*     system at the optimum in different files

      if(debug .gt. 10) write(*,*) 'calling writetas'

      call writetas(ar,nr,ai,ni)

      return
      end

*end of subroutine optimize


      subroutine simulation(ar,nr,ai,ni)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)


*     PARAMETERS AND COMMON BLOCKS

      include 'common3.inc'

        if (dnsim.gt.0) then
         if (restartfile(1:5) .ne. 'dummy') then
           call restart(ar,nr,ai,ni,ar(pxsim))
         endif
        endif

         call Xparser(ar,nr,ai,ni,
     &                'sel',ar(pxsim),dnsim,ai(plsim))

         if(debug .gt. 10) write(*,*) 'calling simulate'
         call simulate(ar,nr,ai,ni,
     &             ar(pxsim),dnsim,ai(plsim),
     &             ar(pdir),dndir,ai(pnsimste),
     &             ar(pxdir),ar(pcoline),dmostste,dnevapoi,
     &             ai(pnobs),dnsets,
     &             ar(pstates),ar(pstates0),dnstatea,dmostobs,dnsavese,
     &             ar(pxdata),dnx,ar(pyest),dmostyda,ai(pnydata),
     &             ar(pgpar),max(dngpar,1),
     &             ar(plpar),max(dnlpar,1))

      return
      end

*end of subroutine simulation

      subroutine sensitiv(ar,nr,ai,ni)

      implicit none
*     ARGUMENTS
      integer*4 nr,ni
      real*8    ar(nr)
      integer*4 ai(ni)


*     PARAMETERS AND COMMON BLOCKS

      include 'common3.inc'

         if (restartfile(1:5) .ne. 'dummy') then
           call restart(ar,nr,ai,ni,ar(pxsen))
         endif

         call Xparser(ar,nr,ai,ni,
     &                'sel',ar(pxsen),dnsen,ai(plsen))

         if(debug .gt. 10) write(*,*) 'calling sensit'
         call sensit(ar,nr,ai,ni,
     &             ar(pxsen),dnsen,ai(plsen),
     &             ai(pnobs),dnsets,dmostobs,
     &             ar(pxdata),dnx,ai(pnydata),
     &             ar(pdir),dndir,ai(pnsimste),
     &             ar(pxdir),ar(pcoline),dmostste,dnevapoi,
     &             ar(poptline),ar(poptgrid),ai(pigrid),dncontou)

      return
      end

*end of subroutine sensitiv

