        subroutine restart(ar,nr,ai,ni,
     &                     xaux)

*       called by Readfile

*       This subroutine is called if restartfile is given and its name is not
*       'dummy' and  optimization is used. It reads the nest, nexy or nopt
*       parameters from a given file and calls then xparser to put them into
*       their proper places. Thus the values of the parameters to be optimized
*       given in nml or binary input file are overrided. It is mostly used to
*       restart an optimization that did not meet the given tolerances.

        implicit none
*       ARGUMENTS
        integer*4 nr,ni
        real*8    ar(nr)
        integer*4 ai(ni)

*       local variables

        integer*4  i
        real*8     xaux(*)

        include 'common3.inc'

        if(debug .gt. 10) write(*,*) 'in restart'

        open(21,file=restartfile,status='old')

        if (task(1:3) .eq. 'est') then

             read(21,*) (xaux(i),i=1,dnest)
             if(debug .gt. 10) write(*,*) 'calling xparser'
             call xparser(ar,nr,ai,ni,'par',xaux,dnest,ai(plest))

        elseif (task(1:3) .eq. 'exp') then

             read(21,*) (xaux(i),i=1,dnexp)
             if(debug .gt. 10) write(*,*) 'calling xparser'
             call xparser(ar,nr,ai,ni,'par',xaux,dnexp,ai(plexp))

        elseif (task(1:3) .eq. 'opt') then

             read(21,*) (xaux(i),i=1,dnopt)
             if(debug .gt. 10) write(*,*) 'calling xparser'
             call xparser(ar,nr,ai,ni,'par',xaux,dnopt,ai(plopt))

        elseif (task(1:3) .eq. 'sen') then

             read(21,*) (xaux(i),i=1,dnsen)
             if(debug .gt. 10) write(*,*) 'calling xparser'
             call xparser(ar,nr,ai,ni,'par',xaux,dnsen,ai(plsen))

        elseif (task(1:3) .eq. 'sim') then

             read(21,*) (xaux(i),i=1,dnsim)
             if(debug .gt. 10) write(*,*) 'calling xparser'
             call xparser(ar,nr,ai,ni,'par',xaux,dnsim,ai(plsim))


        endif

        close(21)

        return
        end
