*     The main program which allocates the global arrays
*     and defines other parameters

      program modmain

      integer*4 rsize, isize, machflag

#ifdef ABSOFT
      parameter (rsize = 100000)
      parameter (isize =   5000)
#else
      parameter (rsize = 3500000)
      parameter (isize =  50000)
#endif

      common /machine/ machflag

      integer*4 ai(isize)                ! all allocatable integer tables
      real*8    ar(rsize)                ! all allocatable real tables

       machflag = 1                      ! pc
*      machflag = 2                      ! decstation
*      machflag = 3                      ! convex
*      machflag = 4                      ! HP

      call modest(ar,rsize,ai,isize)

      stop
      end







