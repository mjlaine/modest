MCMC code for modest

This code is linked with Modest library together with mcmcrun library
to do MCMC after the estimation task in Modest.

Subroutine mdstmcmc(ar,nr,ai,ni) in mdstmcmc.F90 is called from Modest
to do mcmc.

Some routines in libmcmcrun are replaced for Modest alternatives

mdstssfun.f90       - replaces ssfunction in libmcmcrun
mdstcheckbounds.f90 - replaces checkbounds in libmcmcrun
mdstinitialize.F90  - replaces initialize in libmcmcrun

So -lmdstmcmc must appear before -lmcmcrun in the compiler command line.

Example compiler command line on MAC:

gfortran -o boxo boxom.f boxoo.f boxoi.f /usr/local/share/modest/modmain.F
-lmodest -lmdstmcmc -lmcmcrun -lodepack32 -framework Accelerate

To make one common library instead of four:

* copy libmodest.a libmdstmcmc.a libmcmcrun.a libodepack32.a  to a
  temporary directory 

* remove some files from libmcmcrun.a

>> ar d libmcmcrun.a ssfunction0.o checkbounds0.o initialize.o

* extract mcmcrun and mdstmcmc

>> ar x libmcmcrun.a 
>> ar x libmdstmcmc.a 
>> ar x libodepack.a 

* add the files in modest lib

>> ar ruv libmodest.a *.o

* now libmodest.a has most of the code (lapack are still needed).
