# $Id: Makefile,v 1.5 2015/01/29 11:31:37 mjlaine Exp $
## mcmc f90  Makefile for modest library

LIBNAME=mdstmcmc

MCMCF90=../mcmcf90

# Try to guess the system
sys := $(shell uname)$(shell uname -m)

# sys := aix

# default is Win32 with gfortran
system=win32
system=maci
f90compiler=gfortran

ifeq ($(sys),Linuxi686)
 system=linux
 f90compiler=gfortran
endif
ifeq ($(sys),Linuxx86_64)
 system=linux64
 f90compiler=gfortran
endif
ifeq ($(sys),Linuxia64)
 system=linux
 f90compiler=intel
endif
ifeq ($(sys),Darwini386)
 system=maci
 f90compiler=gfortran
endif
ifeq ($(sys),Darwinx86_64)
 system=maci
 f90compiler=gfortran
endif
ifeq ($(sys),aix)
 system=aix
 f90compiler=xlf90
endif

# produce debug code?
debug=no
# debug=yes

### no user modification below this, but check linux.mk or win32.mk ###

ifeq ($(system),win32)
  include ${MCMCF90}/win32.mk
endif
ifeq ($(system),mac)
  include ${MCMCF90}/mac.mk
endif
ifeq ($(system),maci)
  include ${MCMCF90}/maci.mk
endif
ifeq ($(system),linux)
  include ${MCMCF90}/linux.mk
endif
ifeq ($(system),linux64)
  include ${MCMCF90}/linux64.mk
endif
ifeq ($(system),aix)
  include ${MCMCF90}/aix.mk
endif

LIBSRC =  mdst.F90 mdstcov.F mdstlsq.f mdstlsqs.f mdstmcmc.F90 mdstnpar.f \
 mdstpar.f mdstbnd.f mdstufun.f usrfun_dump.F90  dump.F90 \
 mdstexpanneal.f90 mdstexptarget.f mdstexpnpar.f mdstexppar.f \
 mdstdump.f dumpyest.F90 mdstlsq2.f \
 mdstssfun.F90 mdstcheckbounds.f90 mdstinitialize.F90

LIBOBJ := $(patsubst %.f,%.$(O),$(filter %.f,$(LIBSRC))) \
          $(patsubst %.F,%.$(O),$(filter %.F,$(LIBSRC))) \
          $(patsubst %.f90,%.$(O),$(filter %.f90,$(LIBSRC))) \
          $(patsubst %.F90,%.$(O),$(filter %.F90,$(LIBSRC))) \
          $(patsubst %.c,%.$(O),$(filter %.c,$(LIBSRC)))

MFILES = 

SRC = $(LIBSRC) mdst_inc.h common3.inc mcmcinit.nml anneal.nml \
      Makefile Makefile.modest $(MFILES) Readme_modest.txt makelib.sh

LIB = $(HOME)/WinHome/lib/$(LIBNAME).lib

ZIPSRC = $(SRC) $(TARGET) $(LIB)

## how to make these
FMODULES=mcmc util matfiles mcmcprec mcmcrand

.SUFFIXES: .a .$(O) .l$(O) .mod .f90 .F90 .F .f

F90FLAGS += -I${MCMCF90}

all: mcmclib

zip: mdstmcmc.zip

mdstmcmc.zip: $(SRC)
	$(ZIP) -9 mdstmcmc.zip $(SRC)
clean:
	$(RM) *.$(O) *.l$(O) core *~ *.mod
realclean: clean
	$(RM) *.o *.obj *.lo *.il *.lobj *.a *.lib *.mod

mcmclib: libmdstmcmc.a

install: mcmclib
	install -d /usr/local/lib
	install -p -m644 libmdstmcmc.a /usr/local/lib

win32: $(LIBNAME).lib

$(LIBNAME).lib: $(LIBOBJ)
	$(AR) $(LIBNAME).lib $(LIBOBJ)

linux: lib$(LIBNAME).a

lib$(LIBNAME).a: $(LIBOBJ)
	$(AR) ruv lib$(LIBNAME).a $(LIBOBJ)
	$(RANLIB) lib$(LIBNAME).a

lib$(LIBNAME).so: lib$(LIBNAME).so.$(MCMC_VERSION)
	ln -sf lib$(LIBNAME).so.$(MCMC_VERSION) lib$(LIBNAME).so
lib$(LIBNAME).so.$(MCMC_VERSION): $(LIBSOBJ)
	$(CC) -shared -Wl,-soname,lib$(LIBNAME).so.$(MCMC_MAJOR_VERSION) \
         -o lib$(LIBNAME).so.$(MCMC_VERSION) $(LIBSOBJ)


.f90.$(O):
	$(F90) $(F90FLAGS) $(defines) -c $<
.F90.$(O):
	$(F90) $(F90FLAGS) $(defines) -c $<
.f.$(O):
	$(F77) $(F77FLAGS) $(defines) -c $<
.F.$(O):
	$(F77) $(F77FLAGS) $(defines) -c $<
.c.$(O):
	$(CC) $(CCFLAGS) $(defines) -c $<
# Unix ELF shared library objects
.f90.l$(O):
	$(F90) $(F90FLAGS) -KPIC $(defines) -c -o $*.l$(O) $<
.F90.l$(O):
	$(F90) $(F90FLAGS) -KPIC $(defines) -c -o $*.l$(O) $<
.f.l$(O):
	$(F77) $(F77FLAGS) -KPIC $(defines) -c -o $*.l$(O) $<
.F.l$(O):
	$(F77) $(F77FLAGS) -KPIC $(defines) -c -o $*.l$(O) $<

# Cancel .mod.o rule
%.o : %.mod

# dependencies
mdstmcmc.$(O): mdst.$(O)
mdstssfun.$(O): mdst.$(O)
mdstcheckbounds.$(O): mdst.$(O)
mdstinitialize.$(O): mdst.$(O)
