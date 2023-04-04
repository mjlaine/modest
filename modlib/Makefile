# Modest Makefile
# libmodest.a in unix/linux

uname := $(shell uname)
arch  := $(shell uname -m)

# settings
fcompiler=gfortran
include_blas=no
use_mcmc=yes
use_lsode=yes
debug=no

ifeq ($(uname),Linux)
 fcompiler=gfortran
endif
ifeq ($(uname),Darwin)
 fcompiler=gfortran_maci
endif

LIBNAME=modest

.SUFFIXES: .f90 .F90 .F .f .o .lo

# f source
FSRC = clmdiff.f constrai.f csimflex.f cstates.f cstats.f ddrfode.f \
 deuler.f df.f direct.f dirline.f dlsvrr.f dscolr.f eispack.f eulode.f \
 expdes.f exsensit.f falg0.f fimp0.f fimp1.f fode0.f fode1.f hess.f \
 hessup.f inits0.f jac.f jacob.f jacoblm.f linpack.f linpack2.f \
 lmdiff.f \
 machcon.f newfimp.f nrjr.f observat.f \
 odefode.f odeode.f odessa.f optfun0.f optn.f \
 r2.f  readdata.f readfile.f restart.f root0.f salg.f sensit.f \
 simflex.f simp.f simulate.f sstate.f writetas.f writopml.f \
 xdirline.f xerror.f xparser.f
# code that has to be preprocessed
FPRE=modmain.F modest.F modio.F sode.F lsq.F
# f code not used
FEXTRA= lsq1.f pchev.f pchez.f pchqa.f readbini.f ddrode.f
# extra blas routines
FEXTRA2=linpack3.f
FEXTRA3=lodeode.f
FEXTRA4=common3.inc
SRC=$(FSRC) $(FPRE) $(FEXTRA) $(FEXTRA2) $(FEXTRA3) $(FEXTRA4) Makefile

ifeq ($(include_blas),yes)
FSRC+=$(FEXTRA2)
endif
ifeq ($(use_lsode),yes)
FSRC+=$(FEXTRA3)
endif

LIBOBJ := $(patsubst %.f,%.o,$(FSRC))
LIBOBJ += $(patsubst %.F,%.o,$(FPRE))


LIBSOBJ = $(patsubst %.o,%.lo,$(LIBOBJ))

## general f77
FC=f77
ifeq ($(debug),yes)
FFLAGS=-g
else
FFLAGS=-O
endif


ifeq ($(fcompiler),g77) 
FC=g77
ifeq ($(debug),yes)
FFLAGS=-mtune=native -g
else
FFLAGS=-mtune=native -O6 -falign-loops=2 -falign-jumps=2 -falign-functions=2 \
    -fno-f2c -malign-double -fomit-frame-pointer
endif
endif

ifeq ($(fcompiler),gfortran) 
FC=gfortran
ifeq ($(debug),yes)
FFLAGS=-mtune=native -g
else
FFLAGS=-mtune=native -O3 -ftree-vectorize
endif
FFLAGS+=-DGFORTRAN
endif

ifeq ($(fcompiler),gfortran_maci)
FC=gfortran
ifeq ($(debug),yes)
FFLAGS= -g
else
FFLAGS=-O3 -mtune=native -ftree-vectorize -fexternal-blas
endif
FFLAGS+=-DGFORTRAN
endif


ifeq ($(fcompiler),sun)
FC=f90 -f77
ifeq ($(debug),yes)
FFLAGS=-xpentium -g
else
FFLAGS=-xpentium -O5
endif
endif


ifeq ($(fcompiler),intel) 
#FC=ifc
FC=ifort
FFLAGS=-tpp6 -axWMK -fpp -dps -vms -w90 -cm -FI -extend_source
FFLAGS=-tpp7 -ssp -parallel -fpp -dps -vms -w90 -cm -FI
ifeq ($(debug),yes)
FFLAGS+=-g
else
FFLAGS+=-O3 -ip
endif
endif

ifeq ($(fcompiler),absoft) 
FC=f90
# g77 compatibility
FFLAGS=-YEXT_NAMES=LCS -s -B108 -YCFRL=1 -DABSOFT
# do loops once , statics local, checks, 
FFLAGS= -DABSOFT -ej -s -C -N116
FFLAGS= -DABSOFT -s
ifeq ($(debug),yes)
FFLAGS+=-g -m1
FFLAGS_NOOPT:=$(FFLAGS)
else
FFLAGS_NOOPT:=$(FFLAGS)
FFLAGS+=-O -B101
endif
endif

ifeq ($(use_lsode),yes) 
FFLAGS += -DUSE_LSODE
endif
ifeq ($(use_mcmc),yes) 
FFLAGS += -DUSE_MCMC
endif

.F.o:
	$(FC) $(FFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<
.F.lo:
	$(FC) $(FFLAGS) -KPIC -c -o $*.lo $<
.f.lo:
	$(FC) $(FFLAGS) -KPIC -c -o $*.lo $<

all: lib$(LIBNAME).a

lib$(LIBNAME).a: $(LIBOBJ)
	ar ruv lib$(LIBNAME).a $(LIBOBJ)
	ranlib lib$(LIBNAME).a

lib$(LIBNAME).so: lib$(LIBNAME).so.1.0
	ln -sf lib$(LIBNAME).so.1.0 lib$(LIBNAME).so
lib$(LIBNAME).so.1.0: $(LIBSOBJ)
	$(FC) $(FFLAGS) -shared -o lib$(LIBNAME).so.1.0 $(LIBSOBJ) 
#	$(CC) -shared -rdynamic -Wl,-soname,lib$(LIBNAME).so.1 \
#         -o lib$(LIBNAME).so.1.0 $(LIBSOBJ)


install: lib$(LIBNAME).a common3.inc modmain.F
	install -d /usr/local/lib
	install -p -m644 lib$(LIBNAME).a /usr/local/lib
	install -d /usr/local/share/modest
	install -p -m644 common3.inc /usr/local/share/modest
	install -p -m644 modmain.F /usr/local/share/modest

zip: $(SRC)
	zip -j -9 modlib.zip $(SRC)

clean:
	rm -f *.o *.lo *.il *~ core

realclean:
	rm -f *.o *.lo *.il *~ core lib$(LIBNAME).a $(LIBNAME).lib

ifeq ($(fcompiler),absoft) 
odessa.o: odessa.f
	$(FC) $(FFLAGS_NOOPT) -c $<
endif
