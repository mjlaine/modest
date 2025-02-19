## Main Makefile for Modest/MCMC
## Makes modest, odepack, and mcmc libraries and nmlio executable

SUBDIRS=modlib nmlio odepack mcmcf90 mdstmcmc lapack
ISUBDIRS=modlib nmlio


# Try to make this work on Windows, too
ifeq ($(OS),Windows_NT)
   ifneq ($(strip $(filter %sh,$(basename $(realpath $(SHELL))))),)
      POSIXSHELL = 1
   else
      POSIXSHELL =
   endif
else
   POSIXSHELL = 1
endif
ifneq ($(POSIXSHELL),)
    COPY = cp -f
    MKDIR = mkdir -p
    RM = rm
    RMDIR = rm -rf
    CMDSEP = &&
    PSEP = /
else
    COPY = copy /y
    MKDIR = mkdir
    RM = del
    RMDIR = rmdir /s /q
    CMDSEP = &
    PSEP = \\
endif

all: modest

.PHONY: $(SUBDIRS) modest

unzip:
	@for i in $(SUBDIRS); do \
	  (if [ ! -d "$$i" ]; then unzip $$i.zip; fi) \
	done

modlib:
	$(MAKE) -C modlib

nmlio:
	$(MAKE) -C nmlio

odepack:
	$(MAKE) -C odepack

lapack:
	$(COPY) lapack$(PSEP)make.inc.example lapack$(PSEP)make.inc
	$(MAKE) -C lapack/SRC la_xisnan.o
	$(MAKE) -C lapack/SRC double
	$(MAKE) -C lapack/BLAS

mcmcf90:
	$(MAKE) -C mcmcf90

mdstmcmc:
	$(MAKE) -C mdstmcmc


modest: libmodest.a

OBJ = $(notdir $(wildcard combine/*.o))

libmodest.a: extractobj
	cd combine \
          $(CMDSEP) ar -ruv libmodest.a $(OBJ)
	$(COPY) combine$(PSEP)libmodest.a .
	$(RMDIR) combine

extractobj: $(SUBDIRS)
	$(MKDIR) combine
	$(COPY) modlib$(PSEP)libmodest.a combine$(PSEP)
	$(COPY) mdstmcmc$(PSEP)libmdstmcmc.a combine$(PSEP)
	$(COPY) mcmcf90$(PSEP)libmcmcrun.a combine$(PSEP)
	$(COPY) odepack$(PSEP)libodepack.a combine$(PSEP)
	$(COPY) lapack$(PSEP)liblapack.a combine$(PSEP)
	$(COPY) lapack$(PSEP)librefblas.a combine$(PSEP)
	ar d combine$(PSEP)libmcmcrun.a ssfunction0.o checkbounds0.o initialize.o dump.o
	ar d combine$(PSEP)libodepack.a dgesl.o dgefa.o
	cd combine \
          $(CMDSEP) ar x libmcmcrun.a \
          $(CMDSEP) ar x libmdstmcmc.a \
	  $(CMDSEP) ar x libodepack.a \
	  $(CMDSEP) ar x liblapack.a \
	  $(CMDSEP) ar x librefblas.a

install: modest install2
	install -p -m644 libmodest.a /usr/local/lib

install2: $(SUBDIRS)
	@for i in $(ISUBDIRS); do (cd $$i; $(MAKE) install); done

clean:
	@for i in $(SUBDIRS); do (cd $$i; $(MAKE) clean); done

zip:
	zip -r modest.zip $(SUBDIRS) -x *~ *.o *.a */CVS/* *.zip

test:
	zip -r /tmp/t.zip nmlio -x *~ *.o *.a */CVS/* *.zip
