## Main Makefile for Modest/MCMC
## Makes modest, odepack, and mcmc libraries and nmlio executable

SUBDIRS=modlib nmlio odepack mcmcf90 mdstmcmc lapack
ISUBDIRS=modlib nmlio

all: $(SUBDIRS)

.PHONY: $(SUBDIRS)

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
	cp lapack/make.inc.example lapack/make.inc
	$(MAKE) -C lapack/SRC double
	$(MAKE) -C lapack/BLAS

mcmcf90:
	if [ ! -d "mcmcf90" ]; then unzip mcmcf90.zip; fi	
	$(MAKE) -C mcmcf90

mdstmcmc:
	$(MAKE) -C mdstmcmc

modest: all
	mkdir -p combine
	cp modlib/libmodest.a combine/
	cp mdstmcmc/libmdstmcmc.a combine/
	cp mcmcf90/libmcmcrun.a combine/
	cp odepack/libodepack.a combine/
	cp lapack/liblapack.a combine/
	cp lapack/librefblas.a combine/
	ar d combine/libmcmcrun.a ssfunction0.o checkbounds0.o initialize.o dump.o
	ar d combine/libodepack.a dgesl.o dgefa.o
	(cd combine \
          && ar x libmcmcrun.a \
          && ar x libmdstmcmc.a \
	  && ar x libodepack.a \
	  && ar x liblapack.a \
	  && ar x librefblas.a \
	  && ar -ruv libmodest.a *.o)
	cp combine/libmodest.a .
	rm -rf combine


install: install2 modest
	install -p -m644 libmodest.a /usr/local/lib

install2: all
	@for i in $(ISUBDIRS); do (cd $$i; $(MAKE) install); done

clean:
	@for i in $(SUBDIRS); do (cd $$i; $(MAKE) clean); done

zip:
	zip -r modest.zip $(SUBDIRS) -x *~ *.o *.a */CVS/* *.zip

test:
	zip -r /tmp/t.zip nmlio -x *~ *.o *.a */CVS/* *.zip
