###
### modest mcmc Makefile template
###

# name of the project and exe-file
PROJECT=boxo

# fortran source files (.f or .for)
SRC=$(PROJECT)m.f $(PROJECT)o.f $(PROJECT)i.f
INC=$(PROJECT).inc
NMLIO=nmlio
NML=$(PROJECT).nml

# some linkers don't link MAIN from modest library
#SRC+=modmain.F

## compiler settings for gfortran on Linux
FC=gfortran
FFLAGS=-O3
# LIBS=-L/usr/local/lib -lmodest -llapack -lblas
LIBS=-L/usr/local/lib -lmodest

# compile the project
$(PROJECT): $(SRC) $(INC)
	$(FC) $(FFLAGS) -o $(PROJECT) $(SRC) $(LIBS)

.SUFFIXES: .f .for .o

.f.o:
	$(FC) $(FFLAGS) -c $<
.for.o:
	$(FC) $(FFLAGS) -c $<

# run the project
run: $(PROJECT)$(EXE) $(NML) mcmcinit.nml
	$(NMLIO) $(NML)
	./$(PROJECT)

clean:
	rm -f *~ *.o core

realclean:
	rm -f *~ *.o *.obj *.bak *.tmp core $(PROJECT)$(EXE) modtemp.in \
              *.arg *.est *.sta *.dmp *.eva mod.end mcmcnml.txt

realrealclean:
	rm -f *~ *.o *.obj *.bak *.tmp core $(PROJECT)$(EXE) modtemp.in \
              *.arg *.est *.sta *.dmp *.eva mod.end mcmcnml.txt \
              chain.mat sschain.mat s2chain.mat chain.dat sschain.dat s2chain.dat \
              mcmc*0.dat mcmccovf.dat mcmcparf.dat mcmcsigma2f.dat mcmcmean.dat
