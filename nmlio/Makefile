###
### nmlio makefile
###

## $Revision: 1.3 $  $Date: 2007/01/10 09:41:55 $

SRC=smbl.f nmlunix.f nmlpc.f common3.inc Makefile

# gfortran
FC=gfortran
FFLAGS=-O2
RM=rm -f
DEFINES=-DGFORTRAN
EXE=
# intel fortran
#FC=ifort
#FFLAGS=-O2 -fpp -vms -w90 -cm -FI
#RM=rm -f
#DEFINES=
#EXE=
# g77
#FC=g77
#FFLAGS=-O2
#RM=rm -f
#DEFINES=
#EXE=
# Visual fortran
#FC=df
#FFLAGS=-Ox -G5 -fpp -nologo
#RM=del
#DEFINES=-DWIN32
#EXE=.exe

nmlio: smbl.f nmlunix.f
	$(FC) $(FFLAGS) $(DEFINES) -o nmlio$(EXE) smbl.f nmlunix.f 

zip:
	zip -9 nmlio.zip $(SRC)

clean:
	$(RM) *.o *.obj *~ *.bak core

realclean:
	$(RM) nmlio *.o *~ *.bak core modtemp.in

install: nmlio
	install -d /usr/local/bin
	install nmlio /usr/local/bin/
