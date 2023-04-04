#
# odepack Makefile 
#

LIB=odepack
#FC=g77
#FFLAGS= -O

FC=gfortran
# MACI
FFLAGS= -O3 -mtune=native \
 -ftree-vectorize -fomit-frame-pointer -fexternal-blas

# MACI64 if used from Matlab which has 64 bit integers in BLAS
#FFLAGS= -O3 -mtune=native -fdefault-integer-8 \
# -ftree-vectorize -fomit-frame-pointer -fexternal-blas

# for MODEST, make libodepack32.a
#LIB=odepack32
#FFLAGS= -O3 -mtune=native  \
# -ftree-vectorize -fomit-frame-pointer -fexternal-blas


# Linux
#FFLAGS= -O3 -mtune=native \
# -ftree-vectorize -fomit-frame-pointer -fexternal-blas

# ifort
#FC=ifort
#FFLAGS=-FI -w95 -tpp7 -O3 -ssp -parallel


FSRC = $(wildcard src/*.f)
OBJ = $(FSRC:.f=.o)

EXTRA=extra/xerrwd.f extra/xsetf.f extra/iumach.f extra/ixsav.f extra/xsetun.f \
  extra/dgesl.f extra/dgefa.f extra/dgbfa.f extra/dgbsl.f
OBJ += $(EXTRA:.f=.o)

lib$(LIB).a:	$(OBJ)
	ar ru lib$(LIB).a $?
	ranlib lib$(LIB).a

install:	lib$(LIB).a
	install -d /usr/local/lib
	install -p -m644 lib$(LIB).a /usr/local/lib

test: test.o
	g77 test.o -l$(LIB)
	time a.out

zip:
	zip -9r odepack.zip Makefile src extra

clean:
	rm -f *~ src/*.o extra/*.o
 
