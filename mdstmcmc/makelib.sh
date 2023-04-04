#!/bin/sh
# make a combined modest library that contains all the src

mkdir -p /tmp/makelib && cd /tmp/makelib

if [ $? -ne 0 ]; then
    echo 'error creating temporary directory'
    exit 1
fi

cp /usr/local/lib/libmodest.a .
cp /usr/local/lib/libmdstmcmc.a .
cp /usr/local/lib/libmcmcrun.a .
cp /usr/local/lib/libodepack.a .

# remove some files from libmcmcrun.a and libodepack32.a

ar d libmcmcrun.a ssfunction0.o checkbounds0.o initialize.o dump.o
ar d libodepack.a dgesl.o dgefa.o

# extract mcmcrun and mdstmcmc

rm *.o
ar x libmcmcrun.a 
ar x libmdstmcmc.a 
ar x libodepack.a 

# add the files in modest lib

ar -ruv libmodest.a *.o

echo 'libmodest.a is ready in /tmp/makelib'
echo 'copy libmodest.a to /usr/local/lib'

rm -f *.o
