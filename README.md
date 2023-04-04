# MODEST code

Modest modelling and estimation package source code.

Folders in this repository:

|          |                             |
|----------|-----------------------------|
| modlib   | modest library src          |
| nmlio    | nmlio executable src        |
| mdstmcmc | modest mcmc library src     |
| mcmcf90  | copy of mcmcrun library src |
| boxo     | modest test model           |
| odepack  | odepack library from netlib |


Folder `mcmcf90` is a copy from [https://github.com/mjlaine/mcmcf90]. The `odepack` library is a copy from https://netlib.org/odepack with added `Makefile`, modified directory structure, as well as some auxiliary error handling files.


## Building the modest library

You need a linux machine or similar and `gfortran` compiler to build the modest library.

For compiling and running modest executables, you need both `blas` and `lapack` numerical libraries. 

Each folder has a separate `Makefile` for building the code. The [`Makefile`](Makefile) in the main folder builds all libraries. Command `make install` builds everything and copies `libmodest.a` to `/usr/local/lib`.

Folder [`boxo`](boxo) can be used for testing:
```
cd boxo
make
make run
```

Everything in one session:
```
git clone --recurse-submodules http://github.com/mjlaine/modest
cd modest
make install
cd boxo
make run
```

## Using docker container

The provided [`Dockerfile`](Dockerfile) can be used to build a docker container that builds the modest library and can be used to run modest programs.

Building modest container:
```
docker build --rm -t modest .
```

Quick test run:

```
docker run --rm -it modest sh -c 'cp -r /opt/modest/boxo . && make -C boxo run'
```

Start interactive shell:
```
docker run --rm -it -h modest modest bash
```

Everything in one session:
```
wget http://github.com/mjlaine/modest/blob/master/Dockerfile
docker build --rm -t modest .
docker run --rm -it -h modest modest bash
cp -r /opt/modest/boxo .
make -C boxo run
```

---
marko.laine@fmi.fi
