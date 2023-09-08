# MODEST code

Modest modelling and estimation package source code.

Folders in this repository:

|          |                             |
|----------|-----------------------------|
| modlib   | modest library src          |
| nmlio    | nmlio executable src        |
| mdstmcmc | modest mcmc library src     |
| mcmcf90  | copy of mcmcrun library src |
| lapack   | copy of reference lapack    |
| odepack  | odepack library from netlib |
| boxo     | modest test model           |


Folder `mcmcf90` is a copy from [https://github.com/mjlaine/mcmcf90]. The `odepack` library is a copy from https://netlib.org/odepack with added `Makefile`, modified directory structure, as well as some auxiliary error handling files.

This version compiles one main modest library that also contains reference implementation of LAPACK and BLAS from  [Reference LAPACK](https://github.com/Reference-LAPACK/lapack). 


## Building the modest library

You need a linux machine or similar and `gfortran` compiler to build the modest library.

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

## Using a docker container

The provided [`Dockerfile`](Dockerfile) can be used to build a docker container that compiles the modest library and can be used to run modest programs.

Building modest container:
```
docker build --rm -t modest .
```

Or directly from github (no nothing needs to be downloaded before):
```
docker build --rm -t modest https://github.com/mjlaine/modest.git
```

Quick test run:
```
docker run --rm -it modest sh -c 'cp -r /opt/modest/boxo . && make -C boxo run'
```

Start an interactive shell and mount current directory to sub directory `work` inside the container:
```
docker run --rm -it -h modest -v $(pwd):/home/modest/work modest
```

Everything in one session:
```
docker build --rm -t modest https://github.com/mjlaine/modest.git
docker run --rm -it -h modest modest
cp -r /opt/modest/boxo .
make -C boxo run
```

Or, if you want to download the `Dockerfile` first:
```
wget http://github.com/mjlaine/modest/blob/master/Dockerfile
docker build --rm -t modest .
...
```

---
marko.laine@fmi.fi
