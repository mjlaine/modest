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
| doc      | Old documentation           |


Folder `mcmcf90` is a copy from [https://github.com/mjlaine/mcmcf90]. The `odepack` library is a copy from https://netlib.org/odepack with added `Makefile`, modified directory structure, as well as some auxiliary error handling files. There are some old instructions in folder [`doc`](doc).

This version compiles one main modest library that also contains reference implementation of LAPACK and BLAS from  [Reference LAPACK](https://github.com/Reference-LAPACK/lapack). 


## Building the modest library in Linux

You need a linux machine and `gfortran` compiler to build the modest library.

Installing modest library (`libmodest.a`) and `nmlio` executable:
```
git clone --recurse-submodules http://github.com/mjlaine/modest
cd modest
make
sudo make install
```

Testing with the "boxo" example in folder [`boxo`](boxo):
```
cd boxo
make run
```

#### Some details

Each folder has a separate `Makefile` for building the code. The [`Makefile`](Makefile) in the main folder builds all libraries. Command `make install` builds everything and copies `libmodest.a` to `/usr/local/lib` and `nmlio` to `/usr/local/bin`.


## Installing in Windows

If have `git`, `gfortran`, `make` and other build tools already installed in Windows, follow the instructions for Linux. Below are some ideas on installing the needed build tools in Windows, if they are missing.

#### Install git

In PowerShell use command:
```
winget install --id Git.Git -e --source winget
```

#### Install gfortran

Not testesd, some options below.

https://www.msys2.org/

https://cran.r-project.org/bin/windows/Rtools/rtools40.html

#### Install Modest

Follow the instructions for Linux.

## Installing in Windows with Linux on Windows (WSL)

 - install Windows Terminal https://aka.ms/terminal or `winget install --id Microsoft.WindowsTerminal -e`.
 - install WSL `wsl --install`, reboot.
 - Inside WSL: `sudo apt update && sudo apt install gfortran make`
 - Inside WSL install Modest as in Linux.

## Installing using docker container

The provided [`Dockerfile`](Dockerfile) can be used to build a docker container that compiles the modest library and can be used to run modest programs.

Building modest container:
```
docker build --rm -t modest .
```

Or directly from github (nothing needs to be downloaded before):
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
