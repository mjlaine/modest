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
| doc      | old documentation           |


Folder `mcmcf90` is a copy from [https://github.com/mjlaine/mcmcf90]. The `odepack` library is a copy from https://netlib.org/odepack with added `Makefile`, modified directory structure, as well as some auxiliary error handling files. There are some old instructions in the folder [`doc`](doc).

This version creates one main modest library `libmodest.a` that also contains reference implementation of LAPACK and BLAS from  [Reference LAPACK](https://github.com/Reference-LAPACK/lapack). 


## Installing in Linux

You need `gfortran` compiler to build the modest library.

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

Each folder has a separate `Makefile` for building the code. The [`Makefile`](Makefile) in the main folder builds all libraries and makes a combined `libmodest.a` in the main folder.  Command `make install` builds everything and copies `libmodest.a` to `/usr/local/lib` and `nmlio` to `/usr/local/bin`.

## Installing in Mac

Use instruction for Linux. You can use the `gfortran` compiler provided by [Homebrew package manager](https://brew.sh/) and its [gcc formula](https://formulae.brew.sh/formula/gcc).

## Installing in Windows

If you have `git`, `gfortran`, `make` and other build tools already installed in Windows, follow the instructions for Linux. Below are some ideas on installing the needed build tools in Windows. This will probably not work easily, so you would better use Linux (or Linux on Windows, WSL, see below).

#### Install git

In PowerShell use command:
```
winget install --id Git.Git -e --source winget
```

#### Install gfortran and other build tools

Not tested, some options below.

https://www.msys2.org/
After installin msys2, you probably need something similar to:
```
pacman -S mingw-w64-ucrt-x86_64-gcc
pacman -S mingw-w64-ucrt-x86_64-gcc-fortran
pacman -S make
pacman -S git
```

Alternative way to use MSYS2 is given here:
https://cran.r-project.org/bin/windows/Rtools/rtools40.html

Yet another alternative for gfortran and make are:
 - http://www.equation.com/servlet/equation.cmd?fa=fortran
 - http://www.equation.com/servlet/equation.cmd?fa=make
These will not need a separate unix environment, so they can be used inside an usual Windows command prompt.

#### Install Modest

Follow the instructions for Linux.

## Installing in Windows with Linux on Windows (WSL)

 - install [Windows Terminal](https://aka.ms/terminal) using the link or with command `winget install --id Microsoft.WindowsTerminal -e`.
 - install WSL `wsl --install` (needs administrator rights), reboot.
 - Inside Ubuntu WSL: `sudo apt update && sudo apt install gfortran make`
 - Inside Ubuntu WSL install Modest as in Linux.

## Installing using docker container

The provided [`Dockerfile`](Dockerfile) can be used to build a docker container that compiles the modest library and can be used to run modest programs.

Building the modest container:
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

Or, if you want to download and check the `Dockerfile` first:
```
wget http://github.com/mjlaine/modest/blob/master/Dockerfile
docker build --rm -t modest .
...
```


---
marko.laine@fmi.fi
