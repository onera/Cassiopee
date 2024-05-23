# Installation sous msys2 (windows)

## Install msys2
Download msys2 (https://www.msys2.org)
Install it

## Install dependencies
In an msys2 mingw64 terminal:

    pacman -S mingw64/mingw-w64-x86_64-gcc
    pacman -S mingw64/mingw-w64-x86_64-python
    pacman -S mingw64/mingw-w64-x86_64-python-numpy
    pacman -S mingw64/mingw-w64-x86_64-scons
    pacman -S mingw64/mingw-w64-x86_64-python-pip
    pacman -S mingw64/mingw-w64-x86_64-python-pip-tools
    pacman -S mingw64/mingw-w64-x86_64-hdf5
    pacman -S mingw64/mingw-w64-x86_64-msmpi
    pacman -S mingw64/mingw-w64-x86_64-oce

## Install Cassiopee

    export CASSIOPEE=/d/johndo/Cassiopee
    export MACHINE=win64
    
    source $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8
    cd $CASSIOPEE/Cassiopee
    ./install

## Some usefull pacman commands

Update system:

    pacman -Syu

Find package matching keyword:

    pacman -Ss <keyword>

Install package:

    pacman -S <package>

List installed packages:

    pacman -Qe
    
Remove package:

    pacman -Rs <package>
