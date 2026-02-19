# Installation on windows (using msys2)

## Install msys2
Download msys2 (https://www.msys2.org)
and install it.

## Install dependencies
In an msys2 mingw64 terminal:
```shell
pacman -S mingw64/mingw-w64-x86_64-gcc
pacman -S mingw64/mingw-w64-x86_64-gcc-fortran
pacman -S mingw64/mingw-w64-x86_64-python
pacman -S mingw64/mingw-w64-x86_64-python-numpy
pacman -S mingw64/mingw-w64-x86_64-scons
pacman -S mingw64/mingw-w64-x86_64-python-pip
pacman -S mingw64/mingw-w64-x86_64-python-pip-tools
pacman -S mingw64/mingw-w64-x86_64-hdf5
pacman -S mingw64/mingw-w64-x86_64-opencascade
```

For parallel:
First install Microsoft mpi redistribution. Then:
```shell
pacman -S mingw64/mingw-w64-x86_64-msmpi
pacman -S mingw-w64-x86_64-python-mpi4py
```

<!-- For optional tigl:
```shell
pacman -S mingw64/mingw-w64-x86_64-libxml2
pacman -S mingw-w64-x86_64-python-libxstl
pacman -S mingw64/mingw-w64-x86_64-boost
pacman -S mingw64/mingw-w64-x86_64-sympy
``` -->

and export system paths (if not already done):
```shell
export PATH=/mingw64/bin:$PATH
export LD_LIBRARY_PATH=/mingw64/lib:$LD_LIBRARY_PATH
```

## Install Cassiopee
```shell
export CASSIOPEE=<your_path>/Cassiopee
export MACHINE=msys2
    
source $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8
cd $CASSIOPEE/Cassiopee
./install
```

## Some usefull pacman commands

Update system:
```shell
pacman -Syu
```

Find package matching keyword:
```shell
pacman -Ss <keyword>
```

Install package:
```shell
pacman -S <package>
```

List installed packages:
```shell
pacman -Q
```

Remove package:
```shell
pacman -Rs <package>
```