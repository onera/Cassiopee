# Installation on linux (ubuntu)

Tested on ubuntu 24.04.

## Install dependencies
```shell
sudo apt-get install python3-dev
sudo apt-get install python3-numpy
sudo apt-get install python3-distutils-extra
sudo apt-get install pip
sudo apt-get install scons
sudo apt-get install gcc
sudo apt-get install g++
sudo apt-get install gfortran

sudo apt-get install libopenmpi-dev
sudo apt-get install python3-mpi4py
sudo apt-get install libhdf5-openmpi-dev

sudo apt-get install python3-tk
sudo apt-get install libglu1-mesa-dev
sudo apt-get install mesa-common-dev
sudo apt-get install libosmesa6-dev
sudo apt-get install xorg-dev

sudo apt-get install libocct-foundation-dev
sudo apt-get install libocct-modeling-algorithms-dev
sudo apt-get install libocct-data-exchange-dev
sudo apt-get install libocct-modeling-data-dev
sudo apt-get install libocct-draw-dev

```

## Install Cassiopee
```shell
export CASSIOPEE=/d/johndo/Cassiopee
export MACHINE=ubuntu
    
source $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8
cd $CASSIOPEE/Cassiopee
./install
```

## More apt-get commands

Find package from keyword:
```shell
apt-cache search <keyword>
```

Install package:
```shell
sudo apt-get install <package>
```

Remove package:
```shell
sudo apt-get remove <package>
```