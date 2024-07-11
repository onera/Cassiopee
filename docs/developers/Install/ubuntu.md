# Installation on ubuntu (linux)

## Install dependencies
```shell
sudo apt-get install python3-dev
sudo apt-get install python-numpy
sudo apt-get install scons
sudo apt-get install gcc
sudo apt-get install gfortran
sudo apt-get install hdf5
sudo apt-get install xorg-dev 
```

## Install Cassiopee
```shell
export CASSIOPEE=/d/johndo/Cassiopee
export MACHINE=ubuntu
    
source $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8
cd $CASSIOPEE/Cassiopee
./install
```

## Using apt_get

Find package from keyword:
```shell
apt-cache search <keyword>
```

Install package:
```shell
sudo apt-get install <package>
```