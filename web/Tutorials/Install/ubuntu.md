# Installation sous ubuntu (linux)

## Install dependencies

    sudo apt-get install python-dev
    sudo apt-get install python-numpy
    sudo apt-get install scons
    sudo apt-get install gcc
    sudo apt-get install gfortran
    sudo apt-get install xorg-dev 

## Install Cassiopee

    export CASSIOPEE=/d/johndo/Cassiopee
    export MACHINE=ubuntu
    
    source $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8
    cd $CASSIOPEE/Cassiopee
    ./install