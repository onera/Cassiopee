#!/usr/bin/env bash

source $RECIPE_DIR/build-tools.sh

CC=gcc
CXX=g++
FC=gfortran
#export CXXFLAGS="-Wno-deprecated -fp-model=precise -sox -O3 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2"
#export FFLAGS="-Wno-deprecated -fp-model=precise -sox -O3 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2"
#export CXXFLAGS=""
#export FFLAGS=""

cd cassiopee_$version

export ELSAPROD=conda_gnu_flowsim

# used to pass intel compiler license to the modules build system
#export LM_LICENSE_FILE=$INTEL_LICENSE_FILE


# for finding dependencies in standard paths
export LIBPATH=$PREFIX/lib:$LD_LIBRARY_PATH
export CPPPATH=$PREFIX/include
export LD_LIBRARY_PATH=$LIBPATH

# build the modules, order matters
all_modules=(
KCore \
Converter \
Geom \
Connector \
Compressor \
Generator \
Initiator \
Distributor2 \
Intersector \
Dist2Walls \
RigidMotion \
Transform \
Post \
#CPlot \
)
if [[ $version=="2.3" ]]
then
    modules=${all_modules[*]/Intersector/}
fi

export OMP_NUM_THREADS=24

for module in $modules; do
  pushd $module
  echo "PWD ---> " $PWD
  # use out build settings
  cp -f $RECIPE_DIR/installBase.py .
  rm -rf build
  echo "PREFIX---> " $PREFIX
  ./install $PREFIX
  find ./ -name tmp-libkcore.a
  popd
done

# some installed shell scripts have hardcoded non standard locations
# they just wrap python scripts calls, so we just make the python scripts callable
# and put them in the standard location for the executables
scripts=(
validCassiopee \
cplot \
CCC \
tkCassiopee \
cassiopee
)

for wrapper in ${scripts[*]}; do
   if [ -e $PREFIX/$wrapper.py ] || [ -e $PREFIX/$wrapper ]
   then
       mv $PREFIX/$wrapper* $PREFIX/bin/
       echo "------------------>>>>>>>> "$PREFIX/$wrapper*
   fi
done

