#!/bin/bash
# Check if compiler of mpi-wrapper is the same as classic compiler
# CMake standard : 0 for false, 1 for true
# Check intel, gnu, pgi, ibm compilers
# $1 defines kind of compilers
# $2 defines the wrapper
# $3 defines the reference compiler

if [ $1 = 'CC' ]
    then
    if   [ -n "$($2 --version 2>&1 | grep -i "^icc")" ] && [ -n "$($3 --version 2>&1 | grep -i "^icc")" ]; then exit 1
    elif [ -n "$($2 --version 2>&1 | grep -i "^gcc")" ] && [ -n "$($3 --version 2>&1 | grep -i "^gcc")" ]; then exit 1
    elif [ -n "$($2 --version 2>&1 | grep -i "clang")" ] && [ -n "$($3 --version 2>&1 | grep -i "clang")" ]; then exit 1
    elif [ -n "$($2 --version 2>&1 | grep -i "PathScale")" ] && [ -n "$($3 --version 2>&1 | grep -i "PathScale")" ]; then exit 1
    elif [ -n "$($2 -V 2>&1 | grep -i "Cray")" ] && [ -n "$($3 -V 2>&1 | grep -i "Cray")" ]; then exit 1
    elif [ -n "$($2 -qversion 2>&1 | grep "XL C")" ] && [ -n "$($3 -qversion 2>&1 | grep "XL C")" ]; then exit 1
    elif [ -n "$($2 -V 2>&1 | grep "The Portland Group"  )" ] && [ -n "$($3 -V 2>&1 | grep "The Portland Group"  )" ]; then exit 1
    else exit 0
    fi
fi

if [ $1 = 'CXX' ]
    then
    if   [ -n "$($2 --version 2>&1 | grep -i "^icpc")" ] && [ -n "$($3 --version 2>&1 | grep -i "^icpc")" ]; then exit 1
    elif [ -n "$($2 --version 2>&1 | grep -i "^g++")" ] && [ -n "$($3 --version 2>&1 | grep -i "^g++")" ];        then exit 1
    elif [ -n "$($2 --version 2>&1 | grep -i "clang")" ] && [ -n "$($3 --version 2>&1 | grep -i "clang")" ];        then exit 1
    elif [ -n "$($2 --version 2>&1 | grep -i "PathScale")" ] && [ -n "$($3 --version 2>&1 | grep -i "PathScale")" ]; then exit 1
    elif [ -n "$($2 -V 2>&1 | grep -i "Cray")" ] && [ -n "$($3 -V 2>&1 | grep -i "Cray")" ]; then exit 1
    elif [ -n "$($2 -qversion 2>&1 | grep "XL C")" ] && [ -n "$($3 -qversion 2>&1 | grep "XL C")" ];           then exit 1
    elif [ -n "$($2 -V 2>&1 | grep "The Portland Group"  )" ] && [ -n "$($3 -V 2>&1 | grep "The Portland Group"  )" ];  then exit 1
    else exit 0
    fi
fi

if [ $1 = 'FC' ]
    then
    if   [ -n "$($2 --version 2>&1 | grep -i "^ifort")" ] &&  [ -n "$($3 --version 2>&1 | grep -i "^ifort")" ];       then exit 1
    elif [ -n "$($2 --version 2>&1 | grep -i "^GNU Fortran")" ] && [ -n "$($3 --version 2>&1 | grep -i "^GNU Fortran")" ]; then exit 1
    elif [ -n "$($2 -qversion 2>&1 | grep "XL Fortran")" ] && [ -n "$($3 -qversion 2>&1 | grep "XL Fortran")" ];      then exit 1
    elif [ -n "$($2 -V 2>&1 | grep "The Portland Group"  )" ] && [ -n "$($3 -V 2>&1 | grep "The Portland Group"  )" ];   then	exit 1
    elif [ -n "$($2 -V 2>&1 | grep -i "Cray")" ] && [ -n "$($3 -V 2>&1 | grep -i "Cray")" ]; then exit 1
    else exit 0
    fi
fi
