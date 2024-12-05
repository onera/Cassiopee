#!/usr/bin/env python
from distutils.core import setup, Extension
import os

#=============================================================================
# XCore requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Numpy, MPI
# KCore library
#=============================================================================

# Write setup.cfg file
import KCore.Dist as Dist
Dist.writeSetupCfg()

from KCore.config import *

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi(additionalLibPaths,
                                                     additionalIncludePaths)
(mpi4py, mpi4pyIncDir, mpi4pyLibDir) = Dist.checkMpi4py(additionalLibPaths,
                                                        additionalIncludePaths)

prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir]
includeDirs = [numpyIncDir, kcoreIncDir]
libraries = ["kcore"]

ADDITIONALCPPFLAGS = []
if mpi:
    libraryDirs.append(mpiLibDir)
    includeDirs.append(mpiIncDir)
    ADDITIONALCPPFLAGS += ['-D_MPI']
if mpi4py:
    includeDirs.append(mpi4pyIncDir)
if mpi: libraries += mpiLibs

(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Extensions ==================================================================
import srcs
listExtensions = []
listExtensions.append(
    Extension('XCore.xcore',
              sources=['XCore/xcore.cpp']+srcs.cpp_srcs,
              include_dirs=["XCore"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
              extra_link_args=Dist.getLinkArgs()
              ))

# setup ======================================================================
setup(
    name="XCore",
    version="4.0",
    description="Parallel core for *Cassiopee* modules.",
    author="Onera",
    package_dir={"":"."},
    packages=['XCore'],
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
