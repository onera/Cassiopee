from distutils.core import setup, Extension
#from setuptools import setup, Extension
import os

#=============================================================================
# Converter requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore library
#=============================================================================

# Write setup.cfg file
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

from KCore.config import *

# Test if libhdf5 exists ======================================================
(hdf, hdfIncDir, hdfLibDir, hdflibs) = Dist.checkHdf(additionalLibPaths,
                                                     additionalIncludePaths)

# Test if libnetcdf exists ======================================================
(netcdf, netcdfIncDir, netcdfLibDir, netcdflibs) = Dist.checkNetcdf(additionalLibPaths,
                                                                    additionalIncludePaths)

# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi(additionalLibPaths,
                                                     additionalIncludePaths)
(mpi4py, mpi4pyIncDir, mpi4pyLibDir) = Dist.checkMpi4py(additionalLibPaths,
                                                        additionalIncludePaths)

# Compilation des fortrans ====================================================
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir]
includeDirs = [numpyIncDir, kcoreIncDir]
libraries = ["converter", "kcore"]
if hdf:
    libraryDirs.append(hdfLibDir)
    includeDirs.append(hdfIncDir)
if netcdf:
    libraryDirs.append(netcdfLibDir)
    includeDirs.append(netcdfIncDir)

ADDITIONALCPPFLAGS = []
if mpi:
    libraryDirs.append(mpiLibDir)
    includeDirs.append(mpiIncDir)
    ADDITIONALCPPFLAGS += ['-D_MPI']
if mpi4py:
    includeDirs.append(mpi4pyIncDir)

if hdf:
    for l in hdflibs: libraries.append(l)
if netcdf:
    for l in netcdflibs: libraries.append(l)

if mpi: libraries += mpiLibs

(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Extensions ==================================================================
listExtensions = []
listExtensions.append(
    Extension('Converter.converter',
              sources=['Converter/converter.cpp'],
              include_dirs=["Converter"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
              extra_link_args=Dist.getLinkArgs()
              ) )
import srcs
if srcs.EXPRESSION:
    listExtensions.append(
        Extension('Converter.expression',
                  sources=['Converter/Expression/Expression.cpp'],
                  include_dirs=["Converter"]+additionalIncludePaths+includeDirs,
                  library_dirs=additionalLibPaths+libraryDirs,
                  libraries=libraries+additionalLibs,
                  extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
                  extra_link_args=Dist.getLinkArgs() ) )

# setup ======================================================================
setup(
    name="Converter",
    version="4.0",
    description="Converter for *Cassiopee* modules.",
    author="ONERA",
    url="https://cassiopee.onera.fr",
    packages=['Converter'],
    package_dir={"":"."},
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
