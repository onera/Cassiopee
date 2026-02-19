#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#=============================================================================
# Converter requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore library
#=============================================================================

import KCore.Dist as Dist

# Compiler settings must be set in installBase.py / installBaseUser.py
f77compiler = Dist.getf77Compiler()
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibPaths = Dist.getAdditionalLibPaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Test if libhdf5 exists ======================================================
(hdf, hdfIncDir, hdfLibDir, hdflibs) = Dist.checkHdf()

# Test if libpng exists ======================================================
(png, pngIncDir, pngLibDir) = Dist.checkPng()

# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi()
(mpi4py, mpi4pyIncDir, mpi4pyLibDir) = Dist.checkMpi4py()

# Compilation des fortrans ====================================================
if f77compiler is None:
    print("Error: a fortran 77 compiler is required for compiling Converter.")
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT'+str(c)+'='+v+' '
os.system("make -e FC="+f77compiler+" WDIR=Converter/Fortran "+opt)
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir]
includeDirs = [numpyIncDir, kcoreIncDir]
libraries = ["ConverterF", "kcore"]
if hdf:
    libraryDirs.append(hdfLibDir)
    includeDirs.append(hdfIncDir)
if png:
    libraryDirs.append(pngLibDir)
    includeDirs.append(pngIncDir)
ADDITIONALCPPFLAGS = []
if mpi:
    libraryDirs.append(mpiLibDir)
    includeDirs.append(mpiIncDir)
    ADDITIONALCPPFLAGS += ['-D_MPI']
if mpi4py:
    includeDirs.append(mpi4pyIncDir)
if hdf: libraries += hdflibs
if png: libraries.append('png')
if mpi: libraries += mpiLibs
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

ADDITIONALCPPFLAGS = ['-DUSE_C_REGEX'] # for old gcc < 5.0

# Extensions ==================================================================
import srcs
listExtensions = []
listExtensions.append(
    Extension('Converter.converter',
              sources=['Converter/converter.cpp']+srcs.cpp_srcs,
              include_dirs=["Converter"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
              extra_link_args=Dist.getLinkArgs()
              ))
listExtensions.append(
    Extension('Converter.expression',
              sources=['Converter/Expression/Expression.cpp']+srcs.cpp_srcs,
              include_dirs=["Converter"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
              extra_link_args=Dist.getLinkArgs() ) )
# setup ======================================================================
setup(
    name="Converter",
    version="4.1",
    description="Converter for *Cassiopee* modules.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    package_dir={"":"."},
    packages=['Converter'],
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
