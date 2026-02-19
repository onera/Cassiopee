from setuptools import setup, Extension
import os

#=============================================================================
# XCore requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Numpy, MPI
# KCore library
#=============================================================================
# Compiler settings must be set in installBase.py / installBaseUser.py
import KCore.Dist as Dist
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibPaths = Dist.getAdditionalLibPaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi()
(mpi4py, mpi4pyIncDir, mpi4pyLibDir) = Dist.checkMpi4py()

prod = os.getenv("ELSAPROD") or 'xx'

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

(ok, libs, paths) = Dist.checkCppLibs()
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
    version="4.1",
    description="Parallel core for *Cassiopee* modules.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    package_dir={"":"."},
    packages=['XCore'],
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
