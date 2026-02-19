#=============================================================================
# Intersector requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore library
#=============================================================================
import os
from setuptools import setup, Extension
import KCore.Dist as Dist

additionalLibPaths = Dist.getAdditionalLibPaths()
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Test if xcore exists =======================================================
#(xcoreVersion, xcoreIncDir, xcoreLibDir) = Dist.checkModuleCassiopee("XCore")()

# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi()
(mpi4py, mpi4pyIncDir, mpi4pyLibDir) = Dist.checkMpi4py()

# Compilation des fortrans ====================================================
prod = os.getenv("ELSAPROD") or "xx"

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir]
includeDirs = [numpyIncDir, kcoreIncDir]
libraries = ["intersector", "kcore"]
ADDITIONALCPPFLAGS = []
if mpi:
    libraryDirs.append(mpiLibDir)
    includeDirs.append(mpiIncDir)
    ADDITIONALCPPFLAGS += ['-D_MPI']
if mpi4py:
    includeDirs.append(mpi4pyIncDir)

if mpi: libraries += mpiLibs

(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# setup ======================================================================
setup(
    name="Intersector",
    version="4.1",
    description="Mesh-intersection-based services in *Cassiopee*.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Intersector'],
    package_dir={"":"."},
    ext_modules=[Extension('Intersector.intersector',
                           sources=["Intersector/intersector.cpp"],
                           include_dirs=["Intersector"]+additionalIncludePaths+includeDirs,
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
                           extra_link_args=Dist.getLinkArgs()
                           )]
)

