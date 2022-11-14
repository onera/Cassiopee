from distutils.core import setup, Extension
#from setuptools import setup, Extension
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

# Python ======================================================================
(pythonVersion, pythonIncDir, pythonLibDir, pythonLibs) = Dist.checkPython()

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
# pdm supprime ici
libraries = ["xcore", "scotch1", "scotch2", "pdm", "kcore"]

mySystem = Dist.getSystem()
if mySystem[0] == 'mingw': 
  libraries += ["wsock32"]

ADDITIONALCPPFLAGS = ['-fpermissive']
if mpi:
    libraryDirs.append(mpiLibDir)
    includeDirs.append(mpiIncDir)
    ADDITIONALCPPFLAGS += ['-D_MPI']
    libraries += ["ptscotch", "scotch1", "scotch2", "scotch1"]
if mpi4py:
    includeDirs.append(mpi4pyIncDir)
if mpi: libraries += mpiLibs
        
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
    
# Extensions ==================================================================
listExtensions = []
listExtensions.append(
    Extension('XCore.xcore',
              sources=['XCore/xcore.cpp'],
              include_dirs=["XCore"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
              extra_link_args=Dist.getLinkArgs()
              ) )

listExtensionsPyx = []
import srcs_paradigma
from Cython.Build import cythonize
for c in srcs_paradigma.pyx_srcs:
    name = c.replace('.pyx', '')
    names = name.split('/')
    name = names[0]+'.'+names[-1]
    listExtensionsPyx.append(
        Extension(name,
                  sources=[c],
                  include_dirs=["XCore","XCore/paradigma","XCore/paradigma/ppart","XCore/paradigma/struct","XCore/paradigma/pario","XCore/paradigma/mesh"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir, mpiIncDir, mpi4pyIncDir, pythonIncDir],
                  library_dirs=additionalLibPaths+libraryDirs,
                  libraries=libraries+additionalLibs,
                  extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
                  extra_link_args=[],
                  language='c++'
                  ) )

# setup ======================================================================
setup(
    name="XCore",
    version="3.5",
    description="XCore for *Cassiopee* modules.",
    author="ONERA",
    url="http://elsa.onera.fr/Cassiopee",
    packages=['XCore'],
    package_dir={"":"."},
    ext_modules=listExtensions+cythonize(listExtensionsPyx,include_path=["XCore/paradigma"])
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
