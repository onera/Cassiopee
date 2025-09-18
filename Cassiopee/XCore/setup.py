#from distutils.core import setup, Extension
from setuptools import setup, Extension
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

# Get compilers from Distutils ================================================
(cc, cxx, opt, basecflags, ccshared, ldshared, so_ext) = Dist.getDistUtilsCompilers()

# Python ======================================================================
(pythonVersion, pythonIncDir, pythonLibDir, pythonLibs) = Dist.checkPython()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Test if fortran exists
(ok, fortranLibs, fortranLibDir) = Dist.checkFortranLibs([], additionalLibPaths)

# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi(additionalLibPaths,
                                                     additionalIncludePaths)
(mpi4py, mpi4pyIncDir, mpi4pyLibDir) = Dist.checkMpi4py(additionalLibPaths,
                                                        additionalIncludePaths)

prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir]+fortranLibDir
includeDirs = [numpyIncDir, kcoreIncDir]

import srcs
libraries = ["xcore"]
if srcs.ZOLTAN and mpi: libraries += ["zoltan"]
if srcs.SCOTCH: libraries += ["scotch1", "scotch2"]
if srcs.PARADIGMA != 0: libraries += ["pdm"]
libraries += ["kcore"]
libraries += fortranLibs

mySystem = Dist.getSystem()
if mySystem[0] == 'mingw':
    libraries += ["wsock32"]

if cc == 'icc':
    ADDITIONALCPPFLAGS = ['-fpermissive']
elif cc == 'gcc':
    ADDITIONALCPPFLAGS = ['-fpermissive']
else: ADDITIONALCPPFLAGS = []

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
cython = Dist.checkCython(additionalLibPaths, additionalIncludePaths)

if cython and mpi and mpi4py:
    if srcs.PARADIGMA == 2:
        import srcs_paradigma23 as srcs_paradigma
        paradigmaDir = "XCore/paradgima"
    elif srcs.PARADIGMA == 1:
        import srcs_paradigma
        paradigmaDir = "XCore/paradigma23"
    if srcs.PARADIGMA != 0:
        PPATH = srcs_paradigma.PPATH
        from Cython.Build import cythonize
        for c in srcs_paradigma.pyx_srcs:
            name = c.replace('.pyx', '')
            names = name.split('/')
            name = names[0]+'.'+names[-1]
            listExtensionsPyx.append(
                Extension(name,
                          sources=[c],
                          include_dirs=["XCore","%s"%PPATH,"%s/ppart"%PPATH,"%s/struct"%PPATH,"%s/pario"%PPATH,"%s/mesh"%PPATH,"%s/meshgen"%PPATH,"%s/mpi_wrapper"%PPATH, "%s/util"%PPATH]+additionalIncludePaths+[numpyIncDir, kcoreIncDir, mpiIncDir, mpi4pyIncDir, pythonIncDir],
                          library_dirs=additionalLibPaths+libraryDirs,
                          libraries=libraries+additionalLibs,
                          extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
                          extra_link_args=[],
                          language='c++'
                          ) )
    else:
        def cythonize(srcs, include_path): return []
        paradigmaDir = ''
else:
    def cythonize(srcs, include_path): return []
    paradigmaDir = ''

# setup ======================================================================
setup(
    name="XCore",
    version="4.1",
    description="XCore for *Cassiopee* modules.",
    author="ONERA",
    url="https://cassiopee.onera.fr",
    packages=['XCore'],
    package_dir={"":"."},
    ext_modules=listExtensions+cythonize(listExtensionsPyx,include_path=[paradigmaDir])
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
