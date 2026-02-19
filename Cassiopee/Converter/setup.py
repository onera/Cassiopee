#=============================================================================
# Converter requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore library
#=============================================================================
import os
from setuptools import setup, Extension
from importlib.util import spec_from_file_location, module_from_spec
import KCore.Dist as Dist

def loadModuleFromPath(modname):
    # Load a Python file by filesystem path (PEP-517 isolated build requirement)
    helper = os.path.join(os.path.dirname(__file__), modname + ".py")
    spec = spec_from_file_location(modname, helper)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

additionalLibPaths = Dist.getAdditionalLibPaths()
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Test if libhdf5 exists ======================================================
(hdf, hdfIncDir, hdfLibDir, hdflibs) = Dist.checkHdf()

# Test if libnetcdf exists ======================================================
(netcdf, netcdfIncDir, netcdfLibDir, netcdflibs) = Dist.checkNetcdf()

# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi()
(mpi4py, mpi4pyIncDir, mpi4pyLibDir) = Dist.checkMpi4py()

# Compilation des fortrans ====================================================
prod = os.getenv("ELSAPROD") or 'xx'

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
    ADDITIONALCPPFLAGS += ['-D_MPI4PY']

if hdf:
    for l in hdflibs: libraries.append(l)
if netcdf:
    for l in netcdflibs: libraries.append(l)

if mpi: libraries += mpiLibs

(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
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
srcs = loadModuleFromPath('srcs')
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
    version="4.1",
    description="Converter for *Cassiopee* modules.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Converter'],
    package_dir={"":"."},
    ext_modules=listExtensions
)

