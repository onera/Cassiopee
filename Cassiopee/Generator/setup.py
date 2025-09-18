#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#=============================================================================
# Generator requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler
# Numpy
# KCore
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Setting libraryDirs and libraries ===========================================
from KCore.config import *
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["generator", "generator2", "generator3", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Extensions =================================================================
listExtensions = []
listExtensions.append(
    Extension('Generator.generator',
              sources=["Generator/generator.cpp"],
              include_dirs=["Generator"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup =======================================================================
setup(
    name="Generator",
    version="4.1",
    description="*Cassiopee* module of mesh generation.",
    author="ONERA",
    url="https://cassiopee.onera.fr",
    packages=['Generator'],
    package_dir={"":"."},
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
