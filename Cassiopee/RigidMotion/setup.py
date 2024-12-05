import os
from distutils.core import setup, Extension
#from setuptools import setup, Extension

#=============================================================================
# RigidMotion requires:
# C++ compiler
# Fortran compiler: defined in config.py
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

# Compilation des fortrans ===================================================
from KCore.config import *
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["rigidMotion", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

includeDirs=[numpyIncDir, kcoreIncDir]
# Extensions ==================================================================
listExtensions = []
listExtensions.append(
    Extension('RigidMotion.rigidMotion',
              sources=['RigidMotion/rigidMotion.cpp'],
              include_dirs=["RigidMotion"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# Setup ======================================================================
setup(
    name="RigidMotion",
    version="4.0",
    description="Compute/define rigid motion.",
    author="ONERA",
    url="https://cassiopee.onera.fr",
    packages=['RigidMotion'],
    package_dir={"":"."},
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
