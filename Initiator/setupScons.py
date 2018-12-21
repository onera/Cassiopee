#!/usr/bin/env python
from distutils.core import setup, Extension
import os

#=============================================================================
# Initiator requires :
# [ENV] CASSIOPEE, ELSAPROD
# C++ compiler
# Fortran compiler : defined in config.py
# Numpy
# KCore library
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
libraries = ["initiator", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# setup =======================================================================
setup(
    name="Initiator",
    version="2.9",
    description="Initiator for *Cassiopee* modules.",
    author="Onera",
    package_dir={"":"."},
    packages=['Initiator'],
    ext_modules=[Extension('Initiator.initiator',
                           sources=['Initiator/initiator.cpp'],
                           include_dirs=["Initiator"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )
                 ]
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
