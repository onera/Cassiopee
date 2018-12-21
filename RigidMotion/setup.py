#!/usr/bin/env python
import os, sys
from distutils.core import setup, Extension

#=============================================================================
# RigidMotion requires:
# C++ compiler
# Numpy
# KCore
#
# Optional motion from solvers requires:
# Cassiopee/Kernel
# elsA/Kernel
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Setting libraryDirs and libraries ===========================================
libraryDirs = [kcoreLibDir]
libraries = ["kcore"] 
from KCore.config import *
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

includeDirs = [numpyIncDir, kcoreIncDir]
    
# Extensions =================================================================
import srcs
extensions = [
    Extension('RigidMotion.rigidMotion',
              sources=["RigidMotion/rigidMotion.cpp"]+srcs.cpp_srcs,
              include_dirs=["RigidMotion"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
	)
    ]

# Setup ======================================================================
setup(
    name="RigidMotion",
    version="2.9",
    description="Compute/define rigid motions.",
    author="Onera",
    package_dir={"":"."},
    packages=['RigidMotion'],
    ext_modules=extensions
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
