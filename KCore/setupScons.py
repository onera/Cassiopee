#!/usr/bin/env python

import os
from distutils.core import setup, Extension

#=============================================================================
# KCore requires:
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# Scons
#=============================================================================
# Compiler settings must be set in config.py
from config import *

# Write KCore installation path to installPath.py
import Dist
Dist.writeInstallPath()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraries path =====================================================
libraryDirs = ["build/"+prod]
libraries = ["kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Extensions =================================================================
srcs = ['KCore/kcore.cpp']
extensions = [
    Extension('KCore.kcore',
              sources=srcs,
              include_dirs=["KCore"]+additionalIncludePaths+[numpyIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCArgs(),
              extra_link_args=Dist.getLinkArgs()
	)
    ]

# Setup ======================================================================
setup(
    name="KCore",
    version="2.9",
    description="Core for *Cassiopee* modules.",
    author="Onera",
    package_dir={"":"."},
    packages=['KCore'],
    ext_modules=extensions
    )

# Check PYTHONPATH
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
