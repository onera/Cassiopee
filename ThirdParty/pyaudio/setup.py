#!/usr/bin/env python

import os, sys
from distutils.core import setup, Extension

#=============================================================================
# KCore requires:
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
#=============================================================================
# Compiler settings must be set in config.py
from distutils.core import setup, Extension
import os, sys

# Write KCore installation path to installPath.py
import KCore.Dist as Dist
Dist.writeSetupCfg()
from KCore.config import *

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraries path =====================================================
libraryDirs = ["build/"+prod, kcoreLibDir, "../portaudio/lib"]
libraries = ['portaudio']
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# setup =======================================================================
listExtensions = []
listExtensions.append(
    Extension('pyaudio._portaudio',
             sources=['pyaudio/_portaudiomodule.c'],
             include_dirs=["pyaudio", "../portaudio/include"]+additionalIncludePaths+[numpyIncDir],
             library_dirs=additionalLibPaths+libraryDirs,
             libraries=libraries+additionalLibs,
             extra_compile_args=Dist.getCppArgs(),
             extra_link_args=Dist.getLinkArgs()
	) )

# Setup ======================================================================
setup(
    name="pyaudio",
    version="2.5",
    description="Audio binding.",
    author="",
    package_dir={"":"."},
    packages=['pyaudio'],
    ext_modules=listExtensions
    )

# Check PYTHONPATH
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
