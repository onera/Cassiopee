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

# Write KCore installation path to installPath.py
import KCore.Dist as Dist
Dist.writeSetupCfg()

additionalLibPaths = Dist.getAdditionalLibPaths()
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibs = Dist.getAdditionalLibs()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraries path =====================================================
libraryDirs = ["build/"+prod, kcoreLibDir, "../portaudio/lib"]
libraries = ['portaudio']
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# suppress --static
if prod == 'win64' or prod == 'msys64p3' or prod == 'msys64':
    inp = Dist.getLinkArgs()
    linkArgs = []
    for i in inp:
        if i != '--static': linkArgs.append(i) 
else:
    linkArgs = Dist.getLinkArgs()
    
# setup =======================================================================
listExtensions = []
listExtensions.append(
    Extension('pyaudio._portaudio',
             sources=['pyaudio/_portaudiomodule.c'],
             include_dirs=["pyaudio", "../portaudio/include"]+additionalIncludePaths+[numpyIncDir],
             library_dirs=additionalLibPaths+libraryDirs,
             libraries=libraries+additionalLibs,
             extra_compile_args=Dist.getCppArgs(),
             extra_link_args=linkArgs
	) )

# Setup ======================================================================
setup(
    name="pyaudio",
    version="2.9",
    description="Audio binding.",
    author="",
    package_dir={"":"."},
    packages=['pyaudio'],
    ext_modules=listExtensions
    )

# Check PYTHONPATH
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
