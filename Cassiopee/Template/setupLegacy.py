#!/usr/bin/env python
from distutils.core import setup, Extension
import os, sys

#=============================================================================
# Template requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore
#=============================================================================

# Write setup.cfg file
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

from KCore.config import *

# Compilation des fortrans ====================================================
#if f77compiler == "None":
#    print("Error: a fortran 77 compiler is required for compiling Fast.")
#args = Dist.getForArgs(); opt = ''
#for c in xrange(len(args)):
#    opt += 'FOPT'+str(c)+'='+args[c]+' '
#os.system("make -e FC="+f77compiler+" WDIR=Template/Fortran "+opt)
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir]
includeDirs = [numpyIncDir, kcoreIncDir]
libraries = ["kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Extensions ==================================================================
import srcs
listExtensions = []
listExtensions.append(
    Extension('Template.template',
              sources=['Template/template.cpp']+srcs.cpp_srcs,
              include_dirs=["Template"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Template",
    version="2.0",
    description="Template module.",
    author="You",
    package_dir={"":"."},
    packages=['Template'],
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
