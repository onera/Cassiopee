#!/usr/bin/env python
from distutils.core import setup, Extension
import os, sys

#=============================================================================
# Post requires:
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

# Compilation des fortrans ===================================================
from KCore.config import *
if f77compiler == "None":
    print("Error: a fortran 77 compiler is required for compiling Post.")
    sys.exit()
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT'+str(c)+'='+v+' '
os.system("make -e FC="+f77compiler+" F90=true WDIR=Post/Fortran "+opt)
os.system("make -e FC="+f77compiler+" F90=true WDIR=Post/zipper "+opt)
if f90compiler != "None" and os.access('Post/usurp', os.F_OK):
    os.system("(cd Post/usurp; make -e FC="+f77compiler+" F90="+f90compiler+" "+opt+")")
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["PostF", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

if f90compiler != "None" and os.access('Post/usurp', os.F_OK): libraries.append("UsurpF")

import srcs

# extensions =================================================================
listExtensions = []
listExtensions.append(
    Extension('Post.post',
              sources=["Post/post.cpp"]+srcs.cpp_srcs,
              include_dirs=["Post"]+additionalIncludePaths+[numpyIncDir,kcoreIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Post",
    version="4.0",
    description="Post-processing of CFD solutions.",
    author="Onera",
    package_dir={"":"."},
    packages=['Post'],
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
