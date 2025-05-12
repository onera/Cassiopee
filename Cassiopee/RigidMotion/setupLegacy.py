#!/usr/bin/env python
import os
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

# Compilation des fortrans ===================================================
from KCore.config import *
if f77compiler == "None":
    print("Error: a fortran 77 compiler is required for compiling RigidMotion.")
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT'+str(c)+'='+v+' '
os.system("make -e FC="+f77compiler+" WDIR=RigidMotion/Fortran "+opt)
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod,kcoreLibDir]
libraries = ["RigidMotionF","kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Extensions =================================================================
import srcs
setup(
    name="RigidMotion",
    version="4.1",
    description="Rigid motion module.",
    author="Onera",
    package_dir={"":"."},
    packages=['RigidMotion'],
    ext_modules=[Extension('RigidMotion.rigidMotion',
                           sources=["RigidMotion/rigidMotion.cpp"]+srcs.cpp_srcs,
                           include_dirs=["RigidMotion"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )]
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
