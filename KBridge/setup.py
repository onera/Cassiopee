#!/usr/bin/env python
import os, sys
from distutils.core import setup, Extension

#=============================================================================
# KBridge requires:
# C++ compiler
# Numpy
# KCore
# Cassiopee/Kernel
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Is Cartesian solver present? ================================================
(Cassiopee, CassiopeeIncDir, CassiopeeLibDir, CassiopeeUseMpi) = Dist.checkCassiopee()
if (Cassiopee == False):
    print "Error: Cassiopee solver must be installed for compiling KBridge."
    sys.exit()

# Setting libraryDirs and libraries ===========================================
libraryDirs = [kcoreLibDir]
libraries = ["kcore"]
from KCore.config import *
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths ; libraries += libs

includeDirs = [numpyIncDir, kcoreIncDir]
includeDirs += [CassiopeeIncDir]
libraryDirs += [CassiopeeLibDir]
libraries += ["elsAc"]
if (CassiopeeUseMpi == True): libraries += ["mpi"]
    
# Extensions =================================================================
import srcs
extensions = [
    Extension('KBridge.kbridge',
              sources=["KBridge/kbridge.cpp"]+srcs.cpp_srcs,
              include_dirs=["KBridge"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
	)
    ]

# Setup ======================================================================
setup(
    name="KBridge",
    version="2.6",
    description="Bridge to Cassiopee solver.",
    author="Onera",
    package_dir={"":"."},
    packages=['KBridge'],
    ext_modules=extensions
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath() ; Dist.checkLdLibraryPath()
