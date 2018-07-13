#!/usr/bin/env python
from distutils.core import setup, Extension
import os, sys

#=============================================================================
# Distributor requires:
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

# Is Cassiopee present? ======================================================
(Cassiopee, CassiopeeIncDir, CassiopeeLibDir, CassiopeeUseMpi) = Dist.checkCassiopee()
if (Cassiopee == False):
    print "Error: Cassiopee/Kernel is required for the compilation of Distributor."
    sys.exit()

# Setting libraryDirs and libraries ===========================================
libraryDirs = [kcoreLibDir]
if (Cassiopee == True): libraryDirs.append(CassiopeeLibDir)
libraries = ["kcore", "elsAc"]
from KCore.config import *
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths ; libraries += libs
if (CassiopeeUseMpi == True): libraries += ["mpi"]
    
# setup ======================================================================
setup(
    name="Distributor",
    version="2.6",
    description="Distributor for *Cassiopee* modules.",
    author="S. Peron, C. Benoit, G. Jeanfaivre, P. Raud",
    package_dir={"":"."},
    packages=['Distributor'],
    ext_modules=[Extension('Distributor.distributor',
                           sources=['Distributor/distributor.cpp',
                                    'Distributor/distributor1.cpp',
                                    'Distributor/distributor2.cpp'
                                    ],
                           include_dirs=["Distributor"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir, CassiopeeIncDir],
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )
                 ]
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
