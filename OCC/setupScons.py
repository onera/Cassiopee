#!/usr/bin/env python
from distutils.core import setup, Extension
import os, sys

#=============================================================================
# OCC requires:
# ELSAPROD variable defined in environment
# C++ compiler
# KCore library
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Test if generator exists ===================================================
(generatorVersion, generatorIncDir, generatorLibDir) = Dist.checkGenerator()
    
# Compilation des fortrans ===================================================
from KCore.config import *
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir, generatorLibDir]
includeDirs = [kcoreIncDir, generatorIncDir]
libraries = ["occ_cassiopee", 
"TKIGES", "TKXSBase", "TKShHealing", "TKTopAlgo", "TKPrim", "TKBool", "TKBool2", "TKBool3", "TKBool4", "TKGeomAlgo", "TKBRep", "TKBRep2", "TKGeomBase", "TKG3d", "TKMath", "TKernel", "TKG2d",
"TKIGES", "TKXSBase", "TKShHealing", "TKTopAlgo", "TKPrim", "TKBool", "TKBool2", "TKBool3", "TKBool4", "TKGeomAlgo", "TKBRep", "TKBRep2", "TKGeomBase", "TKG3d", "TKMath", "TKernel", "TKG2d", "generator", "kcore"]

(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# setup ======================================================================
setup(
    name="OCC",
    version="2.7",
    description="OpenCascade.",
    author="Onera",
    package_dir={"":"."},
    packages=['OCC'],
    ext_modules=[Extension('OCC.occ',
                           sources=["OCC/occ.cpp"],
                           include_dirs=["OCC"]+additionalIncludePaths+includeDirs,
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs() + ['-lrt']
                           )]
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
