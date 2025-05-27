from distutils.core import setup, Extension
#from setuptools import setup, Extension
from KCore.config import *
import os

#=============================================================================
# OCC requires:
# ELSAPROD variable defined in environment
# C++ compiler
# KCore library
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Test if generator exists ===================================================
(generatorVersion, generatorIncDir, generatorLibDir) = Dist.checkGenerator()

# Test if open-cascade is already installed ==================================
(OCEPresent, OCEIncDir, OCELibDir) = Dist.checkOCE(additionalLibPaths,
                                                   additionalIncludePaths)

if not OCEPresent: os._exit(0)

# Compilation des fortrans ===================================================
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir, generatorLibDir]
includeDirs = [numpyIncDir, kcoreIncDir, generatorIncDir]
#libraries = ["occ_cassiopee", "generator", "converter", "kcore"]
libraries = ["occ_cassiopee", "generator", "kcore"]

if OCEPresent:
    libraryDirs += [OCELibDir]
    includeDirs += [OCEIncDir]

import srcs
libOCE = srcs.allMods
if OCEPresent and Dist.getSystem()[0] == 'mingw':
    libOCE = [i+".dll" for i in libOCE]
libraries += libOCE + libOCE

(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# setup ======================================================================
setup(
    name="OCC",
    version="4.1",
    description="OpenCascade python module.",
    author="ONERA",
    url="https://cassiopee.onera.fr",
    packages=['OCC'],
    package_dir={"":"."},
    ext_modules=[Extension('OCC.occ',
                           sources=["OCC/occ.cpp"],
                           include_dirs=["OCC"]+additionalIncludePaths+includeDirs,
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )]
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
