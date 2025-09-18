#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#=============================================================================
# Dist2Walls requires:
# C++ compiler
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

# Setting libraryDirs and libraries ===========================================
from KCore.config import *
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'
libraryDirs = [kcoreLibDir, 'build/'+prod]
libraries = ["dist2walls", "kcore"]
from KCore.config import *
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Extensions =================================================================
extensions = [
    Extension('Dist2Walls.dist2walls',
              sources=["Dist2Walls/dist2walls.cpp"],
              include_dirs=["Dist2Walls"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              )
]

# Setup ======================================================================
setup(
    name="Dist2Walls",
    version="4.1",
    description="Computation of distance to walls.",
    author="ONERA",
    url="https://cassiopee.onera.fr",
    packages=['Dist2Walls'],
    package_dir={"":"."},
    ext_modules=extensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
