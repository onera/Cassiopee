#=============================================================================
# Dist2Walls requires:
# C++ compiler
# Numpy
# KCore
#=============================================================================
import os
from setuptools import setup, Extension
import KCore.Dist as Dist

additionalLibPaths = Dist.getAdditionalLibPaths()
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Setting libraryDirs and libraries ===========================================
prod = os.getenv("ELSAPROD") or "xx"
libraryDirs = [kcoreLibDir, 'build/'+prod]
libraries = ["dist2walls", "kcore"]

(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkFortranLibs()
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
    url="https://onera.github.io/Cassiopee/",
    packages=['Dist2Walls'],
    package_dir={"":"."},
    ext_modules=extensions
)

