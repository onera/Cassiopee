#=============================================================================
# RigidMotion requires:
# C++ compiler
# Fortran compiler: defined in config.py
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

# Compilation des fortrans ===================================================
prod = os.getenv("ELSAPROD") or "xx"

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["rigidMotion", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

includeDirs=[numpyIncDir, kcoreIncDir]
# Extensions ==================================================================
listExtensions = []
listExtensions.append(
    Extension('RigidMotion.rigidMotion',
              sources=['RigidMotion/rigidMotion.cpp'],
              include_dirs=["RigidMotion"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# Setup ======================================================================
setup(
    name="RigidMotion",
    version="4.1",
    description="Compute/define rigid motion.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['RigidMotion'],
    package_dir={"":"."},
    ext_modules=listExtensions
)

