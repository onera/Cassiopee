#=============================================================================
# Generator requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler
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
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["generator", "generator2", "generator3", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# Extensions =================================================================
listExtensions = []
listExtensions.append(
    Extension('Generator.generator',
              sources=["Generator/generator.cpp"],
              include_dirs=["Generator"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup =======================================================================
setup(
    name="Generator",
    version="4.1",
    description="*Cassiopee* module of mesh generation.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Generator'],
    package_dir={"":"."},
    ext_modules=listExtensions
)

