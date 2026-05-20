#=============================================================================
# Geom requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore library
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

# Setting includeDirs, libraryDirs and libraries ===========================================
ADDITIONALCPPFLAGS = []
includeDirs = [numpyIncDir, kcoreIncDir]
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["geom", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

adolc, adolcIncDir, adolcLibDir, adolcLib = Dist.checkAdolc()
if adolc:
    includeDirs += [adolcIncDir]
    ADDITIONALCPPFLAGS = ["-DE_ADOLC"]
    libraryDirs += [adolcLibDir]
    libraries.append(adolcLib)

# setup ======================================================================
setup(
    name="Geom",
    version="4.2",
    description="Geometry definition for *Cassiopee* modules.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Geom'],
    package_dir={"":"."},
    ext_modules=[Extension('Geom.geom',
                           sources=["Geom/geom.cpp"],
                           include_dirs=["Geom"]+additionalIncludePaths+includeDirs,
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
                           extra_link_args=Dist.getLinkArgs()+ADDITIONALCPPFLAGS
                           )]
)
