#=============================================================================
# Initiator requires:
# [ENV] CASSIOPEE, ELSAPROD
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

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["initiator", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# setup =======================================================================
setup(
    name="Initiator",
    version="4.1",
    description="Initiator for *Cassiopee* modules.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Initiator'],
    package_dir={"":"."},
    ext_modules=[Extension('Initiator.initiator',
                           sources=['Initiator/initiator.cpp'],
                           include_dirs=["Initiator"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )
                 ]
)

