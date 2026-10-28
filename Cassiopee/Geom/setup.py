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

# Write build info to file
Dist.writeBuildInfo("Geom/Geom.py")

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")

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
listExtensions = [
    Extension(
        'Geom.geom',
        sources=["Geom/geom.cpp"],
        include_dirs=["Geom"]+additionalIncludePaths+includeDirs,
        library_dirs=additionalLibPaths+libraryDirs,
        libraries=libraries+additionalLibs,
        extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
        extra_link_args=Dist.getLinkArgs()+ADDITIONALCPPFLAGS
    )
]

setup(
    ext_modules=listExtensions
)
