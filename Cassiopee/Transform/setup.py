#=============================================================================
# Transform requires:
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
Dist.writeBuildInfo("Transform/Transform.py")

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLib) = Dist.checkModuleCassiopee("KCore")

prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLib]
includeDirs = [numpyIncDir, kcoreIncDir]
libraries = ["transform", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# Extensions ==================================================================
listExtensions = [
    Extension(
        'Transform.transform',
        sources=['Transform/transform.cpp'],
        include_dirs=["Transform"]+additionalIncludePaths+includeDirs,
        library_dirs=additionalLibPaths+libraryDirs,
        libraries=libraries+additionalLibs,
        extra_compile_args=Dist.getCppArgs(),
        extra_link_args=Dist.getLinkArgs()
    )
]

# setup ======================================================================
setup(
    ext_modules=listExtensions
)
