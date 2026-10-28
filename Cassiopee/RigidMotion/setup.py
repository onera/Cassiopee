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

# Write build info to file
Dist.writeBuildInfo("RigidMotion/RigidMotion.py")

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["rigidMotion", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

includeDirs=[numpyIncDir, kcoreIncDir]

# Extensions ==================================================================
listExtensions = [
    Extension(
        'RigidMotion.rigidMotion',
        sources=['RigidMotion/rigidMotion.cpp'],
        include_dirs=["RigidMotion"]+additionalIncludePaths+includeDirs,
        library_dirs=additionalLibPaths+libraryDirs,
        libraries=libraries+additionalLibs,
        extra_compile_args=Dist.getCppArgs(),
        extra_link_args=Dist.getLinkArgs()
    )
]

# Setup ======================================================================
setup(
    ext_modules=listExtensions
)
