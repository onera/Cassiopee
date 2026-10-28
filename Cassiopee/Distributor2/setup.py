#=============================================================================
# Distributor2 requires:
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

# Write build info to file
Dist.writeBuildInfo("Distributor2/Distributor2.py")

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Setting libraryDirs and libraries ===========================================
prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")

libraryDirs = ['build/'+prod, kcoreLibDir]
libraries = ["distributor2", "kcore"]
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs

# setup ======================================================================
listExtensions = [
    Extension(
        'Distributor2.distributor2',
        sources=['Distributor2/distributor2.cpp'],
        include_dirs=["Distributor2"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
        library_dirs=additionalLibPaths+libraryDirs,
        libraries=libraries+additionalLibs,
        extra_compile_args=Dist.getCppArgs(),
        extra_link_args=Dist.getLinkArgs()
    )
]

setup(
    ext_modules=listExtensions
)
