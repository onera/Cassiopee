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

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Setting libraryDirs and libraries ===========================================
prod = os.getenv("ELSAPROD") or "xx"

libraryDirs = ['build/'+prod, kcoreLibDir]
libraries = ["distributor2", "kcore"]
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs

# setup ======================================================================
setup(
    name="Distributor2",
    version="4.1",
    description="Distributor for arrays and pyTrees.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Distributor2'],
    package_dir={"":"."},
    ext_modules=[Extension('Distributor2.distributor2',
                           sources=['Distributor2/distributor2.cpp'],
                           include_dirs=["Distributor2"]+additionalIncludePaths+[numpyIncDir,kcoreIncDir],
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )
                 ]
)

