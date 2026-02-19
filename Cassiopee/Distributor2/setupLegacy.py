#from distutils.core import setup, Extension
from setuptools import setup, Extension

#=============================================================================
# Distributor2 requires:
# C++ compiler
# Numpy
# KCore
#=============================================================================

import KCore.Dist as Dist

# Compiler settings must be set in installBase.py / installBaseUser.py
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibPaths = Dist.getAdditionalLibPaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Setting libraryDirs and libraries ===========================================
libraryDirs = [kcoreLibDir]
libraries = ["kcore"]
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# setup ======================================================================
import srcs
setup(
    name="Distributor2",
    version="4.1",
    description="Distributor for arrays and pyTrees.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    package_dir={"":"."},
    packages=['Distributor2'],
    ext_modules=[Extension('Distributor2.distributor2',
                           sources=['Distributor2/distributor2.cpp']+srcs.cpp_srcs,
                           include_dirs=["Distributor2"]+additionalIncludePaths+[numpyIncDir,kcoreIncDir],
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )
                 ]
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
