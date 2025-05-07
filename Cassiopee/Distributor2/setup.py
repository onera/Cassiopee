from distutils.core import setup, Extension
#from setuptools import setup, Extension
import os

#=============================================================================
# Distributor2 requires:
# C++ compiler
# Numpy
# KCore
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Setting libraryDirs and libraries ===========================================
from KCore.config import *
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

libraryDirs = ['build/'+prod, kcoreLibDir]
libraries = ["distributor2", "kcore"]
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# setup ======================================================================
setup(
    name="Distributor2",
    version="4.1",
    description="Distributor for arrays and pyTrees.",
    author="ONERA",
    url="https://cassiopee.onera.fr",
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

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
