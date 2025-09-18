#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#=============================================================================
# Connector requires:
# C++ compiler
# Fortran compiler
# Numpy
# KCore
#=============================================================================

# Write setup.cfg file
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

from KCore.config import *

# Compilation des fortrans ===================================================
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["connector", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
includeDirs = [numpyIncDir, kcoreIncDir]
ADDITIONALCPPFLAGS=[]

# setup =======================================================================
listExtensions = []
listExtensions.append(
    Extension('Connector.connector',
              sources=['Connector/connector.cpp'],
              include_dirs=["Connector"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Connector",
    version="4.1",
    description="Connector for *Cassiopee* modules.",
    author="ONERA",
    url="https://cassiopee.onera.fr",
    packages=['Connector'],
    package_dir={"":"."},
    ext_modules=listExtensions)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
