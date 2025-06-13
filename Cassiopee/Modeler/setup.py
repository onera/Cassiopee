from distutils.core import setup, Extension
#from setuptools import setup
import os

#=============================================================================
# Modeler requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore, Converter, Generator, Transform
#=============================================================================

# Write setup.cfg file
import KCore.Dist as Dist
Dist.writeSetupCfg()

from KCore.config import *

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Test if open-cascade is already installed ==================================
(OCCPresent, OCCIncDir, OCCLibDir) = Dist.checkOCC(additionalLibPaths,
                                                   additionalIncludePaths)

import srcs
if srcs.TIGL and not OCCPresent:
    print("Warning: open cascade not found on your system. Modeler.Tigl not installed.")
    os._exit(0)

prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir]
includeDirs = [numpyIncDir, kcoreIncDir]
libraries = ["kcore", "modeler", "modeler"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

if OCCPresent:
    libraryDirs += [OCCLibDir]
    includeDirs += [OCCIncDir]

libOCC = Dist.getOCCModules()
if OCCPresent and Dist.getSystem()[0] == 'mingw':
    libOCC = [i+".dll" for i in libOCC]
libraries += libOCC + libOCC

import srcs
if srcs.TIXI: libraries += ["curl", "xml2", "xslt"]
if srcs.TIGL: libraries += ["boost_system", "boost_filesystem", "boost_date_time"]

# Extensions ==================================================================
listExtensions = []
listExtensions.append(
    Extension('Modeler.modeler',
              sources=['Modeler/modeler.cpp'],
              include_dirs=["Modeler"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Modeler",
    version="4.1",
    description="Modeler module.",
    author="ONERA",
    url="https://cassiopee.onera.fr",
    packages=['Modeler', 'Modeler.CPACS'],
    package_dir={"":"."},
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
