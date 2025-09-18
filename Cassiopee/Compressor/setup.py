#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#=============================================================================
# Compressor requires:
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
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'
libraryDirs = ['build/'+prod, kcoreLibDir]
libraries = ["compressor", "kcore"]
from KCore.config import *
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

import srcs

# Extensions =================================================================
extensions = [
    Extension('Compressor.compressor',
              sources=["Compressor/compressor.cpp"],
              include_dirs=["Compressor"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+["kzstd"]+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs())]
if srcs.SZ:
    mySystem = Dist.getSystem()
    if mySystem[0] == 'mingw':
        additionalLibs += ["z"] # sometimes zlib1
    else: additionalLibs += ["z"]

    extensions += [
        Extension('Compressor.sz.csz',
                  sources=["Compressor/sz/compressor.cpp"],
                  include_dirs=["Compressor", "Compressor/sz/include"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
                  library_dirs=additionalLibPaths+libraryDirs,
                  libraries=libraries+["ksz", "kzstd"]+additionalLibs,
                  extra_compile_args=Dist.getCppArgs(),
                  extra_link_args=Dist.getLinkArgs())]
if srcs.ZFP:
    extensions += [
        Extension('Compressor.zfp.czfp',
                  sources=["Compressor/zfp/compressor.cpp"],
                  include_dirs=["Compressor", "Compressor/zfp/include"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
                  library_dirs=additionalLibPaths+libraryDirs,
                  libraries=libraries+["kzfp"]+additionalLibs,
                  extra_compile_args=Dist.getCppArgs(),
                  extra_link_args=Dist.getLinkArgs())]

# Setup ======================================================================
setup(
    name="Compressor",
    version="4.1",
    description="Compress CFD solutions.",
    author="ONERA",
    url="https://cassiopee.onera.fr",
    package_dir={"":"."},
    packages=['Compressor', 'Compressor.sz', 'Compressor.zfp'],
    ext_modules=extensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
