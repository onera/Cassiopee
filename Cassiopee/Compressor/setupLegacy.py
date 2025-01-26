#!/usr/bin/env python
from distutils.core import setup, Extension

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
libraryDirs = [kcoreLibDir]
libraries = ["kcore"]
from KCore.config import *
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

includeDirs = [numpyIncDir, kcoreIncDir]

# Extensions =================================================================
import srcs
extensions = [
    Extension('Compressor.compressor',
              sources=["Compressor/compressor.cpp"]+srcs.cpp_srcs,
              include_dirs=["Compressor"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()),
    #Extension('Compressor.sz.csz',
    #          sources=["Compressor/sz/compressor.cpp"],
    #          include_dirs=["Compressor"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
    #          library_dirs=additionalLibPaths+libraryDirs,
    #          libraries=libraries+["SZ",]+additionalLibs,
    #          extra_compile_args=Dist.getCppArgs(),
    #          extra_link_args=Dist.getLinkArgs()),
    #Extension('Compressor.zfp.czfp',
    #          sources=["Compressor/zfp/compressor.cpp"],
    #          include_dirs=["Compressor"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
    #          library_dirs=additionalLibPaths+libraryDirs,
    #          libraries=libraries+["zfp",]+additionalLibs,
    #          extra_compile_args=Dist.getCppArgs(),
    #          extra_link_args=Dist.getLinkArgs()),
]

# Setup ======================================================================
setup(
    name="Compressor",
    version="4.0",
    description="Compress CFD solutions.",
    author="Onera",
    package_dir={"":"."},
    #packages=['Compressor', 'Compressor.zfp', 'Compressor.sz'],
    packages=['Compressor'],
    ext_modules=extensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
