#from distutils.core import setup, Extension
from setuptools import setup, Extension

#=============================================================================
# Compressor requires:
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
    version="4.1",
    description="Compress CFD solutions.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    package_dir={"":"."},
    #packages=['Compressor', 'Compressor.zfp', 'Compressor.sz'],
    packages=['Compressor'],
    ext_modules=extensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
