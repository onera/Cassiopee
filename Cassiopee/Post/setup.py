#=============================================================================
# Post requires:
# C++ compiler
# Fortran compiler
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

prod = os.getenv("ELSAPROD") or "xx"

# Setting includeDirs, libraryDirs and libraries ===========================================
ADDITIONALCPPFLAGS = []
includeDirs = [numpyIncDir, kcoreIncDir]
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["post", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

adolc, adolcIncDir, adolcLibDir, adolcLib = Dist.checkAdolc()
if adolc:
    includeDirs += [adolcIncDir]
    ADDITIONALCPPFLAGS = ["-DE_ADOLC"]
    libraryDirs += [adolcLibDir]
    libraries.append(adolcLib)

# extensions =================================================================
listExtensions = []
listExtensions.append(
    Extension('Post.post',
              sources=["Post/post.cpp"],
              include_dirs=["Post"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
              extra_link_args=Dist.getLinkArgs()+ADDITIONALCPPFLAGS
              ) )

# setup ======================================================================
setup(
    name="Post",
    version="4.2",
    description="Post-processing of CFD solutions.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Post'],
    package_dir={"":"."},
    ext_modules=listExtensions
)
