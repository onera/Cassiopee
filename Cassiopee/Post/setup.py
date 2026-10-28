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

# Write build info to file
Dist.writeBuildInfo("Post/Post.py")

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")

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
listExtensions = [
    Extension(
        'Post.post',
        sources=["Post/post.cpp"],
        include_dirs=["Post"]+additionalIncludePaths+includeDirs,
        library_dirs=additionalLibPaths+libraryDirs,
        libraries=libraries+additionalLibs,
        extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
        extra_link_args=Dist.getLinkArgs()+ADDITIONALCPPFLAGS
    )
]

# setup ======================================================================
setup(
    ext_modules=listExtensions
)
