#=============================================================================
# Post requires:
# C++ compiler
# Fortran compiler
# Numpy
# KCore
#=============================================================================
import os
from setuptools import setup, Extension
from importlib.util import spec_from_file_location, module_from_spec
import KCore.Dist as Dist

def loadModuleFromPath(modname):
    # Load a Python file by filesystem path (PEP-517 isolated build requirement)
    helper = os.path.join(os.path.dirname(__file__), modname + ".py")
    spec = spec_from_file_location(modname, helper)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

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

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["post", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

if Dist.ADOLC:
    (adolc, adolcIncDir, adolcLibDir, adolcLib) = Dist.checkAdolc()
    if adolc:
        libraryDirs += adolcLibDir; libraries += [adolcLib]

srcs = loadModuleFromPath('srcs')

# extensions =================================================================
listExtensions = []
listExtensions.append(
    Extension('Post.post',
              sources=["Post/post.cpp"],
              include_dirs=["Post"]+additionalIncludePaths+[numpyIncDir,kcoreIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Post",
    version="4.1",
    description="Post-processing of CFD solutions.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Post'],
    package_dir={"":"."},
    ext_modules=listExtensions
)

