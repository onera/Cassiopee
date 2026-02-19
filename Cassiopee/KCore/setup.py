#=============================================================================
# KCore requires:
# C++ compiler
# Fortran compiler
# Numpy
# Scons
#=============================================================================
import os
from setuptools import setup, Extension
from importlib.util import spec_from_file_location, module_from_spec

def loadModuleFromPath(modname):
    # Load a Python file by filesystem path (PEP-517 isolated build requirement)
    helper = os.path.join(os.path.dirname(__file__), modname + ".py")
    spec = spec_from_file_location(modname, helper)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

# Compiler settings must be set in installBase.py / installBaseUser.py
Dist = loadModuleFromPath('Dist')
installBase = loadModuleFromPath('installBase')
Dist.setConfigDict(installBase.installDict)
additionalLibPaths = Dist.getAdditionalLibPaths()
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists
numpyVersion, numpyIncDir, numpyLibDir = Dist.checkNumpy()

prod = os.getenv("ELSAPROD") or "xx"
# Setting libraries path
libraryDirs = ["build/" + prod]
libraries = ["kcore"]

ok, libs, paths = Dist.checkFortranLibs()
libraryDirs += paths
libraries += libs

ok, libs, paths = Dist.checkCppLibs()
libraryDirs += paths
libraries += libs

if Dist.ADOLC:
    adolc, adolcIncDir, adolcLibDir, adolcLib = Dist.checkAdolc()
    if adolc:
        libraryDirs += adolcLibDir
        libraries.append(adolcLib)

# Extensions
listExtensions = [
    Extension(
        "KCore.kcore",
        sources=["KCore/kcore.cpp"],
        include_dirs=["KCore"] + additionalIncludePaths + [numpyIncDir],
        library_dirs=additionalLibPaths + libraryDirs,
        libraries=libraries + additionalLibs,
        extra_compile_args=Dist.getCppArgs(),
        extra_link_args=Dist.getLinkArgs(),
    )
]

setup(
    name="KCore",
    version="4.1",
    description="Core for *Cassiopee* modules.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['KCore'],
    ext_modules=listExtensions
)
