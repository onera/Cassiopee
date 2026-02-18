#=============================================================================
# Transform requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore library
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
Dist = loadModuleFromPath('../KCore/Dist')
installBase = loadModuleFromPath('../KCore/installBase')
Dist.setConfigDict(installBase.installDict)
additionalLibPaths = Dist.getAdditionalLibPaths()
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLib) = Dist.checkModuleCassiopee("KCore")

prod = os.getenv("ELSAPROD") or "xx"

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLib]
includeDirs = [numpyIncDir, kcoreIncDir]
libraries = ["transform", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# Extensions ==================================================================
listExtensions = []
listExtensions.append(
    Extension('Transform.transform',
              sources=['Transform/transform.cpp'],
              include_dirs=["Transform"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Transform",
    version="4.1",
    description="Transformations of arrays/pyTrees for *Cassiopee* modules.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    package_dir={"":"."},
    packages=['Transform'],
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
installPath = loadModuleFromPath('../KCore/installPath')
installPathDict = {
    "installPath": installPath.installPath,
    "libPath": installPath.libPath,
    "includePath": installPath.includePath
}
Dist.checkPythonPath(installPathDict)
Dist.checkLdLibraryPath(installPathDict)
