#=============================================================================
# Dist2Walls requires:
# C++ compiler
# Numpy
# KCore
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
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Setting libraryDirs and libraries ===========================================
prod = os.getenv("ELSAPROD") or "xx"
libraryDirs = [kcoreLibDir, 'build/'+prod]
libraries = ["dist2walls", "kcore"]

(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs

# Extensions =================================================================
extensions = [
    Extension('Dist2Walls.dist2walls',
              sources=["Dist2Walls/dist2walls.cpp"],
              include_dirs=["Dist2Walls"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              )
]

# Setup ======================================================================
setup(
    name="Dist2Walls",
    version="4.1",
    description="Computation of distance to walls.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Dist2Walls'],
    package_dir={"":"."},
    ext_modules=extensions
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
