#=============================================================================
# Intersector requires:
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
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Test if xcore exists =======================================================
#(xcoreVersion, xcoreIncDir, xcoreLibDir) = Dist.checkModuleCassiopee("XCore")()

# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi()
(mpi4py, mpi4pyIncDir, mpi4pyLibDir) = Dist.checkMpi4py()

# Compilation des fortrans ====================================================
prod = os.getenv("ELSAPROD") or "xx"

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir]
includeDirs = [numpyIncDir, kcoreIncDir]
libraries = ["intersector", "kcore"]
ADDITIONALCPPFLAGS = []
if mpi:
    libraryDirs.append(mpiLibDir)
    includeDirs.append(mpiIncDir)
    ADDITIONALCPPFLAGS += ['-D_MPI']
if mpi4py:
    includeDirs.append(mpi4pyIncDir)

if mpi: libraries += mpiLibs

(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# setup ======================================================================
setup(
    name="Intersector",
    version="4.1",
    description="Mesh-intersection-based services in *Cassiopee*.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Intersector'],
    package_dir={"":"."},
    ext_modules=[Extension('Intersector.intersector',
                           sources=["Intersector/intersector.cpp"],
                           include_dirs=["Intersector"]+additionalIncludePaths+includeDirs,
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs()+ADDITIONALCPPFLAGS,
                           extra_link_args=Dist.getLinkArgs()
                           )]
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
