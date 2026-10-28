#=============================================================================
# Compressor requires:
# C++ compiler
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

# Write build info to file
Dist.writeBuildInfo("Compressor/Compressor.py")

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Setting libraryDirs and libraries ===========================================
prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")
libraryDirs = ['build/'+prod, kcoreLibDir]
libraries = ["compressor", "kcore"]
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

srcs = loadModuleFromPath('srcs')

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
    ext_modules=extensions
)
