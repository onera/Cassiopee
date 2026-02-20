#=============================================================================
# OCC requires:
# ELSAPROD variable defined in environment
# C++ compiler
# KCore library
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

# Test if generator exists ===================================================
(generatorVersion, generatorIncDir, generatorLibDir) = Dist.checkModuleCassiopee("Generator")

# Test if open-cascade is already installed ==================================
(OCCPresent, OCCIncDir, OCCLibDir) = Dist.checkOCC()

if not OCCPresent:
    print("Warning: open cascade not found on your system. OCC not installed.")
    os._exit(0)

# Compilation des fortrans ===================================================
prod = os.getenv("ELSAPROD") or "xx"

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir, generatorLibDir]
includeDirs = [numpyIncDir, kcoreIncDir, generatorIncDir]
libraries = ["occ_cassiopee", "generator", "kcore"]

if OCCPresent:
    libraryDirs += [OCCLibDir]
    includeDirs += [OCCIncDir]

srcs = loadModuleFromPath('srcs')
libOCC = Dist.getOCCModules()
if OCCPresent and Dist.getSystem()[0] == 'mingw':
    libOCE = [i+".dll" for i in libOCC]
libraries += libOCC + libOCC

(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# setup ======================================================================
setup(
    name="OCC",
    version="4.1",
    description="OpenCascade python module.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['OCC'],
    package_dir={"":"."},
    ext_modules=[Extension('OCC.occ',
                           sources=["OCC/occ.cpp"],
                           include_dirs=["OCC"]+additionalIncludePaths+includeDirs,
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )]
)
