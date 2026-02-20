#=============================================================================
# Modeler requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Numpy
# KCore, Converter, Generator, Transform
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

# Test if open-cascade is already installed ==================================
(OCCPresent, OCCIncDir, OCCLibDir) = Dist.checkOCC()

srcs = loadModuleFromPath('srcs')
if srcs.TIGL and not OCCPresent:
    print("Warning: open cascade not found on your system. Modeler.Tigl not installed.")
    os._exit(0)

prod = os.getenv("ELSAPROD") or "xx"

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir]
includeDirs = [numpyIncDir, kcoreIncDir]
libraries = ["kcore", "modeler", "modeler"]
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

if OCCPresent:
    libraryDirs += [OCCLibDir]
    includeDirs += [OCCIncDir]

libOCC = Dist.getOCCModules()
if OCCPresent and Dist.getSystem()[0] == 'mingw':
    libOCC = [i+".dll" for i in libOCC]
if OCCPresent: libraries += libOCC + libOCC

if srcs.TIXI: libraries += ["curl", "xml2", "xslt"]
if srcs.TIGL: libraries += ["boost_system", "boost_filesystem", "boost_date_time"]
#if srcs.TIGL: libraries += ["boost_filesystem-mt", "boost_date_time-mt"]

# Extensions ==================================================================
listExtensions = []
listExtensions.append(
    Extension('Modeler.modeler',
              sources=['Modeler/modeler.cpp'],
              include_dirs=["Modeler"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Modeler",
    version="4.1",
    description="Modeler module.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Modeler', 'Modeler.CPACS'],
    package_dir={"":"."},
    ext_modules=listExtensions
)
