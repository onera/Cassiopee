#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#=============================================================================
# OCC requires:
# ELSAPROD variable defined in environment
# C++ compiler
# KCore library
# Generator library
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Test if generator exists ===================================================
(generatorVersion, generatorIncDir, generatorLibDir) = Dist.checkGenerator()

# Setting libraryDirs and libraries ==========================================
libraryDirs = ["build/"+prod, kcoreLibDir, generatorLibDir]
includeDirs = [kcoreIncDir, generatorIncDir]
libraries = ["generator", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# setup ======================================================================
import srcs
setup(
    name="OCC",
    version="2.7",
    description="OpenCascade link.",
    author="Onera",
    package_dir={"":"."},
    packages=['OCC'],
    ext_modules=[Extension('OCC.occ',
                           sources=["OCC/occ.cpp"]+srcs.cpp_srcs,
                           include_dirs=[["OCC", "OCC/occ_inc"]]+additionalIncludePaths+[kcoreIncDir, generatorIncDir],
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )]
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
