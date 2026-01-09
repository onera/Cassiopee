#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#=============================================================================
# Transform requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLib) = Dist.checkKCore()

# Compilation des fortrans ===================================================
from KCore.config import *
if f77compiler == "None":
    print("Error: a fortran 77 compiler is required for compiling Transform.")
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT'+str(c)+'='+v+' '
os.system("make -e FC="+f77compiler+" WDIR=Transform/Fortran "+opt)
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLib]
libraries = ["TransformF", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# setup =======================================================================
import srcs
setup(
    name="Transform",
    version="4.1",
    description="Transformations of arrays/pyTrees for *Cassiopee* modules.",
    author="ONERA",
    package_dir={"":"."},
    packages=['Transform'],
    ext_modules=[Extension('Transform.transform',
                           sources=["Transform/transform.cpp"]+srcs.cpp_srcs,
                           include_dirs=["Transform"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )]
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
