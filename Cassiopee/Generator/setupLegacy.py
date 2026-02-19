#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#=============================================================================
# Generator requires :
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler
# Numpy
# KCore
#=============================================================================

import KCore.Dist as Dist

# Compiler settings must be set in installBase.py / installBaseUser.py
f77compiler = Dist.getf77Compiler()
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibPaths = Dist.getAdditionalLibPaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Compilation des fortrans ====================================================
if f77compiler is None:
    print("Error: a fortran 77 compiler is required for compiling Generator.")
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT'+str(c)+'='+v+' '
os.system("make -e FC="+f77compiler+" WDIR=Generator/Fortran "+opt)
prod = os.getenv("ELSAPROD") or 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["GeneratorF", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# Extensions =================================================================
import srcs
listExtensions = []
listExtensions.append(
    Extension('Generator.generator',
              sources=["Generator/generator.cpp"]+srcs.cpp_srcs+srcs.cpp_srcs2,
              include_dirs=["Generator", "Generator/Netgen/include"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup =======================================================================
setup(
    name="Generator",
    version="4.1",
    description="*Cassiopee* module of mesh generation.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    package_dir={"":"."},
    packages=['Generator'],
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
