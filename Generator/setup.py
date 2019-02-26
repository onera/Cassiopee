#!/usr/bin/env python
from distutils.core import setup, Extension
import os

#=============================================================================
# Generator requires :
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler
# Numpy
# KCore
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Compilation des fortrans ====================================================
from KCore.config import *
if f77compiler == "None":
    print("Error: a fortran 77 compiler is required for compiling Generator.")
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT'+str(c)+'='+v+' '
os.system("make -e FC="+f77compiler+" WDIR=Generator/Fortran "+opt)
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["GeneratorF", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Hack: compilation of predicate
#os.system(Cppcompiler+" -c -O0 -fPIC Generator/Tetgen/predicates.cxx -o build/temp.linux-x86_64-2.9/Generator/Tetgen/predicates.o")

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
    version="2.9",
    description="*Cassiopee* module of mesh generation.",
    author="Onera",
    package_dir={"":"."},
    packages=['Generator'],
    ext_modules=listExtensions
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
