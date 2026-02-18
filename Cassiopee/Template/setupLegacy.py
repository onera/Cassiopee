#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#=============================================================================
# Template requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore
#=============================================================================

import KCore.Dist as Dist

# Compiler settings must be set in installBase.py / installBaseUser.py
# f77compiler = Dist.getf77Compiler()
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
#if f77compiler is None:
#    print("Error: a fortran 77 compiler is required for compiling Fast.")
#args = Dist.getForArgs(); opt = ''
#for c in range(len(args)):
#    opt += 'FOPT'+str(c)+'='+args[c]+' '
#os.system("make -e FC="+f77compiler+" WDIR=Template/Fortran "+opt)
prod = os.getenv("ELSAPROD") or 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir]
includeDirs = [numpyIncDir, kcoreIncDir]
libraries = ["kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# Extensions ==================================================================
import srcs
listExtensions = []
listExtensions.append(
    Extension('Template.template',
              sources=['Template/template.cpp']+srcs.cpp_srcs,
              include_dirs=["Template"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Template",
    version="4.1",
    description="Template module.",
    author="You",
    package_dir={"":"."},
    packages=['Template'],
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
