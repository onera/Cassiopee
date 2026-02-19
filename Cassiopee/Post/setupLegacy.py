#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os, sys

#=============================================================================
# Post requires:
# C++ compiler
# Fortran compiler
# Numpy
# KCore
#=============================================================================

import KCore.Dist as Dist

# Compiler settings must be set in installBase.py / installBaseUser.py
f77compiler = Dist.getf77Compiler()
f90compiler = Dist.getFromConfigDict("f90compiler", "gfortran")
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibPaths = Dist.getAdditionalLibPaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Compilation des fortrans ===================================================
if f77compiler is None:
    print("Error: a fortran 77 compiler is required for compiling Post.")
    sys.exit()
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT'+str(c)+'='+v+' '
os.system("make -e FC="+f77compiler+" F90=true WDIR=Post/Fortran "+opt)
os.system("make -e FC="+f77compiler+" F90=true WDIR=Post/zipper "+opt)
if f90compiler is not None and os.access('Post/usurp', os.F_OK):
    os.system("(cd Post/usurp; make -e FC="+f77compiler+" F90="+f90compiler+" "+opt+")")
prod = os.getenv("ELSAPROD") or 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["PostF", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

if f90compiler is not None and os.access('Post/usurp', os.F_OK): libraries.append("UsurpF")

import srcs

# extensions =================================================================
listExtensions = []
listExtensions.append(
    Extension('Post.post',
              sources=["Post/post.cpp"]+srcs.cpp_srcs,
              include_dirs=["Post"]+additionalIncludePaths+[numpyIncDir,kcoreIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Post",
    version="4.1",
    description="Post-processing of CFD solutions.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    package_dir={"":"."},
    packages=['Post'],
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
