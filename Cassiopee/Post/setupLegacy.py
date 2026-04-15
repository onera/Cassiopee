#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os, sys
import subprocess

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
def runMakeFortran(f77compiler, opt, f90compiler=None):
    targets = ["Post/Fortran", "Post/zipper"]
    for wdir in targets:
        subprocess.run(
            ["make", "-e", f"FC={f77compiler}", "F90=true", f"WDIR={wdir}", *opt.split()],
            check=True
        )
    usurp_dir = "Post/usurp"
    if f90compiler is not None and os.access(usurp_dir, os.F_OK):
        subprocess.run(
            ["make", "-e", f"FC={f77compiler}", f"F90={f90compiler}", *opt.split()],
            cwd=usurp_dir,
            check=True
        )

if f77compiler is None:
    print("Error: a fortran 77 compiler is required for compiling Post.")
    sys.exit()
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT'+str(c)+'='+v+' '
runMakeFortran(f77compiler, opt, f90compiler=f90compiler)
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
    version="4.2",
    description="Post-processing of CFD solutions.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    package_dir={"":"."},
    packages=['Post'],
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
