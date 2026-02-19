#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#=============================================================================
# Initiator requires:
# [ENV] CASSIOPEE, ELSAPROD
# C++ compiler
# Fortran compiler : defined in config.py
# Numpy
# KCore library
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

# Compilation des fortrans ===================================================
if f77compiler is None:
    print("Error: a fortran 77 compiler is required for compiling Initiator.")
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT'+str(c)+'='+v+' '
os.system("make -e FC="+f77compiler+" WDIR=Initiator/Fortran "+opt)
prod = os.getenv("ELSAPROD") or 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["InitiatorF", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

import srcs

# setup =======================================================================
setup(
    name="Initiator",
    version="4.1",
    description="Initiator for *Cassiopee* modules.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    package_dir={"":"."},
    packages=['Initiator'],
    ext_modules=[Extension('Initiator.initiator',
                           sources=['Initiator/initiator.cpp']+srcs.cpp_srcs,
                           include_dirs=["Initiator"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )
                 ]
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
