from setuptools import setup, Extension
import os

#=============================================================================
# KCore requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler
# Numpy
#=============================================================================

import Dist

# Compiler settings must be set in installBase.py / installBaseUser.py
f77compiler = Dist.getf77Compiler()
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibPaths = Dist.getAdditionalLibPaths()
additionalLibs = Dist.getAdditionalLibs()

# Write KCore installation path to installPath.py
Dist.writeInstallPath()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Fortran compilation ========================================================
if f77compiler is None:
    print("Error: a fortran 77 compiler is required for compiling KCore.")
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT%d=%s '%(c, v)
os.system("make -e FC="+f77compiler+" WDIR=KCore/Fld "+opt)
os.system("make -e FC="+f77compiler+" WDIR=KCore/Interp "+opt)
os.system("make -e FC="+f77compiler+" WDIR=KCore/Metric "+opt)
os.system("make -e FC="+f77compiler+" WDIR=KCore/CompGeom "+opt)
os.system("make -e FC="+f77compiler+" WDIR=KCore/Loc "+opt)
os.system("make -e FC="+f77compiler+" WDIR=KCore/Linear "+opt)
prod = os.getenv("ELSAPROD") or 'xx'

# Setting libraries path =====================================================
libraryDirs = ["build/"+prod]
libraries = ["Fld", "Interp", "Metric", "CompGeom", "Loc", "Linear"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# Extensions =================================================================
import srcs
extensions = [
    Extension('KCore.kcore',
              sources=['KCore/kcore.cpp']+srcs.cpp_srcs,
              include_dirs=["KCore", "KCore/Metis"]+additionalIncludePaths+[numpyIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              )
]

# Setup ======================================================================
setup(
    name="KCore",
    version="4.1",
    description="Core for *Cassiopee* modules.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    package_dir={"":"."},
    packages=['KCore'],
    ext_modules=extensions
)

# Check PYTHONPATH
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
