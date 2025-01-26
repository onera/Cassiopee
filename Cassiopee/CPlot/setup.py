import os
from distutils.core import setup, Extension
#from setuptools import setup, Extension
import KCore.config

#=============================================================================
# CPlot requires:
# C++ compiler
# Numpy
# KCore
# GL
# optional: MPEG, OSMesa
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

from KCore.config import *
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'
libraryDirs = ["build/"+prod]
includeDirs = [numpyIncDir, kcoreIncDir, kcoreLibDir+'/../include']

libraries = []
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Test if MPEG exists =========================================================
from srcs import MPEG
if MPEG:
    (mpeg, mpegIncDir, mpegLib) = Dist.checkMpeg(additionalLibPaths,
                                                 additionalIncludePaths)
    if mpeg:
        libraries += ["avcodec", "avutil"]
        libraryDirs += [mpegLib]
        includeDirs += [mpegIncDir]

libraryDirs += [kcoreLibDir]
libraries += ["kcore"]

# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi(additionalLibPaths,
                                                     additionalIncludePaths)


# Test if GL exists ==========================================================
#(isGL, GLIncDir, GLLibDir) = Dist.checkGL(additionalLibPaths,
#                                          additionalIncludePaths)
isGL = True; GLIncDir = []; GLLibDir = []

mySystem = Dist.getSystem()
if mySystem[0] == 'mingw' and mySystem[1] == '32':
    libraries += ["wsock32", "winmm", "gdi32"]
    libGL = ['opengl32', 'glu32']
elif mySystem[0] == 'mingw' and mySystem[1] == '64':
    libraries += ["wsock32", "winmm", "gdi32"]
    libGL = ['opengl32', 'glu32']
elif mySystem[0] == 'Darwin':
    libraries += ["X11", "Xmu"]
    libGL = ['GL', 'GLU']
else:
    libraries += ["Xi", "Xmu", "rt"]
    libGL = ['GL', 'GLU']

# Test if OSMesa exists =======================================================
(OSMesa, OSMesaIncDir, OSMesaLibDir, OSMesaLibname) = Dist.checkOSMesa(additionalLibPaths,
                                                                       additionalIncludePaths)

# Extensions =================================================================
EXTRA = ['-D__SHADERS__']
if mpi:
    libraryDirs.append(mpiLibDir)
    includeDirs.append(mpiIncDir)
    libraries += mpiLibs
    EXTRA += ['-D_MPI']
EXTRA += Dist.getCppArgs()

extensions = []
if isGL:
    extensions += [
        Extension('CPlot.cplot',
                  sources=['CPlot/cplot.cpp'],
                  include_dirs=["CPlot", "CPlot/Shaders"]+additionalIncludePaths+includeDirs,
                  library_dirs=additionalLibPaths+libraryDirs,
                  libraries=['cplot', 'cplot1', 'cplot', 'cplot1']+libGL+libraries+additionalLibs,
                  extra_compile_args=EXTRA,
                  extra_link_args=Dist.getLinkArgs())
    ]

if OSMesa:
    extensions += [
        Extension('CPlot.cplotOSMesa',
                  sources=['CPlot/cplotOSMesa.cpp'],
                  include_dirs=["CPlot", "CPlot/Shaders2.0"]+additionalIncludePaths+includeDirs+[OSMesaIncDir],
                  library_dirs=additionalLibPaths+libraryDirs+[OSMesaLibDir],
                  libraries=['cplot', 'cplot2', 'cplot', 'cplot2', OSMesaLibname]+libraries+additionalLibs+libGL,
                  extra_compile_args=EXTRA+['-D__MESA__'],
                  extra_link_args=Dist.getLinkArgs())
    ]


# Setup ======================================================================
setup(
    name="CPlot",
    version="4.0",
    description="A plotter for *Cassiopee* Modules.",
    author="ONERA",
    url="https://cassiopee.onera.fr",
    packages=['CPlot'],
    package_dir={"":"."},
    ext_modules=extensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
