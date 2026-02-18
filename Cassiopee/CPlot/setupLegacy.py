#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#=============================================================================
# CPlot requires:
# C++ compiler
# Numpy
# KCore
# GL
# optional: PNG, MPEG, OSMesa
#=============================================================================

import KCore.Dist as Dist
UseOSMesa = True

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

libraries = ["GLU", "kcore", "Xi", "Xmu", "rt"]
if not UseOSMesa: libraries += ["GL"]

#libraryDirs = [kcoreLibDir]
libraryDirs = [kcoreLibDir, '/usr/X11R6/lib64']
includeDirs = [numpyIncDir, kcoreIncDir]

(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

# Test if PNG exists =========================================================
(png, pngIncDir, pngLib) = Dist.checkPng()
if png:
    libraries += ["png"]
    libraryDirs += [pngLib]
    includeDirs += [pngIncDir]

# Test if MPEG exists =========================================================
(mpeg, mpegIncDir, mpegLib) = Dist.checkMpeg()
if mpeg:
    libraries += ["avcodec"]
    libraryDirs += [mpegLib]
    includeDirs += [mpegIncDir]

# Test if OSMesa exists =======================================================
# Put this to True for using CPlot in batch mode
if UseOSMesa:
    (OSMesa, OSMesaIncDir, OSMesaLib) = Dist.checkOSMesa()
    if OSMesa:
        libraries += ["OSMesa"]
        libraryDirs += [OSMesaLib]
        includeDirs += [OSMesaIncDir]
else: OSMesa = False

# Extensions =================================================================
import srcs

import KCore.installPath
EXTRA = ['-D__SHADERS__']
if OSMesa: EXTRA += ['-D__MESA__']

cppargs = Dist.getCppArgs()
EXTRA += cppargs

extensions = [
    Extension('CPlot.cplot',
              sources=['CPlot/cplot.cpp']+srcs.cpp_srcs,
              include_dirs=["CPlot"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=EXTRA,
              extra_link_args=Dist.getLinkArgs()
              )
]

# Setup ======================================================================
setup(
    name="CPlot",
    version="4.1",
    description="A plotter for *Cassiopee* Modules.",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    package_dir={"":"."},
    packages=['CPlot'],
    ext_modules=extensions
)

# Install shaders + textures ==================================================
os.system("cp CPlot/Shaders/*.vert %s/CPlot/"%KCore.installPath.installPath)
os.system("cp CPlot/Shaders/*.frag %s/CPlot/"%KCore.installPath.installPath)
os.system("cp CPlot/Shaders/*.geom %s/CPlot/"%KCore.installPath.installPath)
os.system("cp CPlot/Shaders/*.tcs %s/CPlot/"%KCore.installPath.installPath)
os.system("cp CPlot/Shaders/*.tes %s/CPlot/"%KCore.installPath.installPath)
os.system("cp CPlot/Textures/*.png %s/CPlot/"%KCore.installPath.installPath)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
