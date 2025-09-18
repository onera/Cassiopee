#from distutils.core import setup, Extension
from setuptools import setup, Extension
import KCore.config
import os

#=============================================================================
# CPlot requires:
# C++ compiler
# Numpy
# KCore
# GL
# optional: PNG, MPEG, OSMesa
#=============================================================================

# If you want to use CPlot as a offscreen plotter (as on clusters)
# set UseOSMesa to True (requires mesa)
UseOSMesa = KCore.config.CPlotOffScreen

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

libraries = ["GLU", "kcore", "Xi", "Xmu", "rt"]
from KCore.config import *
if not UseOSMesa: libraries += ["GL"]

#libraryDirs = [kcoreLibDir]
libraryDirs = [kcoreLibDir, '/usr/X11R6/lib64']
includeDirs = [numpyIncDir, kcoreIncDir]

(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Test if PNG exists =========================================================
(png, pngIncDir, pngLib) = Dist.checkPng(additionalLibPaths,
                                         additionalIncludePaths)
if png:
    libraries += ["png"]
    libraryDirs += [pngLib]
    includeDirs += [pngIncDir]

# Test if MPEG exists =========================================================
(mpeg, mpegIncDir, mpegLib) = Dist.checkMpeg(additionalLibPaths,
                                             additionalIncludePaths)
if mpeg:
    libraries += ["avcodec"]
    libraryDirs += [mpegLib]
    includeDirs += [mpegIncDir]

# Test if OSMesa exists =======================================================
# Put this to True for using CPlot in batch mode
if UseOSMesa:
    (OSMesa, OSMesaIncDir, OSMesaLib) = Dist.checkOSMesa(additionalLibPaths,
                                                         additionalIncludePaths)
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
    author="Onera",
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
