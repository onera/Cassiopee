#!/usr/bin/env python

#=============================================================================
# In installBase.py, you can specify: 
# - the name of the fortran 90 compiler
# fortran 90 compiler is required.
# If you want to specify another C compiler than the one python is
# compiled with, set Cppcompiler. Otherwise, let it to "None".
# When using openMP (useOMP=True), take care of using compatibles Cpp
# and fortran.
# Optional:
# - additional include paths
# - additional library paths
# - additional libraries
# If you want use CPlot as an offscreen plotter on a cluster without a gfx
# card (requires Mesa installed), set CPlotOffScreen to True
#=============================================================================

import platform, re, os

a = platform.uname()
system = a[0] # Linux, Windows
host = a[1]   # host name
prod = os.getenv("ELSAPROD")

#==============================================================================
# Nouveau systeme de configuration par la base d'installation
#==============================================================================
try: import installBase
except: import KCore.installBase as installBase
dict = installBase.installDict
key = ''
# prod est tout d'abord cherche dans le dictionnaire
if prod is not None:
    for i in dict:
        if re.compile(i).search(prod) is not None:
            key = i; break
# puis le uname
if key == '':
    for i in dict:
        if re.compile(i).search(host) is not None:
            key = i; break

if key == '': # not found in install base
    print("Warning: %s was not found in KCore/installBase.py."%host)
    print("Warning: using default compilers and options.")
    print("Warning: to change that, add a block in KCore/installBase.py.")
    key = 'default'

v = dict[key]
#print('%s was found in install base.'%host)
description = v[0]
f77compiler = v[1]
f90compiler = v[2]
Cppcompiler = v[3]
CppAdditionalOptions = v[4]
f77AdditionalOptions = v[5]
useOMP = v[6]
useStatic = v[7]
CPlotOffScreen = v[8]
additionalIncludePaths = v[9]
additionalLibs = v[10]
additionalLibPaths = v[11]
