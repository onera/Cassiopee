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
#=============================================================================

import os
# IN: ELSAPROD only
prod = os.getenv("ELSAPROD")

#==============================================================================
# Nouveau systeme de configuration par la base d'installation
#==============================================================================
try: import KCore.installBase as installBase
except: import installBase

configDict = installBase.installDict
key = ''
# prod est tout d'abord cherche dans le dictionnaire
if prod is not None:
    prod = prod.replace('_i8', '').replace('_DBG', '').split('_b-')[0]
    if prod in configDict: key = prod

if key == '': # not found in install base
    print("Warning: %s were found in KCore/installBase.py."%(prod))
    print("Warning: using default compilers and options.")
    print("Warning: to change that, add a block in KCore/installBase.py.")
    key = 'default'

v = configDict[key]
description = v[0]
f77compiler = v[1]
f90compiler = v[2]
Cppcompiler = v[3]
CppAdditionalOptions = v[4]
f77AdditionalOptions = v[5]
useOMP = v[6]
useStatic = v[7]
additionalIncludePaths = v[8]
additionalLibs = v[9]
additionalLibPaths = v[10]
useCuda = v[11]
NvccAdditionalOptions = v[12]
