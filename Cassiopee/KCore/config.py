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

import platform, re, os

a = platform.uname()
system = a[0] # Linux, Windows
host = a[1]   # host name
prod = os.getenv("ELSAPROD")

#==============================================================================
# Nouveau systeme de configuration par la base d'installation
#==============================================================================
try: import KCore.installBase as installBase
except: import installBase

dict = installBase.installDict
key = ''
# prod est tout d'abord cherche dans le dictionnaire
if prod is not None:
    if prod.endswith('_DBG'): prod = prod[:-4]
    if prod.endswith('_i8'): prod = prod[:-3]
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
additionalIncludePaths = v[8]
additionalLibs = v[9]
additionalLibPaths = v[10]
useCuda = v[11]
NvccAdditionalOptions = v[12]

# this part to add intel compilers to distutils
if Cppcompiler.find('icc') == 0 or Cppcompiler.find('icpc') == 0:
    def new_compiler(plat=None, compiler=None, verbose=0, dry_run=0, force=0):
        from distutils.unixccompiler import UnixCCompiler
        compiler = UnixCCompiler(None, dry_run, force)
        compiler.cc_exe = Cppcompiler
        compiler.set_executables(compiler=Cppcompiler, compiler_cxx=Cppcompiler, compiler_so=Cppcompiler,
                                 linker_exe=Cppcompiler, linker_so=Cppcompiler+' -shared')
        return compiler

    from distutils import ccompiler
    ccompiler.new_compiler = new_compiler
