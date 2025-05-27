# Functions used in *Cassiopee* modules setup.py
import os, sys, platform, glob, subprocess
from distutils import sysconfig
from distutils.core import Extension
#from setuptools import sysconfig
#from setuptools import Extension

# Toggle to True for compiling for debug (valgrind, inspector, sanitizer)
DEBUG = False

# Toggle to True for compiling all integers in i8
EDOUBLEINT = False
# Toggle to True for compiling global index in i8
GDOUBLEINT = False

#==============================================================================
# Check module import
# Write SUCCESS or FAILED (with colored output)
#==============================================================================
def checkModuleImport(moduleName):
    # Remove . from PYTHONPATH
    try: del sys.path[sys.path.index('')]
    except: pass
    try: del sys.path[sys.path.index(os.getcwd())]
    except: pass
    # Try to detect if colored output is supported
    color = False
    if sys.stdout.isatty(): color = True

    # sec / lock
    #moduleBase = moduleName.split('.')[0]
    #os.chmod(moduleBase, 0o700)

    # try import module
    try:
        __import__(moduleName)
        if color:
            print("\033[32m%s correctly installed.\033[0m"%moduleName)
        else: print("%s correctly installed."%moduleName)
    except Exception as inst:
        if color:
            print("\033[31mFAILED: %s\033[0m"%inst)
            print("\033[31mFAILED: %s badly installed.\033[0m"%moduleName)
        else:
            print("FAILED: %s"%inst)
            print("FAILED: %s badly installed."%moduleName)

#==============================================================================
# Return informations on the current operating system
# Return: Unix, Windows, Darwin, Java, mingw + bits of the system ('32', '64')
#==============================================================================
def getSystem():
    # Detection Mingw
    try:
        # On se base en premier sur ELSAPROD
        key = os.environ['ELSAPROD']
        if key == 'win64': return ['mingw', '64']
        elif key == 'win32': return ['mingw', '32']
        elif key == 'msys64': return ['mingw', '64']
        elif key == 'msys64p3': return ['mingw', '64']
        elif key == 'msys32': return ['mingw', '32']
        key = os.environ['MSYSTEM']
        if key == 'MINGW32': return ['mingw', '32']
        elif key == 'MINGW64': return ['mingw', '64']
    except: pass
    # System bits
    bits = '32'
    try:
        bits = platform.architecture()[0]
        bits = bits[0:2]
    except: pass
    # Windows, unix, mac
    name = platform.uname()
    return [name[0], bits]

#==============================================================================
# Get name in environ
# Return '' if name is not in environ
#==============================================================================
def getenv(name):
    if name in os.environ: return os.environ[name]
    else: return ''

#==============================================================================
# Get name of the Data folder
#==============================================================================
def getDataFolderName(name='Data'):
    elsaprod = os.getenv("ELSAPROD")
    if elsaprod is not None:
        if not '_i8' in elsaprod and EDOUBLEINT:
            print("Warning: ELSAPROD {} compiled in i8 but recommended suffix "
                  "'_i8' is missing".format(elsaprod))
        if not '_DBG' in elsaprod and DEBUG:
            print("Warning: ELSAPROD {} compiled in DEBUG but recommended "
                  "suffix '_DBG' is missing".format(elsaprod))
        name += '_' + elsaprod
    else:
        name += '_xx'
    if sys.version_info[0] == 2: name += '2'
    return name

#==============================================================================
# Check all (python, numpy, C++, fortran, hdf, mpi, mpi4py, png, osmesa, mpeg)
#==============================================================================
def checkAll(summary=True):
    from config import additionalLibPaths, additionalIncludePaths, useCuda

    out = []
    try:
        (pythonVersion, pythonIncDir, pythonLibDir, pythonLibs) = checkPython()
        out += ['Python: %s'%pythonVersion]
    except: out += ['Python: include is missing.']
    try:
        (numpyVersion, numpyIncDir, numpyLibDir) = checkNumpy()
        out += ['Numpy: %s'%numpyVersion]
    except: out += ['Numpy: is missing or numpy includes are missing.']
    (ok, CppLibs, CppLibPaths) = checkCppLibs(additionalLibPaths, additionalIncludePaths)
    if ok: out += ['C++: OK (%s, %s).'%(CppLibs, CppLibPaths)]
    else: out += ['C++: Fail.']
    (ok, FLibs, FLibPaths) = checkFortranLibs(additionalLibPaths, additionalIncludePaths)
    if ok: out += ['f77: OK (%s, %s).'%(FLibs, FLibPaths)]
    else: out += ['f77: Fail.']
    (ok, hdfIncDir, hdfLibDir, hdflibs) = checkHdf(additionalLibPaths, additionalIncludePaths)
    if ok: out += ['hdf: OK (%s, %s).'%(hdfIncDir, hdfLibDir)]
    else: out += ['hdf: missing (%s, %s).'%(hdfIncDir, hdfLibDir)]
    (ok, pngIncDir, pngLibDir) = checkPng(additionalLibPaths, additionalIncludePaths)
    if ok: out += ['png: OK (%s, %s).'%(pngIncDir, pngLibDir)]
    else: out += ['png: missing (%s, %s).'%(pngIncDir, pngLibDir)]
    (ok, osmesaIncDir, osmesaLibDir, libname) = checkOSMesa(additionalLibPaths, additionalIncludePaths)
    if ok: out += ['osmesa: OK (%s, %s).'%(osmesaIncDir, osmesaLibDir)]
    else: out += ['osmesa: missing (%s, %s).'%(osmesaIncDir, osmesaLibDir)]
    (ok, mpegIncDir, mpegLibDir) = checkMpeg(additionalLibPaths, additionalIncludePaths)
    if ok: out += ['mpeg: OK (%s, %s).'%(mpegIncDir, mpegLibDir)]
    else: out += ['mpeg: missing (%s, %s).'%(mpegIncDir, mpegLibDir)]
    (ok, mpiIncDir, mpiLibDir, mpiLibs) = checkMpi(additionalLibPaths, additionalIncludePaths)
    if ok: out += ['mpi: OK (%s, %s).'%(mpiIncDir, mpiLibDir)]
    else: out += ['mpi: missing (%s %s).'%(mpiIncDir, mpiLibDir)]
    (ok, mpi4pyIncDir, mpi4pyLibDir) = checkMpi4py(additionalLibPaths, additionalIncludePaths)
    if ok: out += ['mpi4py: OK (%s).'%(mpi4pyIncDir)]
    else: out += ['mpi4py: missing (%s).'%(mpi4pyIncDir)]

    if useCuda:
        (ok, cudaIncDir, cudaLib, cudaBin) = checkCuda(additionalLibPaths, additionalIncludePaths)
        if ok: out += ['cuda: used (%s)'%(cudaIncDir)]
        else: out += ['cuda: missing. Not used (%s).'%(cudaIncDir)]
    if summary:
        print('Summary:')
        print('========')
        for i in out: print(i)

#==============================================================================
# Check python includes / libs
#==============================================================================
def checkPython():
    pythonVersion = sysconfig.get_python_version()
    #vars = sysconfig.get_config_vars()

    pythonIncDir = sysconfig.get_python_inc()
    if not os.path.exists(pythonIncDir):
        raise SystemError("Error: Python includes are required for the compilation of Cassiopee modules.")

    pythonLibDir = sysconfig.get_python_lib()
    try:
        a = sysconfig.get_config_var('LDLIBRARY')
        a = a.replace('lib', '')
        a = a.replace('.a', '')
        a = a.replace('.so', '')
        pythonLibs = [a]
    except: pythonLibs = []
    return (pythonVersion, pythonIncDir, pythonLibDir, pythonLibs)

#=============================================================================
# Check numpy include / module
#=============================================================================
def checkNumpy():
    numpyVersion = False
    numpyIncDir = ''
    try:
        import numpy
        numpyIncDir = numpy.get_include()
        numpyVersion = numpy.__version__
    except ImportError:
        raise SystemError("Error: numpy is required for the compilation of Cassiopee modules.")

    if not os.path.exists(numpyIncDir):
        raise SystemError("Error: numpy includes are required for the compilation of Cassiopee modules.")

    return (numpyVersion, numpyIncDir, '')

#=============================================================================
# Retourne le chemin d'installation des modules comme cree par distUtils
# if type=0, return full install path
# else return dict {'lib': 'lib', 'pyversion': 'python3.12', 'site': 'site-packages'}
#=============================================================================
def getInstallPath(prefix, type=0):
    mySystem = getSystem()[0]; bits = getSystem()[1]

    import site
    a = site.getsitepackages()[0].split('/')[-4:]
    if type == 0:
        if a[0] != 'local':
            installPath = '%s/%s/%s/%s'%(prefix, a[1], a[2], a[3])  # 'prefix/lib/python3.12/site-packages'
        else:
            installPath = '%s/%s/%s/%s/%s'%(prefix, a[0], a[1], a[2], a[3])  # 'prefix/local/lib/python3.12/site-packages'
    else:
        installPath = {'lib': a[1], 'pyversion': a[2], 'site': a[3]}  # {'lib': 'lib', 'pyversion': 'python3.12', 'site': 'site-packages'}
    return installPath

    '''
    # Based on distutils (to be sure)
    if os.environ['ELSAPROD'][0:6] == 'msys64' or os.environ['ELSAPROD'] == 'win64':
        pythonLib = sysconfig.get_python_lib()
        pythonLib = pythonLib.split('/')
        pythonVersion = pythonLib[-2]
        Site = pythonLib[-1]
        Lib = pythonLib[-3]
        installPath = '%s/%s/%s/site-packages'%(prefix, Lib, pythonVersion)
    elif os.environ['ELSAPROD'][0:6] == 'ubuntu': # debian style
        pythonLib = sysconfig.get_python_lib()
        pythonLib = pythonLib.split('/')
        pversion = sys.version_info
        pythonVersion = "python{}.{}".format(pversion[0], pversion[1])
        Site = pythonLib[-1]
        Lib = pythonLib[-3]
        installPath = '%s/local/%s/%s/dist-packages'%(prefix, Lib, pythonVersion)
    elif mySystem == 'Windows' or mySystem == 'mingw':
        installPath = prefix + "/Lib/site-packages"
    elif mySystem == 'Darwin':
        pythonLib = sysconfig.get_python_lib()
        pythonLib = pythonLib.split('/')
        pythonVersion = pythonLib[-2]
        installPath = prefix + '/lib/python'+pythonVersion+'/site-packages'
    else: # standard unix
        pythonLib = sysconfig.get_python_lib()
        pythonLib = pythonLib.split('/')
        # Based on python lib
        #installPath = prefix + '/' + '/'.join(pythonLib[-3:])
        # Python version
        pversion = sys.version_info
        pythonVersion = "python{}.{}".format(pversion[0], pversion[1])
        Site = pythonLib[-1]
        Lib = pythonLib[-3]
        installPath = '%s/%s/%s/site-packages'%(prefix, Lib, pythonVersion)

    # temporary for tests
    if installPath != retn:
        print("WARNING: new installPath is not correct.")
        print("WARNING: old: ", installPath)
        print("WARNING: new: ", retn)

    if type == 0: return installPath
    else: return {'lib': Lib, 'pyversion': pythonVersion, 'site': Site}
    '''

#==============================================================================
# Functions returning the names of the remote repo & branch and the commit hash
#==============================================================================
def getGitOrigin(cassiopeeIncDir):
    mySystem = getSystem()[0]
    if mySystem == 'mingw' or mySystem == 'Windows':
        lpath = cassiopeeIncDir.replace('/', '\\')
        cmd = "cd {} && git config --get remote.origin.url".format(lpath)
    else: # unix
        lpath = cassiopeeIncDir
        cmd = "cd {}; git config --get remote.origin.url 2>/dev/null".format(lpath)
    try:
        origin = subprocess.check_output(cmd, shell=True)
        return origin.decode('utf-8', 'ignore').strip()
    except: return "unknown"

def getGitBranch(cassiopeeIncDir):
    mySystem = getSystem()[0]
    if mySystem == 'mingw' or mySystem == 'Windows':
        lpath = cassiopeeIncDir.replace('/', '\\')
        cmd = "cd {} && git rev-parse --abbrev-ref HEAD".format(lpath)
    else: # unix
        lpath = cassiopeeIncDir
        cmd = "cd {}; git rev-parse --abbrev-ref HEAD 2>/dev/null".format(lpath)
    try:
        branchName = subprocess.check_output(cmd, shell=True)
        return branchName.decode('utf-8', 'ignore').strip()
    except: return "unknown"

def getGitHash(cassiopeeIncDir):
    mySystem = getSystem()[0]
    if mySystem == 'mingw' or mySystem == 'Windows':
        lpath = cassiopeeIncDir.replace('/', '\\')
        cmd = "cd {} && git rev-parse --short HEAD".format(lpath)
    else: # unix
        lpath = cassiopeeIncDir
        cmd = "cd {}; git rev-parse --short HEAD 2>/dev/null".format(lpath)
    try:
        sha = subprocess.check_output(cmd, shell=True)
        return sha.decode('utf-8', 'ignore').strip()
    except: return "unknown"

#=============================================================================
# nom de l'egg cree par setup tools
#=============================================================================
def getEggModuleDirName(moduleName, version):
    import platform
    pversion = sys.version_info
    #mod = __import__(moduleName)
    fullName = moduleName+'-'+version+'-py%d.%d'%(pversion[0],pversion[1])+'-'+platform.system().lower()+'-'+platform.processor()+'.egg'
    return fullName

#=============================================================================
# Write installPath, the installation path of Cassiopee to installPath.py
#=============================================================================
def writeInstallPath():
    import re
    prefix = sys.prefix
    a = sys.argv
    for c, i in enumerate(a):
        if re.compile('--prefix=').search(i) is not None: prefix = i[9:] # setup
        elif re.compile('prefix=').search(i) is not None: prefix = i[7:] # setup
        elif re.compile('--prefix').search(i) is not None: prefix = a[c+1] # pip

    installPath = getInstallPath(prefix)

    p = open('installPath.py', 'w')
    if p is None:
        raise SystemError("Error: can not open file installPath.py for writing.")
    p.write('installPath = \'%s\'\n'%installPath)

    import site
    a = site.getsitepackages()[0].split('/')[-4:]
    if a[0] != 'local': libPath = '%s/%s'%(prefix, a[1])  # 'prefix/lib'
    else: libPath = '%s/%s/%s'%(prefix, a[0], a[1])  # 'prefix/local/lib'
    p.write('libPath = \'%s\'\n'%libPath)

    '''
    mySystem = getSystem()[0]; bits = getSystem()[1]
    if mySystem == 'Windows' or mySystem == 'mingw': Lib = 'Lib'
    elif mySystem == 'Darwin': Lib = 'lib'
    else:
        pythonLib = sysconfig.get_python_lib()
        pythonLib = pythonLib.split('/')
        Lib = pythonLib[-3]
    if os.environ['ELSAPROD'][0:6] == 'ubuntu': # debian style
        libPath = '%s/local/%s'%(prefix,Lib)
    else:
        libPath = '%s/%s'%(prefix,Lib)
    p.write('libPath = \'%s\'\n'%libPath)
    '''

    cwd = os.getcwd()
    p.write('includePath = \'%s\'\n'%(cwd))
    gitOrigin = getGitOrigin(cwd)
    gitBranch = getGitBranch(cwd)
    gitHash = getGitHash(cwd)[:7]
    p.write('gitOrigin = \'%s\'\n'%(gitOrigin))
    p.write('gitBranch = \'%s\'\n'%(gitBranch))
    p.write('gitHash = \'%s\'\n'%(gitHash))
    p.close()

#==============================================================================
# Write env files
# Directement dans le repertoire d'installation
#==============================================================================
def writeEnvs():
    try: import KCore.installPath as K
    except: import installPath as K
    libPath = K.libPath
    installPathLocal = K.installPath
    env = os.environ
    cassiopee = env.get('CASSIOPEE', '')
    elsaprod = env.get('ELSAPROD', '')
    if cassiopee != '': envPath = libPath+'/../../../'
    else: envPath = libPath+'/../'
    cmdPath = libPath+'/..'
    installLD = os.getenv('LD_LIBRARY_PATH')

    # max cores
    try:
        import multiprocessing
        mt = multiprocessing.cpu_count()
    except: mt = 1

    # sh ou bash
    # usage: source $CASSIOPEE/Dist/env_Cassiopee.sh
    with open(envPath+"env_Cassiopee.sh", 'w') as p:
        p.write("ulimit -s unlimited\n")
        if cassiopee != '': p.write("export CASSIOPEE=%s\n"%cassiopee)
        if elsaprod != '': p.write("export ELSAPROD=%s\n"%elsaprod)
        p.write("export OMP_NUM_THREADS=%d\n"%mt)
        p.write("export PATH=%s:%s/bin:$PATH\n"%(cmdPath,cmdPath))
        p.write("if [ \"$PYTHONPATH\" = \"\" ]; then\n")
        p.write("      export PYTHONPATH=%s\n"%installPathLocal)
        p.write("else\n")
        p.write("      export PYTHONPATH=%s:$PYTHONPATH\n"%installPathLocal)
        p.write("fi\n")
        if installLD is None:
            p.write("if [ \"$LD_LIBRARY_PATH\" = \"\" ]; then\n")
            p.write("      export LD_LIBRARY_PATH=%s\n"%libPath)
            p.write("else\n")
            p.write("      export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n"%libPath)
            p.write("fi\n")
        else:
            p.write("if [ \"$LD_LIBRARY_PATH\" = \"\" ]; then\n")
            p.write("      export LD_LIBRARY_PATH=%s:%s\n"%(libPath,installLD))
            p.write("else\n")
            p.write("      export LD_LIBRARY_PATH=%s:%s:$LD_LIBRARY_PATH\n"%(libPath,installLD))
            p.write("fi\n")

    # csh ou tcsh
    # usage: source $CASSIOPEE/Dist/env_Cassiopee.csh
    with open(envPath+"env_Cassiopee.csh", 'w') as p:
        p.write("limit stacksize unlimited\n")
        if cassiopee != '': p.write("setenv CASSIOPEE %s\n"%cassiopee)
        if elsaprod != '': p.write("setenv ELSAPROD %s\n"%elsaprod)
        p.write("setenv OMP_NUM_THREADS %d\n"%mt)
        p.write("set path=(%s %s/bin $path)\n"%(cmdPath,cmdPath))
        p.write("if ($?PYTHONPATH == 0) then\n")
        p.write("     setenv PYTHONPATH %s\n"%installPathLocal)
        p.write("else\n")
        p.write("     setenv PYTHONPATH %s:$PYTHONPATH\n"%installPathLocal)
        p.write("endif\n")
        if installLD is None:
            p.write("if ($?LD_LIBRARY_PATH == 0) then\n")
            p.write("     setenv LD_LIBRARY_PATH %s\n"%libPath)
            p.write("else\n")
            p.write("     setenv LD_LIBRARY_PATH %s:$LD_LIBRARY_PATH\n"%libPath)
            p.write("endif\n")
        else:
            p.write("if ($?LD_LIBRARY_PATH == 0) then\n")
            p.write("     setenv LD_LIBRARY_PATH %s:%s\n"%(libPath,installLD))
            p.write("else\n")
            p.write("     setenv LD_LIBRARY_PATH %s:%s:$LD_LIBRARY_PATH\n"%(libPath,installLD))
            p.write("endif\n")

    # bat
    with open(envPath+"env_Cassiopee.bat", 'w') as p:
        p.write("path = "+libPath+";"+cmdPath+"%PATH%\n")
        p.write("set PYTHONPATH="+installPathLocal+";%PYTHONPATH%\n")
        p.write("set OMP_NUM_THREADS=%NUMBER_OF_PROCESSORS%\n")

    # module
    # usage: module use $CASSIOPEE/Dist
    # module load cassiopee
    with open(envPath+"cassiopee", 'w') as p:
        p.write("#%Module1.0#####################################################################\n")
        p.write("##\n")
        p.write("## CASSIOPEE\n")
        p.write("##\n")
        p.write("module-whatis   \"Set the environment for using Cassiopee\"\n")
        if cassiopee != '': p.write("setenv CASSIOPEE %s\n"%cassiopee)
        if elsaprod != '': p.write("setenv ELSAPROD %s\n"%elsaprod)
        p.write("setenv OMP_NUM_THREADS %d\n"%mt)
        p.write("prepend-path PATH %s\n"%cmdPath)
        p.write("prepend-path PATH %s/bin\n"%cmdPath)
        p.write("prepend-path PYTHONPATH %s\n"%installPathLocal)
        if installLD is not None:
            p.write("prepend-path LD_LIBRARY_PATH %s\n"%installLD)
        p.write("prepend-path LD_LIBRARY_PATH %s\n"%libPath)

#==============================================================================
# Write setup.cfg en fonction du compilateur C++ (si different de None)
# setup.cfg est utilise par setup de python pour choisir le compilo.
#==============================================================================
def writeSetupCfg():
    try: from KCore.config import Cppcompiler
    except: from config import Cppcompiler
    mySystem = getSystem()

    # Windows + mingw
    if mySystem[0] == 'mingw' and mySystem[1] == '32':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=mingw32\n')
        p.close(); return
    if mySystem[0] == 'mingw' and mySystem[1] == '64':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=mingw32\n')
        p.close(); return

    # Unix
    if Cppcompiler == "None" or Cppcompiler == "":
        a = os.access("./setup.cfg", os.F_OK)
        if a: os.remove("./setup.cfg")
    elif Cppcompiler == 'icc' or Cppcompiler == 'icpc':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=unix\n')
        p.close()
    elif Cppcompiler == 'icx' or Cppcompiler == 'icpx':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=unix\n')
        p.close()
    elif Cppcompiler == 'gcc' or Cppcompiler == 'g++':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=unix\n')
        p.close()
    elif Cppcompiler == 'pgcc' or Cppcompiler == 'pgc++':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=unix\n')
        p.close()
    elif Cppcompiler == 'nvc' or Cppcompiler == 'nvc++':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=unix\n')
        p.close()
    elif Cppcompiler == 'craycc' or Cppcompiler == 'craycxx':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=unix\n')
        p.close()
    elif Cppcompiler == 'cc' or Cppcompiler == 'cc':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=unix\n')
        p.close()
    elif Cppcompiler == 'clang' or Cppcompiler == 'clang++':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=unix\n')
        p.close()
    else:
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=%s\n'%Cppcompiler)
        p.close()

#==============================================================================
# Retourne le compilo c, c++ et ses options tels que definis dans distutils
# ou dans config.py (installBase.py)
#==============================================================================
def getDistUtilsCompilers():
    vars = sysconfig.get_config_vars('CC', 'CXX', 'OPT',
                                     'BASECFLAGS', 'CCSHARED',
                                     'LDSHARED', 'SO')
    for i, v in enumerate(vars):
        if v is None: vars[i] = ""

    try: from KCore.config import Cppcompiler
    except: from config import Cppcompiler
    if Cppcompiler != 'None' or Cppcompiler != '':
        if Cppcompiler == 'clang':
            vars[0] = Cppcompiler; vars[1] = Cppcompiler.replace('clang', 'clang++')
        elif Cppcompiler == 'clang++':
            vars[0] = Cppcompiler.replace('clang++', 'clang'); vars[1] = Cppcompiler

        elif Cppcompiler == 'pgcc':
            vars[0] = Cppcompiler; vars[1] = Cppcompiler.replace('pgcc', 'pgc++')
        elif Cppcompiler == 'pgc++':
            vars[0] = Cppcompiler.replace('pgc++', 'pgcc'); vars[1] = Cppcompiler

        elif Cppcompiler == 'craycc':
            vars[0] = Cppcompiler; vars[1] = Cppcompiler.replace('craycc', 'craycxx')
        elif Cppcompiler == 'craycxx':
            vars[0] = Cppcompiler.replace('craycxx', 'craycc'); vars[1] = Cppcompiler

        elif Cppcompiler == 'cc':
            vars[0] = Cppcompiler; vars[1] = Cppcompiler

        elif Cppcompiler == 'nvc':
            vars[0] = Cppcompiler; vars[1] = Cppcompiler.replace('nvc', 'nvc++')
        elif Cppcompiler == 'nvc++':
            vars[0] = Cppcompiler.replace('nvc++', 'nvc'); vars[1] = Cppcompiler

        elif Cppcompiler.find('g++') != -1: # g++-version + mingw-g++-version
            vars[0] = Cppcompiler.replace('g++', 'gcc'); vars[1] = Cppcompiler
        elif Cppcompiler.find('gcc') != -1:
            vars[0] = Cppcompiler; vars[1] = Cppcompiler.replace('gcc', 'g++')

        elif Cppcompiler == 'icpc':
            vars[0] = Cppcompiler.replace('icpc', 'icc'); vars[1] = Cppcompiler
        elif Cppcompiler == 'icc':
            vars[0] = Cppcompiler; vars[1] = Cppcompiler.replace('icc', 'icpc')

        elif Cppcompiler == 'icpx':
            vars[0] = Cppcompiler.replace('icpx', 'icx'); vars[1] = Cppcompiler
        elif Cppcompiler == 'icx':
            vars[0] = Cppcompiler; vars[1] = Cppcompiler.replace('icx', 'icpx')

    (cc, cxx, opt, basecflags, ccshared, ldshared, so_ext) = vars
    cc = cc.split(' ') # enleve les options si mises dans cc
    cc = cc[0]
    cxx = cxx.split(' ')
    cxx = cxx[0]
    return (cc, cxx, opt, basecflags, ccshared, ldshared, so_ext)

#==============================================================================
# Retourne le pre-processeur utilise pour les fichiers fortrans
# IN: config.Cppcompiler
#==============================================================================
def getPP():
    try: from KCore.config import Cppcompiler
    except: from config import Cppcompiler
    sizes = '-DREAL_E="REAL*8"'
    if EDOUBLEINT: sizes += ' -DINTEGER_E="INTEGER*8"'
    else: sizes += ' -DINTEGER_E="INTEGER*4"'
    if GDOUBLEINT: sizes += ' -DINTEGER_G="INTEGER*8"'
    else: sizes += ' -DINTEGER_G="INTEGER*4"'
    sizes += ' -DINTEGER_L="INTEGER*4"'
    if Cppcompiler == 'icl.exe': PP = 'fpp.exe '+sizes+' \\I'
    elif Cppcompiler == "x86_64-w64-mingw32-gcc":
        PP = 'x86_64-w64-mingw32-cpp -P -traditional %s -I'%sizes
    else: PP = 'cpp -P -traditional %s -I'%sizes
    return PP

#==============================================================================
# Retourne l'achiveur pour faire les librairies statiques
# IN: config.Cppcompiler
#==============================================================================
def getAR():
    try: from KCore.config import Cppcompiler
    except: from config import Cppcompiler
    if Cppcompiler == "icl.exe": return 'ar.exe '
    elif Cppcompiler == "x86_64-w64-mingw32-gcc":
        return 'x86_64-w64-mingw32-ar'
    else: return 'ar'

#==============================================================================
# Retourne le prefix pour le repertoire ou on stocke les modules f90
# IN: config.f90compiler
#==============================================================================
def getFortranModDirPrefix():
    try: from KCore.config import f90compiler
    except: from config import f90compiler
    if f90compiler == 'ifort': return '-module'
    elif f90compiler == 'ifort.exe': return '-module'
    elif f90compiler == 'gfortran': return '-J'
    elif f90compiler == 'g95': return '-fmod'
    else: return ''

#==============================================================================
# Retourne 1 si oui
# IN: config.useOMP
#==============================================================================
def useOMP():
    try: from KCore.config import useOMP
    except: from config import useOMP
    if useOMP: return 1
    else: return 0

#==============================================================================
# Retourne 1 si on produit des librairies statiques
# IN: config.useStatic
#==============================================================================
def useStatic():
    try: from KCore.config import useStatic
    except: from config import useStatic
    if useStatic: return 1
    else: return 0

#==============================================================================
# Retourne 1 si on dispose de cuda
# IN: config.useCuda
#==============================================================================
def useCuda():
    try: from KCore.config import useCuda
    except: from config import useCuda
    if useCuda: return 1
    else: return 0

#==============================================================================
# Retourne les versions des compilateurs
#==============================================================================
def getVersion(compiler):
    major = 0; minor = 0
    cmd = [compiler, '--version']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    out = ''
    while True:
        line = proc.stdout.readline().decode()
        if line != '': out += line
        else: break
    out = out.split('\n')
    out = out[0]
    out = out.split(' ')
    for i in out:
        isVersion = i.split('.')
        if len(isVersion)>1: # maybe
            try:
                major = int(isVersion[0])
                minor = int(isVersion[1])
                break
            except: continue
    return (major, minor)

def getCppVersion():
    try: from KCore.config import Cppcompiler
    except: from config import Cppcompiler
    return getVersion(Cppcompiler)

def getForVersion():
    try: from KCore.config import f77compiler
    except: from config import f77compiler
    return getVersion(f77compiler)

#=============================================================================
# Retourne le nbre de lignes du cache
# Se base sur l'option du compilateur C si elle contient -DCACHELINE=XX
#==============================================================================
def getCacheLine():
    opts = getCppArgs()
    for i in opts:
        if i[0:11] == '-DCACHELINE':
            val = int(i[12:]); return val
    return 1

#=============================================================================
# Retourne le niveau de vectorisation simd
# Se base sur l'option du compilateur C si elle contient -DSIMD=XX
#==============================================================================
def getSimd():
    opts = getCppArgs()
    for i in opts:
        if i[0:6] == '-DSIMD':
            val = i[7:]; return val
    return 1

# Retourne les options SIMD pour les compilateurs
# Se base sur les options precedentes qui doivent contenir -DSIMD
def getSimdOptions():
    try: from KCore.config import Cppcompiler
    except: from config import Cppcompiler
    options = getCppAdditionalOptions()
    simd = ''
    for i in options:
        if i[0:6] == '-DSIMD':
            simd = i[7:]; break
    opts = []
    if Cppcompiler.find("icpc") == 0 or Cppcompiler.find("icc") == 0 or Cppcompiler.find("icx") == 0:
        if   simd == 'SSE4.2'  : opts += ['-xSSE4.2']
        elif simd == 'AVX2'    : opts += ['-xCORE-AVX2']
        elif simd == 'AVX512'  : opts += ['-xCORE-AVX512']
        elif simd == 'MIC'     : opts += ['-xMIC-AVX512']
        elif simd == 'AVX2P512': opts += ['-axCORE-AVX2,CORE-AVX512']
        elif simd == 'XHOST'   : opts += ['-xHost']
    elif Cppcompiler.find("gcc") == 0 or Cppcompiler.find("g++") == 0:
        if   simd == 'SSE4.2'  : opts += ['-msse4.2']
        elif simd == 'AVX2'    : opts += ['-mavx2']
        elif simd == 'AVX512'  : opts += ['-mavx512f']
        elif simd == 'MIC'     : opts += ['-mavx512er']
        elif simd == 'AVX2P512': opts += ['-mavx512f']
    #print('simd', opts)
    return opts

# Retourne True si opt est une option Simd
def isSimd(opt):
    if opt[0:7] == '-xCORE-': return True
    if opt[0:7] == '-axCORE': return True
    if opt[0:6] == '-xMIC-' : return True
    if opt[0:6] == '-xHost' : return True
    if opt[0:5] == '-mavx'  : return True
    return False

#=============================================================================
# Retourne le nbre de socket du noeud
# Se base sur l'option du compilateur C si elle contient -DNB_SOCKET=XX
#==============================================================================
def getNbSocket():
    opts = getCppArgs()
    for i in opts:
        if i[0:11] == '-DNB_SOCKET':
            val = int(i[12:]); return val
    return 1

#==============================================================================
# Retourne le nbre de coeur physique par socket
# Se base sur l'option du compilateur C si elle contient -DCORE_PER_SOCK=XX
#==============================================================================
def getCorePerSocket():
    opts = getCppArgs()
    for i in opts:
        if i[0:15] == '-DCORE_PER_SOCK':
            val = int(i[16:]); return val
    return 1

#==============================================================================
# Fonction interne a getFilesOfExt
#==============================================================================
def scanext(args, dir, file):
    givenExts = args[0]
    ret = args[1]
    for f in file:
        (root,ext) = os.path.splitext(f)
        tot = '%s/%s'%(dir,f)
        t = os.path.islink(tot)
        m = True
        if f[len(f)-1] == '~' : m = False
        if 'build' in tot: m = False
        if '.svn' in tot: m = False
        if 'Stubs' in tot: m = False
        if 'test' in tot: m = False
        if 'EXEMPLES' in tot: m = False
        if ext in givenExts and not t and m:
            ret.append(tot)

#==============================================================================
# Parcours un repertoire a la recherche des fichiers d'une certaine extension
# Ex: getFilesOfExt('Converter', ['.cpp'])
#==============================================================================
def getFilesOfExt(rootdir, givenExts):
    ret = []
    os.path.walk(rootdir, scanext, [givenExts, ret])
    return ret

#==============================================================================
# Sort a file list to solve the dependance to "use module" for fortran90
#==============================================================================
def sortFileListByUse(files):
    # Calcul du nombre de USE et le nom des modules utilises
    nbs = {}; mods = {}; defmods = {}
    for f in files:
        nb = 0
        fh = open(f, 'r')
        d = fh.read()
        fh.close()
        d = d.splitlines()
        mod = []; defmod = []
        for j in d:
            tr = j[0:]; tr = tr.lstrip(); tr1 = tr[0:3]
            if tr1 == 'use' or tr1 == 'USE':
                tr2 = tr[3:]; tr2 = tr2.lstrip()
                tr2 = tr2.split(' ')
                tr2 = tr2[0]
                tr2 = tr2.split(',')
                tr2 = tr2[0]
                tr2 = tr2.split('!')
                tr2 = tr2[0]
                if tr2 not in mod:
                    mod.append(tr2)
                    nb += 1
            tr3 = tr[0:6]
            if tr3 == 'module' or tr3 == 'MODULE':
                tr2 = tr[6:]; tr2 = tr2.lstrip()
                tr2 = tr2.split(' ')
                tr2 = tr2[0]
                tr2 = tr2.split(',')
                tr2 = tr2[0]
                tr2 = tr2.split('!')
                tr2 = tr2[0]
                defmod.append(tr2)
        nbs[f] = nb
        mods[f] = mod
        defmods[f] = defmod

    # tri
    files = sorted(files, key=lambda x: nbs[x])

    out = []
    lf0 = 0
    lf = len(files)

    while lf-lf0 != 0:
        lf0 = lf
        # selectionne les fichiers a 0 dependance
        c = 0
        for f in files:
            if nbs[f] > 0: break
            c += 1

        l0 = files[0:c]
        out += l0
        files = files[c:]
        lf = len(files)

        # remplace les 0 use
        for f in files:
            for l in l0:
                for j in defmods[l]:
                    if j in mods[f]:
                        nbs[f] -= 1; mods[f].remove(j)
        # tri
        files = sorted(files, key=lambda x: nbs[x])
    if lf > 0: out += files
    #print('MODULES')
    #for f in out:
    #    if len(defmods[f])>0: print(f, defmods[f], nbs[f], '->', mods[f])
    return out

#==============================================================================
# Retourne les options additionelles des compilos definies dans config
#==============================================================================
def getCppAdditionalOptions():
    try: from KCore.config import CppAdditionalOptions
    except: from config import CppAdditionalOptions
    return CppAdditionalOptions

def getf77AdditionalOptions():
    try: from KCore.config import f77AdditionalOptions
    except: from config import f77AdditionalOptions
    return f77AdditionalOptions

#==============================================================================
# Retourne les arguments pour le compilateur C (utilise aussi pour Cpp)
# IN: config.Cppcompiler, config.useStatic, config.useOMP,
# config.CppAdditionalOptions
#==============================================================================
def getCArgs():
    try: from KCore.config import Cppcompiler
    except: from config import Cppcompiler
    mySystem = getSystem()
    compiler = Cppcompiler.split('/')
    l = len(compiler)-1
    Cppcompiler = compiler[l]
    if Cppcompiler == "None": return []
    options = getCppAdditionalOptions()[:]
    if EDOUBLEINT: options += ['-DE_DOUBLEINT']
    if GDOUBLEINT: options += ['-DG_DOUBLEINT']
    if Cppcompiler == "icpc" or Cppcompiler == "icc":
        v = getCppVersion()
        if DEBUG:
            options += ['-g', '-O0', '-wd47', '-wd1224', '-fp-trap=divzero,overflow,invalid']
        else: options += ['-DNDEBUG', '-O2', '-wd47', '-wd1224']

        # hack pour intel 19
        if v[0] == 19:
            for c, o in enumerate(options):
                if o == '-O2': options[c] = '-O1'

        if v[0] < 15:
            options += ['-fp-speculation=strict']
        else:
            options += ['-fp-model=precise'] # modif 2.6
        if useOMP() == 1:
            if v[0] < 15: options += ['-openmp']
            else: options += ['-qopenmp']
        if useStatic() == 1: options += ['-static']
        else: options += ['-fPIC']
        options += getSimdOptions()
        return options
    elif Cppcompiler.find("gcc") == 0 or Cppcompiler.find("g++") == 0:
        if DEBUG:
            options += ['-g', '-O0', '-Wall', '-pedantic', '-D_GLIBCXX_DEBUG_PEDANTIC']
            options += ['-ggdb']
            if mySystem[0] != 'mingw': # no asan on mingw
                options += ['-fsanitize=address']
                #options += ['-fsanitize=thread']
        else: options += ['-DNDEBUG', '-O3', '-Wall', '-Werror=return-type']
        if useOMP() == 1: options += ['-fopenmp']
        if useStatic() == 1: options += ['--static', '-static-libstdc++', '-static-libgcc']
        else: options += ['-fPIC']
        if mySystem[0] == 'mingw' and mySystem[1] == '32':
            options.remove('-fPIC')
            options += ['-large-address-aware']
        options += getSimdOptions()
        return options
    elif Cppcompiler == "icl.exe":
        options += ['/EHsc', '/MT']
        if useOMP() == 1: options += ['/Qopenmp']
        return options
    elif Cppcompiler == "icx" or Cppcompiler == "icpx":
        if DEBUG: options += ['-g', '-O0', '-fp-trap=divzero,overflow,invalid']
        else: options += ['-DNDEBUG', '-O2',]
        options += ['-fp-model=precise'] # existe encore?
        if useOMP() == 1: options += ['-qopenmp']
        if useStatic() == 1: options += ['-static']
        else: options += ['-fPIC']
        options += getSimdOptions()
        return options
    elif Cppcompiler == "pgcc" or Cppcompiler == "pgc++":
        if DEBUG: options += ['-g', '-O0']
        else: options += ['-DNDEBUG', '-O3']
        if useOMP() == 1: options += ['-mp=multicore']
        else: options += ['-nomp']
        if useStatic() == 1: options += []
        else: options += ['-fPIC']
        options += getSimdOptions()
        if useCuda(): options += ['-acc=gpu', '-Minfo:accel']
        return options
    elif Cppcompiler == "nvc" or Cppcompiler == "nvc++":
        if DEBUG: options += ['-g', '-O0']
        else: options += ['-DNDEBUG', '-O3']
        if useOMP() == 1: options += ['-mp=multicore']
        else: options += ['-nomp']
        if useStatic() == 1: options += ['-static']
        else: options += ['-fPIC']
        options += getSimdOptions()
        if useCuda(): options += ['-acc=gpu', '-Minfo:accel']
        return options
    elif Cppcompiler == "x86_64-w64-mingw32-gcc" or Cppcompiler == "x86_64-w64-mingw32-g++":
        options += ['-DMS_WIN64', '-fpermissive', '-D__USE_MINGW_ANSI_STDIO=1']
        if DEBUG: options += ['-g', 'O0', '-D_GLIBCXX_DEBUG_PEDANTIC']
        else: options += ['-DNDEBUG', '-O3']
        if useOMP() == 1: options += ['-fopenmp']
        if useStatic() == 1: options += ['--static', '-static-libstdc++', '-static-libgcc']
        else: options += ['-fPIC']
        options += getSimdOptions()
        return options
    elif Cppcompiler == "clang" or Cppcompiler == "clang++":
        if DEBUG: options += ['-g', '-O0', '-Wall', '-D_GLIBCXX_DEBUG_PEDANTIC']
        else: options += ['-DNDEBUG', '-O3', '-Wall']
        if useOMP() == 1: options += ['-fopenmp']
        if useStatic() == 1: options += ['--static', '-static-libstdc++', '-static-libgcc']
        else: options += ['-fPIC']
        options += getSimdOptions()
        return options
    elif Cppcompiler == "craycc" or Cppcompiler == "craycxx":
        if DEBUG: options += ['-g', '-O0', '-Wall', '-D_GLIBCXX_DEBUG_PEDANTIC']
        else: options += ['-DNDEBUG', '-O3', '-Wall']
        if useOMP() == 1: options += ['-fopenmp']
        if useStatic() == 1: options += ['--static', '-static-libstdc++', '-static-libgcc']
        else: options += ['-fPIC']
        options += getSimdOptions()
        return options
    elif Cppcompiler == "cc":
        if DEBUG: options += ['-g', '-O0', '-Wall', '-D_GLIBCXX_DEBUG_PEDANTIC']
        else: options += ['-DNDEBUG', '-O3', '-Wall']
        if useOMP() == 1: options += ['-fopenmp']
        if useStatic() == 1: options += ['--static', '-static-libstdc++', '-static-libgcc']
        else: options += ['-fPIC']
        options += getSimdOptions()
        return options
    else: return options

# Options pour le compilateur C++
def getCppArgs():
    opt = getCArgs()
    try: from KCore.config import Cppcompiler
    except: from config import Cppcompiler
    if Cppcompiler == "icl.exe": opt += ["/std=c++11"]
    else: opt += ["-std=c++11"]
    return opt

#==============================================================================
# Retourne les arguments pour le compilateur Cuda
# IN: config.Cppcompiler, config.useStatic, config.useOMP,
# config.CppAdditionalOptions
#==============================================================================
def getCudaArgs():
    try: from KCore.config import NvccAdditionalOptions
    except: from config import NvccAdditionalOptions
    options = NvccAdditionalOptions
    if DEBUG: options += ['-g', '-O0']
    else: options += ['-DNDEBUG', '-O2']
    return options

#==============================================================================
# Retourne les arguments pour le compilateur Fortran
# IN: config.f77compiler
#==============================================================================
def getForArgs():
    try: from KCore.config import f77compiler
    except: from config import f77compiler
    mySystem = getSystem()
    compiler = f77compiler.split('/')
    l = len(compiler)-1
    f77compiler = compiler[l]
    if f77compiler == "None": return []
    options = getf77AdditionalOptions()
    if f77compiler == "gfortran":
        if DEBUG: options += ['-Wall', '-g', '-O0', '-fbacktrace', '-fbounds-check', '-ffpe-trap=zero,overflow,invalid']
        else: options += ['-Wall', '-O3']
        if useOMP() == 1: options += ['-fopenmp']
        if useStatic() == 1: options += ['--static']
        else: options += ['-fPIC']
        if mySystem[0] == 'mingw' and mySystem[1] == '32':
            options.remove('-fPIC')
            options += ['-large-address-aware']
        options += getSimdOptions()
        if EDOUBLEINT: options += ['-fdefault-integer-8']
        options += ['-fdefault-real-8', '-fdefault-double-8']
        return options
    elif f77compiler == "ifort":
        if DEBUG: options += ['-g', '-O0', '-CB', '-traceback', '-fpe0']
        else: options += ['-O3']
        v = getForVersion()
        if v[0] < 15: options += ['-fp-speculation=strict']
        else: options += ['-fp-model=precise']
        if useOMP() == 1:
            v = getForVersion()
            if v[0] < 15: options += ['-openmp']
            else: options += ['-qopenmp']
        if useStatic() == 1: options += ['-static']
        else: options += ['-fPIC']
        options += getSimdOptions()
        if EDOUBLEINT: options += ['-i8']
        options += ['-r8']
        return options
    elif f77compiler == "ifx":
        if DEBUG: options += ['-g', '-O0', '-CB', '-fpe0']
        else: options += ['-O3']
        options += ['-fp-model=precise']
        if useOMP() == 1: options += ['-qopenmp']
        if useStatic() == 1: options += ['-static']
        else: options += ['-fPIC']
        options += getSimdOptions()
        if EDOUBLEINT: options += ['-i8']
        options += ['-r8']
        return options
    elif f77compiler == "pgfortran":
        if useStatic() == 1: options += ['-static']
        else: options += ['-fPIC']
        if DEBUG: options += ['-g', '-O0']
        else: options += ['-O3']
        if useOMP() == 1: options += ['-mp=multicore']
        options += getSimdOptions()
        if EDOUBLEINT: options += ['-i8']
        options += ['-r8']
        return options
    elif f77compiler == "nvfortran":
        if useStatic() == 1: options += ['-static']
        else: options += ['-fPIC']
        if DEBUG: options += ['-g', '-O0']
        else: options += ['-O3']
        if useOMP() == 1: options += ['-mp=multicore']
        options += getSimdOptions()
        if EDOUBLEINT: options += ['-i8']
        options += ['-r8']
        return options
    elif f77compiler == "x86_64-w64-mingw32-gfortran":
        if DEBUG: options += ['-g', '-O0']
        else: options += ['-O3']
        if useOMP() == 1: options += ['-fopenmp']
        if useStatic() == 1: options += ['--static']
        else: options += ['-fPIC']
        options += getSimdOptions()
        if EDOUBLEINT: options += ['-fdefault-integer-8']
        options += ['-fdefault-real-8', '-fdefault-double-8']
        return options
    elif f77compiler == "ifort.exe":
        if useOMP() == 1: return ['/names:lowercase', '/assume:underscore', '/Qopenmp']
        else: return ['/names:lowercase', '/assume:underscore']
    elif f77compiler == "crayftn":
        if useStatic() == 1: options += ['-static']
        else: options += ['-fPIC']
        if DEBUG: options += ['-g', '-O0']
        else: options += ['-O3']
        if useOMP() == 1: options += ['-fopenmp']
        options += getSimdOptions()
        if EDOUBLEINT: options += ['-i8']
        options += ['-r8']
        return options
    elif f77compiler == "ftn":
        if useStatic() == 1: options += ['-static']
        else: options += ['-fPIC']
        if DEBUG: options += ['-g', '-O0']
        else: options += ['-O3']
        if useOMP() == 1: options += ['-fopenmp']
        options += getSimdOptions()
        if EDOUBLEINT: options += ['-i8']
        options += ['-r8']
        return options
    elif f77compiler == "flang":
        if useStatic() == 1: options += ['-static']
        else: options += ['-fPIC']
        if DEBUG: options += ['-g', '-O0']
        else: options += ['-O3']
        if useOMP() == 1: options += ['-fopenmp']
        options += getSimdOptions()
        if EDOUBLEINT: options += ['-fdefault-integer-8']
        options += ['-fdefault-real-8', '-fdefault-double-8']
        return options
    else: return options

#==============================================================================
# Retourne les arguments pour le linker
# IN: config.Cppcompiler, config.useStatic, config.useOMP
#==============================================================================
def getLinkArgs():
    try: from KCore.config import Cppcompiler
    except: from config import Cppcompiler
    out = []
    if Cppcompiler == 'gcc' or Cppcompiler == 'g++':
        if useStatic() == 1: out += ['--static']
    elif Cppcompiler == 'icc' or Cppcompiler == 'icpc':
        if useStatic() == 1: out += ['-static']
    elif Cppcompiler == 'icx' or Cppcompiler == 'icpx':
        if useStatic() == 1: out += ['-static']
        else: out += ['-shared']
    elif Cppcompiler == "x86_64-w64-mingw32-gcc":
        if useStatic() == 1: out += ['--static']
    elif Cppcompiler == 'pgcc' or Cppcompiler == 'pgc++':
        if useStatic() == 1: out += ['-static']
        else: out += ['-shared']
        if useOMP() == 1: out += ['-mp=multicore']
        if useCuda() == 1: out += ['-acc=gpu', '-Minfo:accel']
    elif Cppcompiler == 'nvc' or Cppcompiler == 'nvc++':
        if useStatic() == 1: out += ['-static']
        else: out += ['-shared']
        if useOMP() == 1: out += ['-mp=multicore']
        if useCuda() == 1: out += ['-acc=gpu', '-Minfo:accel']
    elif Cppcompiler == 'craycc' or Cppcompiler == 'craycxx':
        if useStatic() == 1: out += ['-static']
        else: out += ['-shared']
        if useOMP() == 1: out += ['-fopenmp']
    elif Cppcompiler == 'cc':
        if useStatic() == 1: out += ['-static']
        else: out += ['-shared']
        if useOMP() == 1: out += ['-fopenmp']
    elif Cppcompiler == 'clang' or Cppcompiler == 'clang++':
        if useStatic() == 1: out += ['-static']
        else: out += ['-shared']
        if useOMP() == 1: out += ['-fopenmp']
    mySystem = getSystem()[0]
    if mySystem == 'Darwin':
        if useStatic() == 0: out += ['-dynamiclib']
    return out

#=============================================================================
# Check PYTHONPATH
# Verifie que installPath est dans PYTHONPATH
#=============================================================================
def checkPythonPath():
    import re
    try: import KCore.installPath as K
    except: import installPath as K
    installPathLocal = K.installPath
    a = os.getenv("PYTHONPATH")
    if a is None:
        print('Warning: to use the module, please add: %s to your PYTHONPATH.'%installPathLocal)
    else:
        if re.compile(installPathLocal).search(a) is None:
            print('Warning: to use the module, please add: %s to your PYTHONPATH.'%installPathLocal)

#=============================================================================
# Check LD_LIBRARY_PATH
#=============================================================================
def checkLdLibraryPath():
    import re
    try: import KCore.installPath as K
    except: import installPath as K
    libPath = K.libPath
    a = os.getenv("LD_LIBRARY_PATH")
    b = os.getenv("LIBRARY_PATH")
    if a is None and b is None:
        print("Warning: to use the module, please add: %s to your LD_LIBRARY_PATH (unix) or PATH (windows)."%libPath)
    else:
        if a is not None: ret = a
        else: ret = b
        if re.compile(libPath).search(ret) is None:
            print("Warning: to use the module, please add: %s to your LD_LIBRARY_PATH (unix) or PATH (windows)."%libPath)

#=============================================================================
# Check for KCore (Cassiopee core)
#=============================================================================
def checkKCore():
    try:
        import KCore
        import KCore.installPath
        kcoreIncDir = KCore.installPath.includePath
        kcoreIncDir = os.path.join(kcoreIncDir, 'KCore')
        kcoreLibDir = KCore.installPath.libPath
        return (KCore.__version__, kcoreIncDir, kcoreLibDir)

    except ImportError:
        raise SystemError("Error: kcore library is required for the compilation of this module.")

#=============================================================================
# Check for XCore (Parallel core)
#=============================================================================
def checkXCore():
    try:
        import XCore
        import KCore.installPath
        xcoreIncDir = KCore.installPath.includePath
        xcoreIncDir = os.path.dirname(xcoreIncDir)
        xcoreIncDir = os.path.join(xcoreIncDir, 'XCore/XCore')
        xcoreLibDir = KCore.installPath.libPath
        return (XCore.__version__, xcoreIncDir, xcoreLibDir)

    except ImportError:
        raise SystemError("Error: xcore library is required for the compilation of this module.")

#=============================================================================
# Check for Generator
#=============================================================================
def checkGenerator():
    try:
        import Generator
        import KCore.installPath
        generatorIncDir = KCore.installPath.includePath
        generatorIncDir = os.path.dirname(generatorIncDir)
        generatorIncDir = os.path.join(generatorIncDir, 'Generator/Generator')
        generatorLibDir = KCore.installPath.libPath
        return (Generator.__version__, generatorIncDir, generatorLibDir)

    except ImportError:
        raise SystemError("Error: generator library is required for the compilation of this module.")

#=============================================================================
# Check for FastC module
#=============================================================================
def checkFastC():
    try:
        import FastC.installPath
        import KCore.installPath
        fastcIncDir = FastC.installPath.includePath
        fastcLibDir = KCore.installPath.libPath
        return (FastC.__version__, fastcIncDir, fastcLibDir)

    except ImportError:
        raise SystemError("Error: fastc library is required for the compilation of this module.")

#=============================================================================
# Check for FastS module
#=============================================================================
def checkFastS():
    try:
        import FastS
        import FastC.installPath
        import KCore.installPath
        fastsIncDir = FastC.installPath.includePath
        fastsIncDir = os.path.dirname(fastsIncDir)
        fastsIncDir = os.path.join(fastsIncDir, 'FastS')
        fastsLibDir = KCore.installPath.libPath
        return (FastS.__version__, fastsIncDir, fastsLibDir)

    except ImportError:
        raise SystemError("Error: fasts library is required for the compilation of this module.")

#=============================================================================
# Check for FastP module
#=============================================================================
def checkFastP():
    try:
        import FastP
        import FastC.installPath
        import KCore.installPath
        fastpIncDir = FastC.installPath.includePath
        fastpIncDir = os.path.dirname(fastpIncDir)
        fastpIncDir = os.path.join(fastpIncDir, 'FastP')
        fastpLibDir = KCore.installPath.libPath
        return (FastP.__version__, fastpIncDir, fastpLibDir)

    except ImportError:
        raise SystemError("Error: fastp library is required for the compilation of this module.")

#=============================================================================
# Check for FastLBM module
#=============================================================================
def checkFastLBM():
    try:
        import FastLBM
        import FastC.installPath
        import KCore.installPath
        fastlbmIncDir = FastC.installPath.includePath
        fastlbmIncDir = os.path.dirname(fastlbmIncDir)
        fastlbmIncDir = os.path.join(fastlbmIncDir, 'FastLBM')
        fastlbmLibDir = KCore.installPath.libPath
        return (FastLBM.__version__, fastlbmIncDir, fastlbmLibDir)

    except ImportError:
        raise SystemError("Error: fastlbm library is required for the compilation of this module.")

#=============================================================================
# Check for FastASLBM module
#=============================================================================
def checkFastASLBM():
    try:
        import FastASLBM
        import FastC.installPath
        import KCore.installPath
        fastaslbmIncDir = FastC.installPath.includePath
        fastaslbmIncDir = os.path.dirname(fastaslbmIncDir)
        fastaslbmIncDir = os.path.join(fastaslbmIncDir, 'FastASLBM')
        fastaslbmLibDir = KCore.installPath.libPath
        return (FastASLBM.__version__, fastaslbmIncDir, fastaslbmLibDir)

    except ImportError:
        raise SystemError("Error: fastaslbm library is required for the compilation of this module.")

#=============================================================================
# Check for Connector module
#=============================================================================
def checkConnector():
    try:
        import Connector
        import KCore.installPath
        ConnectorIncDir = KCore.installPath.includePath
        ConnectorIncDir = os.path.dirname(ConnectorIncDir)
        ConnectorIncDir = os.path.join(ConnectorIncDir, 'Connector/Connector')
        ConnectorLibDir = KCore.installPath.libPath
        return (Connector.__version__, ConnectorIncDir, ConnectorLibDir)

    except ImportError:
        raise SystemError("Error: Connector library is required for the compilation of this module.")

#=============================================================================
# Check for Cassiopee Kernel in Dist/ or Kernel/
# La prod de Cassiopee ne doit pas contenir mpi
#=============================================================================
def checkCassiopee():
    Cassiopee = False
    CassiopeeLibDir = ""
    CassiopeeIncDir = ""
    CassiopeeUseMpi = False # si le Kernel utilise mpi

    kvar = os.getenv("CASSIOPEE")
    pvar = os.getenv("ELSAPROD") # prod de Cassiopee
    if kvar is None:
        return (Cassiopee, CassiopeeIncDir, CassiopeeLibDir, CassiopeeUseMpi)
    if pvar is None:
        return (Cassiopee, CassiopeeIncDir, CassiopeeLibDir, CassiopeeUseMpi)

    # Cassiopee must exists and be installed
    # si Kernel existe, on l'utilise en priorite
    # sinon, on utilise Dist
    a1 = os.access(kvar+"/Kernel/include", os.F_OK)
    a2 = os.access(kvar+"/Kernel/lib/"+pvar+'/elsAc.so', os.F_OK)
    b1 = os.access(kvar+"/Dist/include", os.F_OK)
    b2 = os.access(kvar+"/Dist/bin/"+pvar+'/elsAc.so', os.F_OK)

    if a1 and a2:
        Cassiopee = True
        CassiopeeIncDir = kvar + "/Kernel/include"
        CassiopeeLibDir = kvar + "/Kernel/lib/" + pvar

    elif b1 and b2:
        Cassiopee = True
        CassiopeeIncDir = kvar + "/Dist/include"
        CassiopeeLibDir = kvar + "/Dist/bin/" + pvar

    if not Cassiopee: # essai avec une production mpi
        pvar = pvar.split('_')
        if len(pvar) >= 2:
            pvar = pvar[0]+'_mpi_'+pvar[1]
            a1 = os.access(kvar+"/Kernel/include", os.F_OK)
            a2 = os.access(kvar+"/Kernel/lib/"+pvar+'/elsAc.so', os.F_OK)
            b1 = os.access(kvar+"/Dist/include", os.F_OK)
            b2 = os.access(kvar+"/Dist/bin/"+pvar+'/elsAc.so', os.F_OK)
            if a1 and a2:
                Cassiopee = True
                CassiopeeUseMpi = True
                CassiopeeIncDir = kvar + "/Kernel/include"
                CassiopeeLibDir = kvar + "/Kernel/lib/" + pvar

            elif b1 and b2:
                Cassiopee = True
                CassiopeeUseMpi = True
                CassiopeeIncDir = kvar + "/Dist/include"
                CassiopeeLibDir = kvar + "/Dist/bin/" + pvar

    if Cassiopee:
        print('Info: Cassiopee Kernel detected at '+CassiopeeLibDir+'.')
        print('Info: .Cassiopee extension will be built.')

    return (Cassiopee, CassiopeeIncDir, CassiopeeLibDir, CassiopeeUseMpi)

#=============================================================================
# Check for elsA Kernel in Dist/ or Kernel/
#=============================================================================
def checkElsa():
    elsA = False
    elsALibDir = ""
    elsAIncDir = ""
    elsAUseMpi = False

    kvar = os.getenv("ELSA")
    if kvar is None: return (elsA, elsAIncDir, elsALibDir, elsAUseMpi)
    pvar = os.getenv("ELSAPROD")
    if pvar is None: return (elsA, elsAIncDir, elsALibDir, elsAUseMpi)

    a1 = os.access(kvar+"/Kernel/include", os.F_OK)
    a2 = os.access(kvar+"/Kernel/lib/"+pvar+'/elsA.x', os.F_OK)
    b1 = os.access(kvar+"/Dist/include", os.F_OK)
    b2 = os.access(kvar+"/Dist/bin/"+pvar+'/elsA.x', os.F_OK)
    if a1 and a2:
        elsA = True
        elsAIncDir = kvar + "/Kernel/include"
        elsALibDir = kvar + "/Kernel/lib/" + pvar

    elif b1 and b2:
        elsA = True
        elsAIncDir = kvar + "/Dist/include"
        elsALibDir = kvar + "/Dist/bin/" + pvar

    if not elsA: # essai avec une production mpi
        pvar = pvar.split('_')
        if len(pvar) >= 2:
            pvar = pvar[0]+'_mpi_'+pvar[1]
            a1 = os.access(kvar+"/Kernel/include", os.F_OK)
            a2 = os.access(kvar+"/Kernel/lib/"+pvar+'/elsA.x', os.F_OK)
            b1 = os.access(kvar+"/Dist/include", os.F_OK)
            b2 = os.access(kvar+"/Dist/bin/"+pvar+'/elsAc.x', os.F_OK)
            if a1 and a2:
                elsA = True
                elsAUseMpi = True
                elsAIncDir = kvar + "/Kernel/include"
                elsALibDir = kvar + "/Kernel/lib/" + pvar

            elif b1 and b2:
                elsA = True
                elsAUseMpi = True
                elsAIncDir = kvar + "/Dist/include"
                elsALibDir = kvar + "/Dist/bin/" + pvar

    if elsA:
        print('Info: elsA Kernel detected at '+elsALibDir+'.')
        print('Info: .Elsa extension will be built.')
    return (elsA, elsAIncDir, elsALibDir, elsAUseMpi)

#=============================================================================
# Check for GL (libGL)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
#=============================================================================
def checkGL(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libGL.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libGL.a', additionalLibPaths)
    i = checkIncFile__('GL/gl.h', additionalIncludePaths)
    if i is not None and l is not None:
        print('Info: openGL detected at '+l+'.')
        return (True, i, l)
    else:
        print('Info: openGL or GL/gl.h was not found on your system.')
        return (False, '', '')

#=============================================================================
# Check for Glut (libglut)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
#=============================================================================
def checkGlut(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libglut.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libglut.a', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libfreeglut.a', additionalLibPaths)
    i = checkIncFile__('GL/glut.h', additionalIncludePaths)
    if i is not None and l is not None:
        print('Info: glut detected at '+l+'.')
        return (True, i, l)
    else:
        print('Info: libglut or GL/glut.h was not found on your system.')
        return (False, '', '')

#=============================================================================
# Check for Glew (libglew)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkGlew(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libGLEW.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libGLEW.a', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libglew32.a', additionalLibPaths)
    i = checkIncFile__('GL/glew.h', additionalIncludePaths)

    if i is not None and l is not None:
        print('Info: glew detected at '+l+'.')
        return (True, i, l)
    else:
        print('Info: libglew or GL/glew.h was not found on your system. No shader support for CPlot.')
        return (False, '', '')

#=============================================================================
# Check for osmesa (offline rendering)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkOSMesa(additionalLibPaths=[], additionalIncludePaths=[]):
    libname = 'OSMesa'
    l = checkLibFile__('libOSMesa.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libOSMesa.a', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libosmesa.a', additionalLibPaths)
        if l is not None: libname = 'osmesa'
    if l is None:
        l = checkLibFile__('osmesa.dll.a', additionalLibPaths)
        if l is not None: libname = 'osmesa'

    i = checkIncFile__('GL/osmesa.h', additionalIncludePaths)
    if i is not None and l is not None:
        print('Info: libOSmesa detected at %s.'%l)
        return (True, i, l, libname)
    elif l is None and i is not None:
        print('Info: libOSMesa was not found on your system. No OSMESA offscreen support for CPlot.')
        return (False, i, l, libname)
    elif l is not None and i is None:
        print('Info: GL/osmesa.h was not found on your system. No OSMESA offscreen support for CPlot.')
        return (False, i, l, libname)
    else:
        print('Info: libOSMesa and GL/osmesa.h was not found on your system. No OSMESA offscreen support for CPlot.')
        return (False, i, l, libname)

#=============================================================================
# Check for OCE (open cascade edition library)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkOCE(additionalLibPaths=[], additionalIncludePaths=[]):
    #print("INFO: dependance to OCE STUBED.")
    #return (False, None, None)
    l = checkLibFile__('libTKernel.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libTKernel.a', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libTKernel.dll.a', additionalLibPaths)
    i = checkIncFile__('oce/TopTools.hxx', additionalIncludePaths)
    if i is not None: i = i+'/oce'
    if i is None:
        i = checkIncFile__('opencascade/TopTools.hxx', additionalIncludePaths)
        if i is not None: i = i+'/opencascade'
    if i is not None and l is not None:
        print('Info: libOCE detected at %s.'%l)
        return (True, i, l)
    else:
        # On n'affiche pas ici le message, car il peut y avoir un installation locale de OCE
        #print('Info: libOCE or oce/*.hxx was not found on your system. No IGES/STEP support.')
        return (False, i, l)

#=============================================================================
# Check for png (libpng)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkPng(additionalLibPaths=[], additionalIncludePaths=[]):
    #print("INFO: dependance to PNG STUBED.")
    #return (False, None, None)
    l = checkLibFile__('libpng.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libpng.a', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libpng.dll.a', additionalLibPaths)
    i = checkIncFile__('png.h', additionalIncludePaths)
    if i is not None and l is not None:
        print('Info: png detected at %s.'%l)
        return (True, i, l)
    elif l is None and i is not None:
        print('Info: libpng was not found on your system. No png support.')
        return (False, i, l)
    elif l is not None and i is None:
        print('Info: png.h was not found on your system. No png support.')
        return (False, i, l)
    else:
        print('Info: libpng and png.h was not found on your system. No png support.')
        return (False, i, l)

#=============================================================================
# Check for mpeg (ffmpeg)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkMpeg(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libavcodec.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libavcodec.a', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libavcodec.dll.a', additionalLibPaths)
    i = checkIncFile__('libavcodec/avcodec.h', additionalIncludePaths)
    if i is not None:
        i = checkIncFile__('libavutil/mem.h', additionalIncludePaths)
    if i is not None:
        i = checkIncFile__('libavutil/imgutils.h', additionalIncludePaths)
    if i is not None and l is not None:
        print('Info: mpeg detected at %s.'%l)
        return (True, i, l)
    elif l is None and i is not None:
        print('Info: libavcodec was not found on your system. No mpeg support.')
        return (False, i, l)
    elif l is not None and i is None:
        print('Info: libavcodec/avcodec.h,  libavutil/mem.h or libavutil/imgutils.h was not found on your system. No mpeg support.')
        return (False, i, l)
    else:
        print('Info: libavcodec and libavcodec/avcodec.h,  libavutil/mem.h or libavutil/imgutils.h was not found on your system. No mpeg support.')
        return (False, i, l)

#=============================================================================
# Check for Adf
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkAdf(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libcgns.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libcgns.a', additionalLibPaths)
    i = checkIncFile__('adf/ADF.h', additionalIncludePaths)
    if i is not None and l is not None:
        print('Info: Adf detected at %s.'%l)
        return (True, i, l)
    else:
        print('Info: libadf or adf/ADF.h was not found on your system. No adf support.')
        return (False, i, l)

#=============================================================================
# Check for Hdf
# additionalPaths: chemins d'installation non standards : ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie,
# liste des noms des libs)
#=============================================================================
def checkHdf(additionalLibPaths=[], additionalIncludePaths=[]):
    libnames = []
    l = checkLibFile__('libhdf5.so', additionalLibPaths)
    if l is not None: libnames.append('hdf5')
    if l is None:
        l = checkLibFile__('libhdf5.a', additionalLibPaths)
        if l is not None: libnames.append('hdf5')
    if l is None:
        l = checkLibFile__('libhdf5_openmpi.so', additionalLibPaths)
        if l is not None: libnames.append('hdf5_openmpi')
    if l is None:
        l = checkLibFile__('libhdf5_serial.so', additionalLibPaths)
        if l is not None: libnames.append('hdf5_serial')
    i = checkIncFile__('hdf5.h', additionalIncludePaths)
    if i is not None and l is not None:
        print('Info: Hdf5 detected at %s.'%l)
        return (True, i, l, libnames)
    elif l is None and i is not None:
        print('Info: libhdf5 was not found on your system. No hdf5 support.')
        return (False, i, l, libnames)
    elif l is not None and i is None:
        print('Info: hdf5.h was not found on your system. No hdf5 support.')
        return (False, i, l, libnames)
    else:
        print('Info: libhdf5 and hdf5.h was not found on your system. No hdf5 support.')
        return (False, i, l, libnames)

#=============================================================================
# Check for netcdf
# additionalPaths: chemins d'installation non standards : ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie,
# liste des noms des libs)
#=============================================================================
def checkNetcdf(additionalLibPaths=[], additionalIncludePaths=[]):
    libnames = []
    l = checkLibFile__('libnetcdf.so', additionalLibPaths)
    if l is not None: libnames.append('netcdf')
    if l is None:
        l = checkLibFile__('libnetcdf.a', additionalLibPaths)
        if l is not None: libnames.append('netcdf')
    i = checkIncFile__('netcdf.h', additionalIncludePaths)
    if i is not None and l is not None:
        print('Info: netcdf detected at %s.'%l)
        return (True, i, l, libnames)
    elif l is None and i is not None:
        print('Info: libhnetcdf was not found on your system. No netcdf support.')
        return (False, i, l, libnames)
    elif l is not None and i is None:
        print('Info: netcdf.h was not found on your system. No netcdf support.')
        return (False, i, l, libnames)
    else:
        print('Info: libnetcdf and netcdf.h was not found on your system. No netcdf support.')
        return (False, i, l, libnames)

#=============================================================================
# Check for Mpi
# additionalPaths: chemins d'installation non standards : ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie, nom des libs)
#=============================================================================
def checkMpi(additionalLibPaths=[], additionalIncludePaths=[]):
    #print("INFO: dependance to MPI STUBED.")
    #return (False, None, None, None)
    libnames = []
    l = checkLibFile__('libmpi.so', additionalLibPaths)
    if l is not None: libnames.append('mpi')
    if l is None:
        l = checkLibFile__('libmpi.a', additionalLibPaths)
        if l is not None: libnames.append('mpi')
    if l is None:
        l = checkLibFile__('msmpi.lib', additionalLibPaths)
        if l is not None: libnames.append('msmpi')
    if l is None:
        l = checkLibFile__('libmsmpi.dll.a', additionalLibPaths)
        if l is not None: libnames.append('msmpi')

    # pour open mpi
    o = checkLibFile__('libmpi_cxx.so', additionalLibPaths)
    if o is not None and o == l: libnames.append('mpi_cxx')

    i = checkIncFile__('mpi.h', additionalIncludePaths)

    if i is not None and l is not None:
        print('Info: Mpi detected at %s.'%l)
        return (True, i, l, libnames)
    elif l is None and i is not None:
        print('Info: libmpi/msmpi was not found on your system. No Mpi support.')
        return (False, i, l, libnames)
    elif l is not None and i is None:
        print('Info: mpi.h was not found on your system. No Mpi support.')
        return (False, i, l, libnames)
    else:
        print('Info: libmpi/msmpi and mpi.h was not found on your system. No Mpi support.')
        return (False, i, l, libnames)

#=============================================================================
# Check for Mpi4py
# additionalPaths: chemins d'installation non standards : ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkMpi4py(additionalLibPaths=[], additionalIncludePaths=[]):
    #print("INFO: dependance to MPI STUBED.")
    #return (False, None, None, None)
    try: import mpi4py
    except:
        print('Info: mpi4py or mpi4py.MPI.h was not found on your system. No Mpi support.')
        return (False, '', '')

    incPaths = []
    try: import KCore.installPath as K
    except: import installPath as K
    incPaths += [K.installPath+'/mpi4py/include']
    fileN = mpi4py.__file__
    incPaths += [os.path.dirname(fileN)+'/include']
    i = checkIncFile__('mpi4py/mpi4py.MPI.h', additionalIncludePaths+incPaths)
    if i is None:
        i = checkIncFile__('mpi4py/mpi4py.h', additionalIncludePaths+incPaths)

    if i is not None:
        print('Info: Mpi4py detected at %s.'%i)
        return (True, i, '')
    else:
        print('Info: mpi4py.MPI.h or mpi4py.h was not found on your system. No Mpi support.')
        return (False, i, '')

#=============================================================================
# Check for Cuda
# additionalPaths : chemins d'installation non standards : ['/home/toto', ...]
# Retourne : (True/False, chemin des includes, chemin de la librairie, chemin
#             executable nvcc)
#=============================================================================
def checkCuda(additionalLibPaths=[], additionalIncludePaths=[]):
    # Check if the user want cuda supported
    # -------------------------------------
    try: from KCore.config import useCuda
    except: from config import useCuda
    if not useCuda:
        #print('Info: cuda is not activated. No cuda support.')
        return (False, None, None, None, None)
    # Check for executable
    # --------------------
    cuda_root = os.getenv("CUDA_ROOT")
    if cuda_root is None:
        print('Info: CUDA_ROOT environment variable not set. No cuda support.')
        return (False, None, None, None, None)
    has_nvcc_exe = os.access(cuda_root+"/bin/nvcc.exe", os.F_OK)
    has_nvcc     = os.access(cuda_root+"/bin/nvcc", os.F_OK)
    nvcc_exec    = None
    if not has_nvcc and not has_nvcc_exe :
        print('Info: nvcc not found at %s/bin. No cuda support.'%cuda_root)
        return (False, None, None, None, None)
    elif has_nvcc: nvcc_exec = cuda_root+"/bin/nvcc"
    else: nvcc_exec = cuda_root+"/bin/nvcc.exe"

    # Check for library :
    # ------------------
    libPaths = [cuda_root+'/lib/x64', cuda_root+'/lib/win32',
                cuda_root+'/lib64'] + additionalLibPaths
    libnames = []
    l = checkLibFile__('libcuda.so', libPaths)
    if l is None:
        l = checkLibFile__('cuda.lib', libPaths)
        if l is not None:
            libnames.append('cuda.lib')
    else:
        libnames.append('libcuda.so')
    # Check for include :
    # -------------------
    incPaths = [cuda_root+'/include',] + additionalIncludePaths
    incnames = []
    i = checkIncFile__("cuda.h", incPaths)
    if i is not None and l is not None:
        print('Info: cuda detected at %s.'%l)
        return (True, i, l, libnames, nvcc_exec)
    elif l is None and i is not None:
        print('Info: libcuda/cuda.lib was not found on your system. No cuda support.')
        return (False, i, l, libnames, nvcc_exec)
    elif l is not None and i is None:
        print('Info: cuda.h was not found on your system. No cuda support.')
        return (False, i, l, libnames, nvcc_exec)
    else:
        print('Info: libcuda/cuda.lib and cuda.h was not found on your system. No cuda support.')
        return (False, i, l, libnames, nvcc_exec)

    return (False, i, l, libnames, nvcc_exec)

#=============================================================================
# Check for paradigma python binding
# additionalPaths: chemins d'installation non standards : ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkPyParadigma(additionalLibPaths=[], additionalIncludePaths=[]):
    try: import XCore.Pypdm
    except:
        print('Info: python module for paradigma was not found on your system.')
        return (False, '', '')

    i = XCore.Pypdm.__file__
    print('Info: python module for paradigma detected at %s.'%i)
    return (True, i, '')

#=============================================================================
# Check for paradigma
# additionalPaths: chemins d'installation non standards : ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkParadigma(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('lipdm.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libpdm.a', additionalLibPaths)
    i = checkIncFile__('pdm.h', additionalIncludePaths)
    if i is not None and l is not None:
        #print('Info: Paradigma detected at %s.'%l)
        return (True, i, l)
    else:
        #print('Info: libpdm or pdm.h was not found on your system. No paradigma support.')
        return (False, i, l)

#=============================================================================
# Check for BLAS
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie,
# option de compile, nom de la librarie)
#=============================================================================
def checkBlas(additionalLibPaths=[], additionalIncludePaths=[]):
    try: from KCore.config import Cppcompiler
    except: from config import Cppcompiler
    if Cppcompiler == 'icc' or Cppcompiler == 'icpc' or Cppcompiler == 'icx': # intel - cherche dans MKL
        libPrefix = 'libmkl_'; includePrefix = 'mkl_'; compOpt = '-mkl'
    else: # cherche std
        libPrefix = 'lib'; includePrefix = ''; compOpt = ''

    # Search for openblas than contains blas generally
    libnames = []
    libname = 'openblas'
    l = checkLibFile__(libPrefix+'openblas*.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__(libPrefix+'openblas*.a', additionalLibPaths)
    if l is not None: libnames.append(libname)

    if l is None: # search for blas and cblas
        libname = 'blas'; foundBlas = None
        l = checkLibFile__(libPrefix+'blas*.so', additionalLibPaths)
        if l is None:
            l = checkLibFile__(libPrefix+'blas*.a', additionalLibPaths)
        if l is not None: libnames.append(libname); foundBlas = l
        libname = 'cblas'; foundCBlas = None
        l = checkLibFile__(libPrefix+'cblas*.so', additionalLibPaths)
        if l is None:
            l = checkLibFile__(libPrefix+'blas*.a', additionalLibPaths)
        if l is not None: libnames.append(libname); foundCBlas = l
        l = None
        if foundBlas is not None: l = foundBlas
        elif foundCBlas is not None: l = foundCBlas

    # force of mkl
    if libPrefix == 'libmkl_': libnames = ['mkl']; l = ''

    # Chemin des includes
    i = checkIncFile__(includePrefix+'cblas.h', additionalIncludePaths)
    if i is not None and l is not None:
        print('Info: Blas detected at %s.'%l)
        return (True, i, l, compOpt, libnames)
    elif l is None and i is not None:
        print('Info: libblas was not found on your system. No Blas support.')
        return (False, i, l, compOpt, libnames)
    elif l is not None and i is None:
        print('Info: cblas.h was not found on your system. No Blas support.')
        return (False, i, l, compOpt, libnames)
    else:
        print('Info: libblas and cblas.h was not found on your system. No Blas support.')
        return (False, i, l, compOpt, libnames)

#=============================================================================
# Check for LAPACK
# additionalPaths: chemins d'installation non standards : ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie,
# option de compile, nom de la librarie)
#=============================================================================
def checkLapack(additionalLibPaths=[], additionalIncludePaths=[]):
    try: from KCore.config import Cppcompiler
    except: from config import Cppcompiler
    if Cppcompiler == 'icc' or Cppcompiler == 'icpc' or Cppcompiler == 'icx': # intel - cherche dans MKL
        libPrefix = 'libmkl_'; includePrefix = 'mkl_'; compOpt = '-mkl'
    else: # cherche std
        libPrefix = 'lib'; includePrefix = ''; compOpt = ''
    libnames = []

    # Check for openblas, it contains lapack generally
    libname = 'openblas'
    l = checkLibFile__(libPrefix+'openblas*.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__(libPrefix+'openblas*.a', additionalLibPaths)
    if l is not None: libnames.append(libname) # no further search

    else: # Check lapack and lapacke
        libname = 'lapack'; foundLapack = None
        if l is None:
            l = checkLibFile__(libPrefix+'lapack*.so', additionalLibPaths)
        if l is None:
            l = checkLibFile__(libPrefix+'lapack*.a', additionalLibPaths)
        if l is not None: libnames.append(libname); foundLapack = l; l = None
        libname = 'lapacke'; foundLapacke = None
        if l is None:
            l = checkLibFile__(libPrefix+'lapacke*.so', additionalLibPaths)
        if l is None:
            l = checkLibFile__(libPrefix+'lapacke*.a', additionalLibPaths)
        if l is not None: libnames.append(libname); foundLapacke = l

        l = None
        if foundLapack is not None: l = foundLapack
        elif foundLapacke is not None: l = foundLapacke

    i = checkIncFile__(includePrefix+'lapack.h', additionalIncludePaths)
    if i is None:
        i = checkIncFile__(includePrefix+'lapacke.h', additionalIncludePaths)
    if i is not None and l is not None:
        print('Info: Lapack detected at %s.'%l)
        return (True, i, l, compOpt, libnames)
    elif l is None and i is not None:
        print('Info: liblapack was not found on your system. No Lapack support.')
        return (False, i, l, compOpt, libnames)
    elif l is not None and i is None:
        print('Info: lapack.h or lapacke.h was not found on your system. No Lapack support.')
        return (False, i, l, compOpt, libnames)
    else:
        print('Info: liblapack and lapack.h or lapacke.h was not found on your system. No Lapack support.')
        return (False, i, l, compOpt, libnames)

#=============================================================================
# Check for Cython
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkCython(additionalLibPaths=[], additionalIncludePaths=[]):
    try:
        import Cython.Compiler.Main as cython_compiler
        try:
            import Cython
            cythonVersion = Cython.__version__
        except: cythonVersion = True
    except: cythonVersion = False
    if cythonVersion != False:
        print('Info: found Cython version '+cythonVersion)
    return cythonVersion

def cythonize(src):
    import Cython.Compiler.Main as cython_compiler
    sys.stderr.write("cythonize: %r\n" % (src,))
    cython_compiler.compile([src], cplus=True)

#=============================================================================
# Check fortran libs
# additionalLibs: noms des libraries utilises par fortran si non conventionnels
# additionalLibPaths: chemins d'installation des libraries non standards: ['/home/toto',...]
# Retourne (True, [librairies utiles pour le fortran], [paths des librairies])
#=============================================================================
def checkFortranLibs(additionalLibs=[], additionalLibPaths=[],
                     f77compiler=None, useOMP=None):
    if f77compiler is None:
        try: from KCore.config import f77compiler
        except:
            try: from config import f77compiler
            except: f77compiler = 'gfortran'
    if useOMP is None:
        try: from KCore.config import useOMP
        except:
            try: from config import useOMP
            except: useOMP = True
    ret = True; libs = []; paths = []

    # librairies speciales (forcees sans check)
    libs += additionalLibs
    paths += additionalLibPaths

    # gfortran (gfortran, gomp)
    if f77compiler == 'gfortran':
        l = checkLibFile__('libgfortran.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libgfortran.a', additionalLibPaths)

        if l is not None:
            libs += ['gfortran']; paths += [l]

        if useOMP:
            l = checkLibFile__('libgomp.so*', additionalLibPaths)
            if l is None:
                l = checkLibFile__('libgomp.a', additionalLibPaths)
            if l is not None:
                libs += ['gomp']; paths += [l]
            else: ret = False

    # ifort (ifcore, svml, irc, guide, iomp5)
    elif f77compiler == 'ifort':
        l = checkLibFile__('libifcore.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libifcore.a', additionalLibPaths)

        if l is not None:
            libs += ['ifcore']; paths += [l]

        l = checkLibFile__('libsvml.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libsvml.a', additionalLibPaths)
        if l is not None:
            libs += ['svml']; paths += [l]

        l = checkLibFile__('libirc.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libirc.a', additionalLibPaths)
        if l is not None:
            libs += ['irc']; paths += [l]

        if useOMP:
            # check guide
            l = checkLibFile__('libguide.so*', additionalLibPaths)
            if l is None:
                l = checkLibFile__('libguide.a', additionalLibPaths)
            if l is not None:
                libs += ['guide']; paths += [l]
            # check iomp5
            if l is None:
                l = checkLibFile__('libiomp5.so*', additionalLibPaths)
                if l is None:
                    l = checkLibFile__('libiomp5.a', additionalLibPaths)
                if l is not None:
                    libs += ['iomp5']; paths += [l]
                else: ret = False

    # ifx (ifcore, svml, irc, guide, iomp5)
    elif f77compiler == 'ifx':
        l = checkLibFile__('libifcore.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libifcore.a', additionalLibPaths)

        if l is not None:
            libs += ['ifcore']; paths += [l]

        l = checkLibFile__('libsvml.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libsvml.a', additionalLibPaths)
        if l is not None:
            libs += ['svml']; paths += [l]

        l = checkLibFile__('libirc.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libirc.a', additionalLibPaths)
        if l is not None:
            libs += ['irc']; paths += [l]

        if useOMP:
            # check guide
            l = checkLibFile__('libguide.so*', additionalLibPaths)
            if l is None:
                l = checkLibFile__('libguide.a', additionalLibPaths)
            if l is not None:
                libs += ['guide']; paths += [l]
            # check iomp5
            if l is None:
                l = checkLibFile__('libiomp5.so*', additionalLibPaths)
                if l is None:
                    l = checkLibFile__('libiomp5.a', additionalLibPaths)
                if l is not None:
                    libs += ['iomp5']; paths += [l]
                else: ret = False

    # pgfortran
    elif f77compiler == 'pgfortran':
        l = checkLibFile__('libnvf.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libnvf.a', additionalLibPaths)
        if l is not None:
            libs += ['nvf', 'rt']; paths += [l]

        if useOMP:
            l = checkLibFile__('libnvomp.so*', additionalLibPaths)
            if l is None:
                l = checkLibFile__('libnvomp.a', additionalLibPaths)
            if l is not None:
                libs += ['nvomp']; paths += [l]
            else: ret = False

    # nvfortran
    elif f77compiler == 'nvfortran':
        l = checkLibFile__('libnvf.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libnvf.a', additionalLibPaths)
        if l is not None:
            libs += ['nvf', 'rt']; paths += [l]

        if useOMP:
            l = checkLibFile__('libnvomp.so*', additionalLibPaths)
            if l is None:
                l = checkLibFile__('libnvomp.a', additionalLibPaths)
            if l is not None:
                libs += ['nvomp']; paths += [l]
            else: ret = False

    # crayftn
    elif f77compiler == 'crayftn':
        l = checkLibFile__('libf.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libf.a', additionalLibPaths)
        if l is not None:
            libs += ['f', 'rt']; paths += [l]

        #if useOMP:
        #    l = checkLibFile__('libomp.so*', additionalLibPaths)
        #    if l is None:
        #        l = checkLibFile__('libomp.a', additionalLibPaths)
        #    if l is not None:
        #        libs += ['omp']; paths += [l]
        #    else: ret = False

    return (ret, libs, paths)

#=============================================================================
# Check Cpp libs
# additionalLibs: si les libs requises ont des noms non conventionnels
# additionalLibPaths: si chemins des libraries necessaires non standards : ['/home/toto',...]
# Retourne (True, [librairies utiles a cpp], [paths des librairies])
#=============================================================================
def checkCppLibs(additionalLibs=[], additionalLibPaths=[], Cppcompiler=None,
                 useOMP=None):

    if Cppcompiler is None:
        try: from KCore.config import Cppcompiler
        except:
            try: from config import Cppcompiler
            except: Cppcompiler = 'gcc'
    if useOMP is None:
        try: from KCore.config import useOMP
        except:
            try: from config import useOMP
            except: useOMP = True

    ret = True; libs = []; paths = []

    # librairies additionales (forcees sans check)
    libs += additionalLibs
    paths += additionalLibPaths

    # gcc (stdc++, gomp)
    if Cppcompiler.find('gcc') == 0 or Cppcompiler.find('g++') == 0:
        os.environ['CC'] = 'gcc'
        os.environ['CXX'] = 'g++'
        cflags = sysconfig.get_config_var('CFLAGS')
        sysconfig._config_vars['CFLAGS'] = '' # kill setup flags for CC
        sysconfig._config_vars['LDFLAGS'] = '' # kill setup flags for LD

        l = checkLibFile__('libstdc++.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libstdc++.a', additionalLibPaths)
        if l is not None:
            libs += ['stdc++']; paths += [l]

        if DEBUG:
            l = checkLibFile__('libasan.so*', additionalLibPaths)
            if l is None:
                l = checkLibFile__('libasan.a', additionalLibPaths)
            if l is not None: libs += ["asan"]
            #l = checkLibFile__('libtsan.so*', additionalLibPaths)
            #if l is None:
            #    l = checkLibFile__('libtsan.a', additionalLibPaths)
            #if l is not None: libs += ["tsan"]

        if useOMP:
            l = checkLibFile__('libgomp.so*', additionalLibPaths)
            if l is None:
                l = checkLibFile__('libgomp.a', additionalLibPaths)
            if l is not None:
                libs += ['gomp']; paths += [l]
            else: ret = False

    # icc (stdc++, guide ou iomp5)
    if Cppcompiler.find('icc') == 0 or Cppcompiler.find('icpc') == 0:
        os.environ['CC'] = 'icc' # forced to overide setup.cfg
        os.environ['CXX'] = 'icpc'
        cflags = sysconfig.get_config_var('CFLAGS')
        sysconfig._config_vars['CFLAGS'] = '' # kill setup flags for CC
        sysconfig._config_vars['LDFLAGS'] = '' # kill setup flags for LD
        l = checkLibFile__('libstdc++.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libstdc++.a', additionalLibPaths)
        if l is not None:
            libs += ['stdc++']; paths += [l]
        #if DEBUG:
        #    l = checkLibFile__('libchkpwrap.a', additionalLibPaths)
        #    if l is not None:
        #        libs += ['chkpwrap', 'chkp']; paths += [l]

        if useOMP:
            l = checkLibFile__('libguide.so*', additionalLibPaths)
            if l is None:
                l = checkLibFile__('libguide.a', additionalLibPaths)
            if l is not None:
                libs += ['guide']; paths += [l]
            if l is None:
                l = checkLibFile__('libiomp5.so*', additionalLibPaths)
                if l is None:
                    l = checkLibFile__('libiomp5.a', additionalLibPaths)
                if l is not None:
                    libs += ['iomp5']; paths += [l]
                else: ret = False

    # icx
    if Cppcompiler == 'icx' or Cppcompiler == 'icpx':
        os.environ['CC'] = 'icx' # forced to overide setup.cfg
        os.environ['CXX'] = 'icpx'
        os.environ['LDSHARED'] = 'ifx'
        cflags = sysconfig.get_config_var('CFLAGS')
        sysconfig._config_vars['CFLAGS'] = '' # kill setup flags for CC
        sysconfig._config_vars['LDFLAGS'] = '' # kill setup flags for LD
        l = checkLibFile__('libstdc++.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libstdc++.a', additionalLibPaths)
        if l is not None:
            libs += ['stdc++']; paths += [l]

        if useOMP:
            l = checkLibFile__('libguide.so*', additionalLibPaths)
            if l is None:
                l = checkLibFile__('libguide.a', additionalLibPaths)
            if l is not None:
                libs += ['guide']; paths += [l]
            if l is None:
                l = checkLibFile__('libiomp5.so*', additionalLibPaths)
                if l is None:
                    l = checkLibFile__('libiomp5.a', additionalLibPaths)
                if l is not None:
                    libs += ['iomp5']; paths += [l]
                else: ret = False

    # pgcc
    if Cppcompiler == 'pgcc' or Cppcompiler == 'pgc++':
        os.environ['CC'] = 'pgc++' # forced to overide setup.cfg
        os.environ['CXX'] = 'pgc++'
        os.environ['LDSHARED'] = 'pgfortran'
        cflags = sysconfig.get_config_var('CFLAGS')
        sysconfig._config_vars['CFLAGS'] = '' # kill setup flags for CC
        sysconfig._config_vars['LDFLAGS'] = '' # kill setup flags for LD

        l = checkLibFile__('libnvc.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libnvc.a', additionalLibPaths)
        if l is not None:
            libs += ['nvc']; paths += [l]

        if useOMP:
            l = checkLibFile__('libnvomp.so*', additionalLibPaths)
            if l is None:
                l = checkLibFile__('libnvomp.a', additionalLibPaths)
            if l is not None:
                libs += ['nvomp']; paths += [l]
            else: ret = False

    # nvc
    if Cppcompiler == 'nvc' or Cppcompiler == 'nvc++':
        os.environ['CC'] = 'nvc++' # forced to overide setup.cfg
        os.environ['CXX'] = 'nvc++'
        os.environ['LDSHARED'] = 'nvfortran'
        cflags = sysconfig.get_config_var('CFLAGS')
        sysconfig._config_vars['CFLAGS'] = '' # kill setup flags for CC
        sysconfig._config_vars['LDFLAGS'] = '' # kill setup flags for LD

        l = checkLibFile__('libnvc.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libnvc.a', additionalLibPaths)
        if l is not None:
            libs += ['nvc']; paths += [l]

        if useOMP:
            l = checkLibFile__('libnvomp.so*', additionalLibPaths)
            if l is None:
                l = checkLibFile__('libnvomp.a', additionalLibPaths)
            if l is not None:
                libs += ['nvomp']; paths += [l]
            else: ret = False

    # craycc
    if Cppcompiler == 'craycc' or Cppcompiler == 'craycxx':
        os.environ['CC'] = 'craycc' # forced to overide setup.cfg
        os.environ['CXX'] = 'craycxx'
        os.environ['LDSHARED'] = 'crayftn'
        cflags = sysconfig.get_config_var('CFLAGS')
        sysconfig._config_vars['CFLAGS'] = '' # kill setup flags for CC
        sysconfig._config_vars['LDFLAGS'] = '' # kill setup flags for LD

        #l = checkLibFile__('libcraymp.so*', additionalLibPaths)
        #if l is None:
        #    l = checkLibFile__('libcraymp.a', additionalLibPaths)
        #if l is not None:
        #    libs += ['craymp']; paths += [l]
        l = checkLibFile__('libsci_cray.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libsci_cray.a', additionalLibPaths)
        if l is not None:
            libs += ['sci_cray']; paths += [l]
        #if useOMP:
        #    l = checkLibFile__('libomp.so*', additionalLibPaths)
        #    if l is None:
        #        l = checkLibFile__('libomp.a', additionalLibPaths)
        #    if l is not None:
        #        libs += ['omp']; paths += [l]
        #    else: ret = False
        # craycc

    if Cppcompiler == 'cc':
        os.environ['CC'] = 'cc' # forced to overide setup.cfg
        os.environ['CXX'] = 'cc'
        os.environ['LDSHARED'] = 'ftn'
        cflags = sysconfig.get_config_var('CFLAGS')
        sysconfig._config_vars['CFLAGS'] = '' # kill setup flags for CC
        sysconfig._config_vars['LDFLAGS'] = '' # kill setup flags for LD
        l = checkLibFile__('libsci_cray.so*', additionalLibPaths)
        if l is None:
            l = checkLibFile__('libsci_cray.a', additionalLibPaths)
        if l is not None:
            libs += ['sci_cray']; paths += [l]

    return (ret, libs, paths)

#==============================================================================
# Check if file exists LD_LIBRARY_PATH dirs and specified lib dirs
# additionalPaths est en premier pour pouvoir forcer une librairie par config.
#==============================================================================
def checkLibFile__(file, additionalLibPaths):
    p = []
    for i in additionalLibPaths: p.append(i)
    mySystem = getSystem()
    env = os.environ

    if mySystem[0] == 'Windows':
        p1 = env.get('PATH', None)
        if p1 is not None: p += p1.split(';')
    elif mySystem[0] == 'mingw':
        p1 = env.get('PATH', None)
        if p1 is not None:
            p += p1.split(';')
    else: # unix
        p1 = env.get('LD_LIBRARY_PATH', None)
        if p1 is not None: p += p1.split(':')
        p1 = env.get('PATH', None)
        if p1 is not None: p += p1.split(';')
    p1 = env.get('CMAKE_PREFIX_PATH', None)
    if p1 is not None: p += [path+'/lib' for path in p1.split(':')]
    #p += ['/usr/local/lib', '/opt/lib', '/usr/lib', '/opt/local/lib']
    for i in p:
        a = glob.glob(i+'64/'+file)
        if a != []: return i+'64'
        a = glob.glob(i+'/'+file)
        if a != []: return i
    return None

#==============================================================================
# Check if file exists in path dirs and specified include dirs
#==============================================================================
def checkIncFile__(file, additionalIncludePaths):
    p = []
    for i in additionalIncludePaths: p.append(i)
    mySystem = getSystem()
    env = os.environ
    pp = []
    if mySystem[0] == 'Windows':
        p1 = env.get('PATH', None)
        if p1 is not None: pp += p1.split(';')
    elif mySystem[0] == 'mingw':
        p1 = env.get('PATH', None)
        if p1 is not None:
            pp += p1.split(';')
    else: # unix
        p1 = env.get('LD_LIBRARY_PATH', None)
        if p1 is not None: pp += p1.split(':')
        p1 = env.get('PATH', None)
        if p1 is not None: pp += p1.split(':')
    #p += ['/usr/local/include', '/opt/include', '/usr/include', '/opt/local/include']
    p1 = env.get('CMAKE_PREFIX_PATH', None)
    if p1 is not None: pp += [path+'/lib' for path in p1.split(':')]
    for i, v in enumerate(pp):
        s = v.split('/'); ls = len(s)
        if ls > 0 and s[-1] == 'lib': s[-1] = 'include'
        if ls > 1 and s[-2] == 'lib': s[-2] = 'include'; s[-1] = ''
        if ls > 0 and s[-1] == 'lib64': s[-1] = 'include'
        if ls > 1 and s[-2] == 'lib64': s[-2] = 'include'; s[-1] = ''
        out = ''
        for j in s: out += j+'/'
        pp[i] = out[:-1]
    p += pp
    for i in p:
        pf = i+'/'+file
        a = os.access(pf, os.F_OK)
        if a: return i
    return None

#==============================================================================
# Ecrit les infos de build dans un fichier buildInfo.py
#==============================================================================
def writeBuildInfo():
    p = open("buildInfo.py", 'w')
    if p is None:
        raise SystemError("Error: can not open file buildInfo.py for writing.")
    try: import KCore.config as config
    except: import config

    dict = {}
    # Date
    import time
    execTime = time.strftime('%d/%m/%y %Hh%M', time.localtime())
    dict['date'] = execTime

    # Check python
    (pythonVersion, pythonIncDir, pythonLibDir, pythonLibs) = checkPython()
    if pythonVersion != False: dict['python'] = pythonVersion
    else: dict['numpy'] = "None"

    # Check numpy
    (numpyVersion, numpyIncDir, numpyLibDir) = checkNumpy()
    if numpyVersion != False: dict['numpy'] = numpyVersion
    else: dict['numpy'] = "None"

    # Check png
    #(png, pngIncDir, pngLib) = checkPng(config.additionalLibPaths,
    #                                    config.additionalIncludePaths)
    #if png: dict['png'] = pngLib
    #else: dict['png'] = "None"

    # Check ffmpeg
    (mpeg, mpegIncDir, mpegLib) = checkMpeg(config.additionalLibPaths,
                                            config.additionalIncludePaths)
    if mpeg: dict['mpeg'] = mpegLib
    else: dict['mpeg'] = "None"

    # Check hdf5
    (hdf, hdfIncDir, hdfLib, hdflibnames) = checkHdf(config.additionalLibPaths,
                                                     config.additionalIncludePaths)
    if hdf: dict['hdf'] = hdfLib
    else: dict['hdf'] = "None"

    # Check netcdf
    (netcdf, netcdfIncDir, netcdfLib, netcdflibnames) = checkNetcdf(config.additionalLibPaths,
                                                                    config.additionalIncludePaths)
    if netcdf: dict['netcdf'] = netcdfLib
    else: dict['netcdf'] = "None"

    # Check mpi
    (mpi, mpiIncDir, mpiLib, mpiLibs) = checkMpi(config.additionalLibPaths,
                                                 config.additionalIncludePaths)
    if mpi: dict['mpi'] = mpiLib
    else: dict['mpi'] = "None"

    # Check cuda
    (cuda, cudaIndDir, cudaLib, cudalibNames, cudaexec) = checkCuda(config.additionalLibPaths,
                                                                    config.additionalIncludePaths)
    if cuda: dict['cuda'] = cudaLib
    else: dict['cuda'] = "None"

    # Write dictionnary
    p.write("# This file is generated by Cassiopee installer.\n")
    p.write("buildDict = "+str(dict))
    p.close()

#==============================================================================
# Ecrit la base d'installation (ancien config.py) dans le fichier
# installBase.py
# IN: dict: dictionnaire d'install
#==============================================================================
def writeInstallBase(dict):
    p = open("installBase.py", 'w')
    if p is None:
        raise SystemError("Error: can not open file installBase.py for writing.")

    # Write doc
    p.write("# This is the dictionary keeping track of installation.\n# The key is the machine name or ELSAPROD name. For each key a list is stored.\n# [description, \n# f77compiler, libfortdir, libfort, f90compiler, libf90dir, libf90, \n# Cppcompiler, libCpp, useOMP, \n# pngPath, mpegPath, adfPath, hdfPath].\n# Path are list of strings. useOMP, static are booleans. \n# Others are strings.\n")

    # Write dictionary
    #p.write("installDict = "+str(dict))

    # Pretty print dict
    p.write("installDict = {\n")
    kc = 0
    for k in dict:
        p.write("###############################################################################\n")
        if isinstance(k, str): kstr = "\'%s\'"%k
        else: kstr = str(k)
        p.write("%s: [ "%kstr)
        list = dict[k]
        lc = 0
        for l in list:
            lc += 1
            if isinstance(l, str): lstr = "\'%s\'"%l
            else: lstr = str(l)
            if lc == 1:  p.write("%s,\n"%lstr)
            if lc == 2:  p.write("%s, # f77compiler\n"%lstr)
            elif lc == 3: p.write("%s, # f90compiler\n"%lstr)
            elif lc == 4: p.write("%s, # Cppcompiler\n"%lstr)
            elif lc == 5: p.write("%s, # CppAdditionalOptions\n"%lstr)
            elif lc == 6: p.write("%s, # f77AdditionalOptions\n"%lstr)
            elif lc == 7: p.write("%s, # useOMP\n"%lstr)
            elif lc == 8: p.write("%s, # static\n"%lstr)
            elif lc == 9: p.write("%s, # additionalIncludePaths\n"%lstr)
            elif lc == 10: p.write("%s, # additionalLibs\n"%lstr)
            elif lc == 11: p.write("%s, # additionalLibPaths\n"%lstr)
            elif lc == 12: p.write("%s, # useCuda\n"%lstr)
            elif lc == 13: p.write("%s  # NvccAdditionalOptions\n"%lstr)
        kc += 1
        if kc == len(dict): p.write("]\n")
        else: p.write("], \n")
    p.write("}\n")
    p.close()

#==============================================================================
# Sur certains unix, le chemin d'installation contient un lib64
# Cree le symlink pour que lib et lib64 soit equivalent
#==============================================================================
def symLinks():
    system = getSystem()[0]
    bits = getSystem()[1]
    if bits == '64':
        try: import KCore.installPath as K
        except: import installPath as K
        libPath1 = K.libPath
        spls = libPath1.rsplit('/',1)
        if spls[1] == 'lib': libPath2 = spls[0]+'/lib64'
        elif spls[1] == 'lib64': libPath2 = spls[0]+'/lib'
        else: return
        ex1 = os.path.exists(libPath1)
        ex2 = os.path.exists(libPath2)
        lex1 = os.path.lexists(libPath1)
        lex2 = os.path.lexists(libPath2)
        if not ex1 and lex1: # broken link 1
            os.remove(libPath1); ex1 = False
        if not ex2 and lex2: # broken link 2
            os.remove(libPath2); ex2 = False
        if ex1 and not ex2:
            try: # may fail under windows disk system
                os.symlink(libPath1, libPath2)
            except: pass
        elif not ex1 and ex2:
            try: # may fail under windows disk system
                os.symlink(libPath2, libPath1)
            except: pass

#==============================================================================
# Cree les extensions d'un Module
# Cree une extension C/Api : Module.module (a partir de module.cpp)
# Des extensions Cython si necessaire: Module.f (a partir de la liste des pyx
# IN: module: module name ('Converter')
# IN: srcs: le module des sources
#==============================================================================
def createExtensions(module, srcs, includeDirs, libraryDirs, libraries,
                     extraCompileArgs=[], extraLinkArgs=[]):
    listExtensions = []
    minor = module.lower()
    # C/Api module
    Extension(module+'.'+minor,
              sources=[module+'/'+minor+'.cpp'],
              include_dirs=[module]+includeDirs,
              library_dirs=libraryDirs,
              libraries=libraries,
              extra_compile_args=extraCompileArgs,
              extra_link_args=extraLinkArgs)
    # Cython extensions
    try: pyx_srcs = srcs.pyx_srcs
    except: pyx_srcs = []
    for i in pyx_srcs:
        f = i.replace('.pyx', '')
        Extension(module+'.'+f,
                  sources=[module+'/'+f+'.cpp'],
                  include_dirs=[module]+includeDirs,
                  library_dirs=libraryDirs,
                  libraries=libraries,
                  extra_compile_args=extraCompileArgs,
                  extra_link_args=extraLinkArgs)
    return listExtensions

#==============================================================================
# Retourne un dictionnaire des variables de l'environnement que scons
# doit voir
#==============================================================================
def getEnvForScons():
    #return dict(os.environ) # for adastra
    return {'PATH': getenv('PATH'),
            'LD_LIBRARY_PATH': getenv('LD_LIBRARY_PATH'),
            'LM_LICENSE_FILE': getenv('LM_LICENSE_FILE'),
            'INTEL_LICENSE_FILE': getenv('INTEL_LICENSE_FILE'),
            'TMP':getenv('TMP'),
            'TEMP':getenv('TEMP')}

#==============================================================================
# Fortran builder
#==============================================================================
# Ajoute le fortran builder a env
# IN: dirs: include paths
def createFortranBuilder(env, dirs=[], additionalPPArgs='', additionalFortranArgs=[]):
    import SCons
    from SCons.Builder import Builder
    # Pre-processing
    path = ''
    for i in dirs: path += '"%s" -I'%i
    if path != '': path = path[:-3]
    PP = getPP()
    if additionalPPArgs != '': PP = PP[0:-2]+' '+additionalPPArgs+' -I'
    bld = Builder(action=PP+'%s $SOURCES $TARGETS'%path, suffix='.f',
                  src_suffix='.for')
    env.Append(BUILDERS={'FPROC': bld})
    # Fortran compiler
    fortran_builder = Builder(action='$FORTRANCOM',
                              suffix='.o', src_suffix='.f')
    env.Append(BUILDERS={'Fortran': fortran_builder})
    env.Replace(FORTRANCOM='$FORTRAN $FORTRANFLAGS -c -o $TARGET $SOURCE')
    env.Replace(FORTRANSUFFIXES=['.f', '.F', '.f90', '.F90'])
    #env.Replace(FORTRANFLAGS=getForArgs())
    #env.Replace(F90FLAGS=getForArgs())
    #env.Replace(F95FLAGS=getForArgs())
    #env.Replace(SHF90=f90compiler)
    #env.Replace(SHF95=f90compiler)
    pref = getFortranModDirPrefix()
    env.Replace(FORTRANMODDIRPREFIX=pref) # doesnt work
    env.Replace(FORTRANMODDIR='MODS') # doesnt work
    args = getForArgs()
    args += additionalFortranArgs
    if pref != '': args += [pref,'build'] # trick
    env.Replace(FORTRANFLAGS=args)
    return env

# Scan des fichiers (fortran par defaut) et retourne un dict de dependances de
# premier niveau
def findImmediateDeps(parentFolder, searchFolder,
                      depPattern=r'^#include\s*["\'](.+?)["\']',
                      fileExtensions=['.f', '.f90', '.for']):
    import re
    searchFolder = os.path.join(parentFolder, searchFolder)
    # Fortran dependency dict mapping a source file to a list of includes
    deps = {"parentFolder": parentFolder}
    # Regexpr to find include statements
    regInclude = re.compile(depPattern, re.IGNORECASE)
    # Find all fortran files in root directory
    for root, _, files in os.walk(searchFolder):
        for infile in files:
            fileExt = os.path.splitext(infile)[1]
            if fileExt.lower() not in fileExtensions: continue
            filePath = os.path.join(root, infile)
            filePathRel = filePath.replace(parentFolder, "")
            if not filePathRel[0].isalpha(): filePathRel = filePathRel[1:]
            deps[filePathRel] = []
            # Search for include statements and add dependence to dict
            with open(filePath, 'r') as f:
                for line in f:
                    incFound = regInclude.search(line)
                    if incFound is not None:
                        includeFile = incFound.group(1)
                        includePath = os.path.join(parentFolder, includeFile)
                        if os.path.exists(includePath):
                            deps[filePathRel].append(includePath)
    return deps

# Find all dependencies (include) of a file recursively
def findAllDeps(filename, deps={}, cache=None):
    if cache is None: cache = {}
    # Use memoization for dependencies that have already been established
    if filename in cache: return cache[filename]
    includes = deps.get(filename, [])
    # Store all dependencies in a set for unicity
    allIncludes = set()
    for inc in includes:
        allIncludes.add(inc)
        relInc = inc.replace(deps["parentFolder"], "")
        if not relInc[0].isalpha(): relInc = relInc[1:]
        nestedDeps = findAllDeps(relInc, deps, cache)
        allIncludes.update(nestedDeps)
    cache[filename] = allIncludes
    return sorted(allIncludes) # sorting is important, recompiles all otherwise

# Ajoute les dependances au Fortran builder
def envFortranWithDeps(env, filename, deps={}):
    if filename.endswith('90'): target = filename
    else: target = env.FPROC(target=filename)
    includes = findAllDeps(filename, deps=deps)
    for inc in includes: env.Depends(target, env.File(inc))
    return env.Fortran(target=target)

# Cree les noms des fichiers
def createFortranFiles(env, srcs, deps={}):
    for_srcs = []
    try:
        for_srcs.extend(srcs.for_srcs[:])
    except: pass
    try:
        for_srcs.extend(srcs.f90_srcs[:])
    except: pass
    ppf = []
    for f in for_srcs:
        ofile = envFortranWithDeps(env=env, filename=f, deps=deps)
        ppf.append(ofile[0])
    return ppf

# Decoupe une liste de fichiers object en morceaux de taille egale a chunkSize
def chunkObjectFiles(ppf, chunkSize=100):
    chunked_ppf = []
    for i in range(0, len(ppf), chunkSize):
        chunked_ppf.append(ppf[i:i+chunkSize])
    return chunked_ppf

# Scan les .f pour faire les dependences (include)
def fortranScan(node, env, path, arg=None):
    import re
    # scan file to extract all possible includes.
    contents = node.get_text_contents()
    reg = re.compile(r'include\s+(\S+)$', re.M)
    names = reg.findall(contents)
    names = [name.strip() for name in names]
    names = [name.replace('"', '') for name in names]

    # remove duplications
    names = set(names)
    # remove false dep (+in KCore?)
    names = [n for n in names if not 'omp_lib.h' in n]
    return env.File(list(names))

# Cree le scanner Fortran dans env
def createFortranScanner(env):
    import SCons
    fortranscanner = SCons.Scanner.Scanner(function=fortranScan, skeys=['.for'],
                                           recursive=True)
    env.Append(SCANNERS=fortranscanner)
    return env

# Cree le scanner Cuda dans env
def createCudaScanner(env):
    import SCons
    CudaScanner = SCons.Scanner.C.CScanner()
    SCons.Tool.SourceFileScanner.add_scanner(['.cu'], CudaScanner)
    return env

def addCommonNvccVariables(env):
    """
    Add underlying common "NVIDIA CUDA compiler" variables that
    are used by multiple builders.
    """

    # "NVCC common command line"
    if not env.has_key('_NVCCCOMCOM'):
        # nvcc needs '-I' prepended before each include path, regardless of platform
        env['_NVCCWRAPCPPPATH'] = '${_concat("-I", CPPPATH, "", __env__)}'
        # prepend -Xcompiler before each flag
        # assemble the common command line
        env['_NVCCCOMCOM'] = '$_NVCCWRAPCPPPATH'
    return env

def createCudaBuilders(env, dirs=[]):
    import SCons
    # create a builder that makes PTX files from .cu files
    (ok, incCuda, libCude, libNameCuda, binCuda) = checkCuda()
    opts = getCppArgs()
    addCommonNvccVariables(env)
    path = ''
    for i in dirs: path += '-I"%s" '%i
    action_cuda = binCuda + ' -ptx '
    #for o in opts :
    #    action_cuda += o + " "
    action_cuda += path + '$NVCCFLAGS $_NVCCCOMCOM $SOURCES -o $TARGET'
    ptx_builder = SCons.Builder.Builder(action=action_cuda,
                                        emitter={},
                                        suffix='.ptx',
                                        src_suffix=['.cu'])
    env['BUILDERS']['PTXFile'] = ptx_builder

    # create builders that make static & shared objects from .cu files
    static_obj, shared_obj = SCons.Tool.createObjBuilders(env)

    NVCCCOM   = '"' + binCuda + '"' + ' -o $TARGET -c '
    SHNVCCCOM = '"' + binCuda + '"' + ' -shared -o $TARGET -c '
    #for o in opts :
    #    NVCCCOM   += o + " "
    #    SHNVCCCOM += o + " "

    NVCCCOM   += '$NVCCFLAGS $_NVCCWRAPCFLAGS $NVCCWRAPCCFLAGS $_NVCCCOMCOM $SOURCES'
    SHNVCCCOM += '$NVCCFLAGS $_NVCCWRAPSHCFLAGS $_NVCCWRAPSHCCFLAGS $_NVCCCOMCOM $SOURCES'
    # Add this suffix to the list of things buildable by Object
    static_obj.add_action('.cu', NVCCCOM)
    shared_obj.add_action('.cu', SHNVCCCOM)
    static_obj.add_emitter('.cu', SCons.Defaults.StaticObjectEmitter)
    shared_obj.add_emitter('.cu', SCons.Defaults.SharedObjectEmitter)

    return env

#==============================================================================
# Builder Cython
#==============================================================================
def cythonSuffixEmitter(env, source):
    return "$CYTHONCFILESUFFIX"

# Cree le builder Cython dans env
def createCythonBuilder(env):
    import SCons
    from SCons.Builder import Builder
    from SCons.Action import Action
    incs = " -a --cplus "
    for i in env['CPPPATH']: incs += ' -I"%s" '%i

    cypath = ''
    if "CYTHONCOMPATH" in env and env["CYTHONCOMPATH"] != "":
        SYSPATH = ""
        if "PYTHONPATH" in env: SYSPATH=env['PYTHONPATH']
        cypath = env["CYTHONCOMPATH"]
        if cypath != "": cypath += "/"
        #pypath = "PYTHONPATH=%s:%s:%s "%(RESOURCELIBPATH,CYTHONLIBPATH,SYSPATH)
        #cypath = pypath+cypath
    env["CYTHONCOM"] = cypath+"cython"+incs+" -o $TARGET $SOURCE"
    env["CYTHONCFILESUFFIX"] = ".cpp"

    c_file, cxx_file = SCons.Tool.createCFileBuilders(env)
    c_file.suffix['.pyx'] = cythonSuffixEmitter
    c_file.add_action('.pyx', Action("$CYTHONCOM"))

    cython = Builder(
        action=Action("$CYTHONCOM"),
        emitter={},
        suffix='.cpp', src_suffix='.pyx',
        single_source=1)
    env.Append(BUILDERS={'GenerateCython': cython})
    #env = createCythonScanner(env)
    return env

# Scan les .pyx pour faire les dependences (.pxd, .pxi)
def cythonScan(node, env, path, arg=None):
    import itertools
    import re
    # scan file to extract all possible cimports.
    contents = node.get_text_contents()
    names = [reo.findall(contents) for reo in [
        re.compile(r'^\s*from\s+(.+?)\s+cimport\s.*$', re.M),
        re.compile(r'^\s*cimport\s+(.+?)$', re.M),
    ]]
    names = itertools.chain(*names)
    # keep characters before " as ".
    names = [name.split(' as ')[0] for name in names]
    # split each cimport.
    names = itertools.chain(*[name.split(',') for name in names])
    names = [name.strip() for name in names]
    # remove duplications
    names = set(names)
    # prepend with the directory of the original pyx file.
    prefix = os.path.dirname(env.GetBuildPath(node))
    names = [os.path.join(prefix, '%s.pxd'%name) for name in names]
    # only take local pxd file and ignore anything unfound.
    names = [name for name in names if os.path.exists(name)]
    return [env.File(name) for name in names]

# Cree le scanner Cython dans env
def createCythonScanner(env):
    import SCons
    pyxscanner = SCons.Scanner.Scanner(function=cythonScan, skeys=['.pyx'], name='PYX')
    env.Append(SCANNERS=pyxscanner)
    return env

# Cree les targets des fichiers CPP
def createCythonFiles(env, srcs):
    cythonCpp = env.GenerateCython(srcs.pyx_srcs)
    deps = []
    for i in cythonCpp:
        base = os.path.dirname(str(i))
        deps += env.Install('../../'+base, i)
    return deps

#==============================================================================
# Create static library and copy built files to installPath
#==============================================================================
def createStaticLibrary(env, ppf, parentFolder, moduleName):
    """
    Create a static library for a list of pre-processed Fortran and cpp files,
    and return the name of the static library
    """
    if isinstance(ppf[0], list):
        nchunks = len(ppf)
        elsaprod = os.getenv("ELSAPROD")
        staticLib = "lib{}.a".format(moduleName.lower())
        staticLibPath = os.path.join("build", elsaprod, staticLib)
        chunkedStaticLib = "lib{}{:d}.a"
        mergeL = "create {}\n".format(staticLibPath)
        for c in range(nchunks):
            mergeL += "addlib build/{}/lib{}{:d}.a\n".format(
                elsaprod, moduleName, c+1)
        mergeL += "save\nend"
        filename = os.path.join(parentFolder, 'merge.l')
        with open(filename, 'w') as f: f.write(mergeL)
        env.Command(
            staticLib,
            [chunkedStaticLib.format(moduleName, c+1) for c in range(nchunks)] + ['merge.l'],
            "ar -M < merge.l"
        )
        for c in range(nchunks):
            env.StaticLibrary(chunkedStaticLib.format(moduleName, c+1), ppf[c])
    else:
        staticLib = env.StaticLibrary(ppf)
    return staticLib

def copyBuiltFiles(env, staticLib, moduleName, installPath):
    # Copy built files and python files to the install folder
    modDir = os.path.join(installPath, moduleName)
    dp1 = env.Install(modDir, staticLib)
    dp2 = env.Install(modDir, glob.glob('{}/*.py'.format(moduleName)))
    env.Alias(target="install", source=[dp1, dp2])

#==============================================================================
if __name__ == "__main__":
    checkAll()
