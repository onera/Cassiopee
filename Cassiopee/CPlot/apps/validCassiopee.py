# *Cassiopee* GUI for validation and tests
import os, sys, re, glob, signal, platform
import socket
import json
import numpy as np
import subprocess
import threading
import time
import KCore
import KCore.Dist as Dist

# System
mySystem = Dist.getSystem()[0]

# Support MPI?
try:
    import mpi4py
    isMpi = True
except: isMpi = False

# Regexprs
regDiff = re.compile('DIFF')
regFailed = re.compile('FAILED')
regError = re.compile("|".join(['Error', 'Erreur', 'Aborted', 'Abandon', 'Segmentation', 'ERROR: AddressSanitizer', 'getFromArray', 'incoherency']), re.UNICODE)
regLeakError = re.compile('ERROR: LeakSanitizer')
separator = ':'
separatorl = separator+' '
expTest1 = re.compile("_t[0-9]+") # normal tests
expTest2 = re.compile(".py")
expTest3 = re.compile(r"\~")
expTest4 = re.compile("_m[0-9]+") # distributed tests

# Liste des tous les tests obtenue en listant les repertoires
# Un element de la liste est une string comme affichee dans la listbox
TESTS = []
# Test filter: 0 (no filter), 1: sequential tests only, 2: distributed tests only
TESTS_FILTER = 0
# Test metadata: ref info (date, hashes), tag, duration
TESTMETA = {}
# Have test metadata been updated?
TESTMETA_UPDATE = False

# Name of the data folders
DATA = None
# CFD Base
CFDBASEPATH = os.path.join('Validation', 'Cases')
# Paths for 'module' source tests and test/Data folder
MODULESDIR = {'LOCAL': {}, 'GLOBAL': {}}
# Paths to the local and global ValidData folders
VALIDDIR = {'LOCAL': None, 'GLOBAL': None}
# Base used for comparisons. Default is the global base.
BASE4COMPARE = 'GLOBAL'

# User settings dictionary
PREFS = {
    "editor": "emacs"
}

# Si THREAD est None, les test unitaires ne tournent pas
# Sinon, THREAD vaut le thread lance
THREAD = None
# Le process lance sinon None
PROCESS = None
# True if the GUI is used (interactive) or False (command line execution)
INTERACTIVE = len(sys.argv) == 1
# Use Address Sanitizer (ASan) and Leak Sanitizer (LSan) in DEBUG mode only
USE_ASAN = [False, False]

# Est egal a 1 si on doit s'arreter
STOP = 0

# WIDGETS dict
WIDGETS = {}

# Sort test strings
SORT_CATEGORIES = ['Name', 'CPU time', 'Ref. CPU time', 'Date', 'Coverage', 'Tag', 'Status', 'CPU time relDiff.']
SORT_BY = 'Name'
REV_SORT = False


#==============================================================================
# Classes used to by-pass tkinter in the absence of a display environment or
# for command-line execution
#==============================================================================
class NoDisplayListbox:
    def __init__(self, *args, **kwargs):
        self._data = []
        self._active = [] # Mask associated with the data
    def grid(self, *args, **kwargs): pass
    def config(self, *args, **kwargs): pass
    def update(self, *args, **kwargs): pass
    def yview(self, *args, **kwargs): pass

    def insert(self, pos, entry):
        if isinstance(pos, int):
            self._data.insert(pos, entry)
            self._active.insert(pos, False)
        else:
            self._data.append(entry)
            self._active.append(False)

    def delete(self, pos, posEnd=None):
        if not self._data: return
        pos = int(pos)
        if posEnd is None or posEnd == pos:
            # Deleting a single entry
            del self._data[pos]
            del self._active[pos]
            return
        # Deleting entries ...
        ndata = len(self._data)
        if isinstance(posEnd, str):
            # ... from pos till the end
            posEnd = ndata
        else:
            # ... given a range of indices
            posEnd = min(int(posEnd), ndata)
        delIds = list(range(pos, posEnd))
        self._data = [self._data[i] for i in range(ndata) if i not in delIds]
        self._active = [self._active[i] for i in range(ndata) if i not in delIds]

    def selection_set(self, pos, posEnd=None):
        pos = int(pos)
        ndata = len(self._data)
        if posEnd is None:
            # 1 new active entry
            self._active[pos] = True
            return
        if isinstance(posEnd, str):
            # Active entries from pos till the end
            posEnd = ndata
        else:
            # A range of new active entries
            posEnd = min(int(posEnd), ndata)
        for i in range(pos, posEnd): self._active[i] = True

    def curselection(self):
        return [i for i, state in enumerate(self._active) if state]
    def get(self, pos): return self._data[pos]

class NoDisplayIntVar:
    def __init__(self, value, *args, **kwargs): self._value = int(value)
    def set(self, value): self._value = int(value)
    def get(self): return self._value

class NoDisplayStringVar:
    def __init__(self, *args, **kwargs): self._filters = ""
    def set(self, filters): self._filters = str(filters)
    def get(self): return self._filters

class NoDisplayLabel:
    def __init__(self, *args, **kwargs): pass
    def grid(self, *args, **kwargs): pass
    def config(self, *args, **kwargs): pass
    def update(self, *args, **kwargs): pass

class NoDisplayButton:
    def __init__(self, *args, **kwargs): pass
    def configure(self, *args, **kwargs): pass

class NoDisplayEntry:
    def __init__(self, *args, **kwargs): pass
    def grid(self, *args, **kwargs): pass
    def bind(self, *args, **kwargs): pass
    def update(self, *args, **kwargs): pass

class NoDisplayScrollbar:
    def __init__(self, *args, **kwargs): pass
    def grid(self, *args, **kwargs): pass
    def config(self, *args, **kwargs): pass
    def set(self): pass


#==============================================================================
# Get installation paths of Cassiopee, Fast and all PModules
#==============================================================================
def isDBAdmin():
    import getpass
    return getpass.getuser() == "cassiope"

def getInstallPaths():
    try:
        # Check installPath
        import KCore.installPath
        cassiopeeIncDir = KCore.installPath.includePath
        cassiopeeIncDir = os.path.dirname(cassiopeeIncDir)
    except ImportError:
        raise SystemError("Error: KCore module is required to use this script.")
    try:
        import FastC.installPath
        fastIncDir = FastC.installPath.includePath
        fastIncDir = os.path.dirname(fastIncDir)
    except:
        fastIncDir = None
    parentDirname = os.path.dirname(os.path.dirname(cassiopeeIncDir))
    pmodulesDir = os.path.join(parentDirname, 'PModules')
    if not os.path.isdir(pmodulesDir): pmodulesDir = None
    fspluginsDirs = []
    if '_coda' in os.getenv("ELSAPROD"):
        fspluginsDir = os.path.join(parentDirname, 'FSPlugins')
        if not os.path.isdir(fspluginsDir):
            fspluginsDir = os.getenv('FSPLUGINS', None)
        if fspluginsDir is not None:
            for d in os.listdir(fspluginsDir):
                absd = os.path.join(fspluginsDir, d)
                if d.startswith('FS') and os.path.isdir(absd):
                    pyPath = os.path.join(absd, 'py')
                    if os.path.isdir(pyPath):
                        fspluginsDirs.append(os.path.abspath(pyPath))
    return cassiopeeIncDir, fastIncDir, pmodulesDir, *fspluginsDirs

def checkEnvironment():
    # Check environment
    cassiopee = os.getenv('CASSIOPEE')
    if cassiopee is None or cassiopee == '':
        print('Error: CASSIOPEE must be present in your environment.')
        sys.exit()

    # Cannot work because of symbolic links to prods on juno
    #if os.path.join(cassiopee, "Cassiopee") != getInstallPaths()[0]:
    #    print("Error: Path mismatch between $CASSIOPEE and KCore/installPath")
    #    sys.exit()

#==============================================================================
# Simulate check_output since it doesn't exist for early version of python
# Retourne le resultat de cmd comme une string
#==============================================================================
def check_output(cmd, shell, stderr):
    global PROCESS
    version = sys.version_info
    version0 = version[0]
    version1 = version[1]
    mode = 4

    #if (version0 == 2 and version1 >= 7) or (version0 == 3 and version1 >= 2) or version0 > 3:

    if mode == 0: # avec check_output
        out = subprocess.check_output(cmd, shell=shell, stderr=stderr)
        return out
    elif mode == 1: # avec run
        PROCESS = subprocess.run(cmd, check=True, shell=shell, stderr=stderr,
                                 stdout=subprocess.PIPE)
        return PROCESS.stdout
    elif mode == 2: # avec Popen + python 2.7
        import shlex
        cmd = shlex.split(cmd)

        wdir = '.'
        # modifie cd en working dir
        if cmd[0] == 'cd': wdir = cmd[1]; cmd = cmd[3:]
        PROCESS = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, cwd=wdir)
        out = ''
        while True:
            line = PROCESS.stdout.readline()
            if line != '': out += line
            else: break
        ti = ''
        while True:
            line = PROCESS.stderr.readline()
            if line != '': ti += line
            else: break
        # change le retour de time pour etre identique a celui du shell
        i1 = ti.find('elapsed')
        i2 = ti.find('system')
        if i1 != -1 and i2 != -1:
            ti = 'real '+ti[i2+7:i1]
            ti = ti.replace(':', 'm')
            ti += 's'
            out += ti
        return out

    elif mode == 3: # avec Popen + python 3
        cmd = cmd.split(' ')
        wdir = '.'
        # modifie cd en working dir
        if cmd[0] == 'cd': wdir = cmd[1]
        if mySystem == 'windows' or mySystem == 'mingw': cmd = cmd[3:]
        else: cmd = cmd[2:]
        if wdir[-1] == ';': wdir = wdir[:-1]
        PROCESS = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, cwd=wdir, shell=shell)
        out = b''
        while True:
            line = PROCESS.stdout.readline()
            if line != b'': out += line
            else: break
        ti = b''
        while True:
            line = PROCESS.stderr.readline()
            if line != b'': ti += line
            else: break
        # change le retour de time pour etre identique a celui du shell
        i1 = ti.find(b'elapsed')
        i2 = ti.find(b'system')
        if i1 != -1 and i2 != -1:
            ti = b'real '+ti[i2+7:i1]
            ti = ti.replace(b':', b'm')
            ti += b's'
            out += ti
        return out

    elif mode == 4: # inspire de python
        wdir = '.'; ossid = None
        if mySystem == 'windows' or mySystem == 'mingw': ossid = None
        else: ossid = os.setsid
        PROCESS = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, cwd=wdir,
                                   shell=shell, preexec_fn=ossid)

        # max accepted time is between 2 to 6 minutes
        nthreads = float(Threads.get())
        timeout = (100. + 120.*Dist.DEBUG)*(1. + 4.8/nthreads)
        stdout, stderr = PROCESS.communicate(None, timeout=timeout)

        if PROCESS.wait() != 0: stderr += b'\nError: process FAILED (Segmentation Fault or floating point exception).'
        PROCESS = None # fini!
        return stdout+stderr

# retourne une chaine justifiee en fonction de la font et
# d'une taille voulue
def ljust(text, size):
    if generalFontFixed == 1:
        # for mono fonts (faster)
        form = '{:%d}'%size
        return form.format(text)
    else:
        l = generalFont.measure(text)*1.
        l = int(round((size*generalFontA-l)/generalFontS))
        if l > 0:
            form = '{}{:%d}'%l
            return form.format(text, ' ')
        else: return text

#==============================================================================
# build a test string:
# 0       1         2         3             4     5         6    7
# module, testname, CPU time, ref CPU time, date, coverage, tag, status
# IN: module, test: test concerne
# IN: CPUtime: nouveau CPU time
# IN: coverage: nouveau coverage
# IN: tag: nouveau tag
# IN: status: nouveau status
# Recupere les anciennes donnees dans les fichiers time & star
#==============================================================================
def buildString(module, test, CPUtime='...', coverage='...%', status='...',
                tag=' '):
    global TESTMETA_UPDATE
    testName = module+'/'+test
    refDate = '...'
    if testName not in TESTMETA:
        TESTMETA_UPDATE = True
        refCoverage = '...%'
        refCPUtime = '...'
        refTag = ' '
    else:
        if TESTMETA[testName].get('ref', []):
            refDate = TESTMETA[testName]['ref'][0].get('date', '...')
        refCoverage = TESTMETA[testName].get('coverage', '...%')
        refCPUtime = TESTMETA[testName].get('time', '...')
        refTag = TESTMETA[testName].get('tag', ' ')

    execTime = '../../.. ..h..'
    if status != '...': # Not First call
        execTime = time.strftime('%d/%m/%y %Hh%M', time.localtime())

    if coverage == '...%': coverage = refCoverage
    if tag == ' ': tag = refTag

    s = ljust(module, 16)+separatorl+ljust(test, 40)+separatorl+\
        ljust(CPUtime, 10)+separatorl+ljust(refCPUtime, 10)+separatorl+\
        ljust(refDate, 16)+separatorl+ljust(coverage, 5)+separatorl+\
        tag.ljust(2)+separatorl+' '+ljust(status, 10)
    return s

#==============================================================================
# Set paths in MODULESDIR for each module of Cassiopee, Fast, PModules and
# CFDBase both locally and globally. Set paths for ValidData directories
#==============================================================================
def setPaths():
    def _setModuleDirs(*args, **kwargs):
        global MODULESDIR
        loc = kwargs.get('loc', 'LOCAL')
        if loc not in ['LOCAL', 'GLOBAL']: loc = 'LOCAL'
        notTested = ['Upmost', 'FastP']
        cassiopeeIncDir = args[0]
        paths = list(args)[1:]
        for path in paths:
            if path is None: continue
            print('Info: getting {} module tests in: {}.'.format(loc.lower(), path))
            try: mods = os.listdir(path)
            except: mods = []
            for mod in mods:
                if mod not in notTested and mod not in MODULESDIR[loc]:
                    if loc == 'GLOBAL' and mod not in MODULESDIR['LOCAL']:
                        # Skip modules which aren't found locally - would hang
                        continue
                    pathMod = os.path.join(path, mod)
                    a = os.access(os.path.join(pathMod, 'test'), os.F_OK) # PModules svn
                    if a: MODULESDIR[loc][mod] = path
                    else:
                        a = os.access(os.path.join(pathMod, mod, 'test'), os.F_OK) # PModules git
                        if a: MODULESDIR[loc][mod] = pathMod

        print('Info: getting {} module names in: {}.'.format(loc.lower(), cassiopeeIncDir))
        try: mods = os.listdir(cassiopeeIncDir)
        except: mods = []
        for mod in mods:
            if mod not in MODULESDIR[loc]:
                a = os.access(os.path.join(cassiopeeIncDir, mod, 'test'), os.F_OK)
                if a: MODULESDIR[loc][mod] = cassiopeeIncDir

        # Validation CFD
        MODULESDIR[loc]['CFDBase'] = os.path.dirname(os.path.dirname(cassiopeeIncDir))

    global VALIDDIR
    # Module paths when the local base is used
    allPackageDirs = getInstallPaths()
    _setModuleDirs(*allPackageDirs, loc='LOCAL')

    # Local valid paths
    VALIDDIR['LOCAL'] = os.path.join(os.getenv('CASSIOPEE'), 'Cassiopee',
                                     'Valid{}'.format(DATA))
    if not os.access(VALIDDIR['LOCAL'], os.W_OK):
        VALIDDIR['LOCAL'] = os.path.join(os.getcwd(), "Valid{}".format(DATA))
        if not os.path.isdir(VALIDDIR['LOCAL']): os.mkdir(VALIDDIR['LOCAL'])

    # Module paths when the global base is used
    parentDirname = os.path.join('/stck', 'cassiope', 'git')
    if getDBInfo():  # check if a global base exists
        cassiopeeIncDir = os.path.join(parentDirname, 'Cassiopee', 'Cassiopee')
        fastIncDir = os.path.join(parentDirname, 'Fast', 'Fast')
        if not os.path.isdir(fastIncDir): fastIncDir = None
        pmodulesIncDir = os.path.join(parentDirname, 'PModules')
        if not os.path.isdir(pmodulesIncDir): pmodulesIncDir = None
        _setModuleDirs(cassiopeeIncDir, fastIncDir, pmodulesIncDir, loc='GLOBAL')

        # Global valid paths
        VALIDDIR['GLOBAL'] = os.path.join(parentDirname, 'Cassiopee',
                                          'Cassiopee', 'Valid{}'.format(DATA))

#==============================================================================
# Retourne la liste des modules situes dans Cassiopee, Fast et PModules
# Eventuellement peut ajouter "CFDBase", nom referencant les tests
# de validation des solveurs (CFDBase)
#==============================================================================
def getModules():
    return sorted(MODULESDIR[BASE4COMPARE].keys())

#==============================================================================
# Retourne la liste des tests unitaires d'un module
# si module == 'CFDBase', retourne la liste des cas de validation CFD
#==============================================================================
def getTests(module):
    if module == 'CFDBase': return getCFDBaseTests()
    return getUnitaryTests(module)

#==============================================================================
# Retourne la liste des tests unitaires pour un module donne
# Les tests unitaires doivent etre dans module/test
# La variable globale TESTS_FILTER permet de filtrer tests sequentiels et tests
# distribues
#==============================================================================
def getUnitaryTests(module):
    modulesDir = MODULESDIR[BASE4COMPARE][module]
    path = os.path.join(modulesDir, module, 'test')
    files = os.listdir(path)
    tests = []
    for f in files:
        m2 = expTest2.search(f)
        if m2 is None: continue
        m3 = expTest3.search(f)
        if m3 is not None: continue
        if f[0] == '#': continue
        m1 = expTest1.search(f)
        m4 = expTest4.search(f)
        if m1 is not None and TESTS_FILTER != 2: tests.append(f) # test seq
        elif isMpi and m4 is not None and TESTS_FILTER != 1: tests.append(f) # test mpi
    return sorted(tests)

#==============================================================================
# Retourne la liste des cas de validation CFD (CFDBase)
# Il doivent etre dans Validation/Cases
#==============================================================================
def getCFDBaseTests():
    path = os.path.join(MODULESDIR[BASE4COMPARE]['CFDBase'], CFDBASEPATH)
    try: reps = os.listdir(path)
    except: reps = []
    tests = []
    for r in reps: # a terme a supprimer
        if r == 'NACA': tests.append(r) # MB 2D Euler
        elif r == 'NACA_IBC': tests.append(r) # IBC 2D Euler
        elif r == 'DAUPHIN': tests.append(r) # MB 3D Euler
        elif r == 'FLATPLATE': tests.append(r) # MB 3D SA
        elif r == 'RAE2822': tests.append(r) # MB 2D SA
        elif r == 'RAE2822_IBC': tests.append(r) # IBC 2D SA
        elif r == 'CUBE_IBC': tests.append(r) # IBC 3D SA
    return sorted(tests)

#==============================================================================
# Ecrit un fichier contenant date, CPUtime, coverage
#==============================================================================
def writeTime(fileTime, CPUtime, coverage):
    try:
        execTime = time.strftime('%d/%m/%y %Hh%M', time.localtime())
        with open(fileTime, 'w') as f:
            f.write(execTime + '\n')
            f.write(CPUtime + '\n')
            f.write(coverage + '\n')
    except: pass

#==============================================================================
# Ecrit un fichier contenant date, machine, nbre de threads, git info
# et logTxt
#==============================================================================
def writeFinal(filename, gitInfo="", logTxt=None, append=False):
    execTime = time.strftime('%d/%m/%y %Hh%M', time.localtime())
    machine = platform.uname()
    if len(machine) > 1: machine = machine[1]
    else: machine = 'Unkwown'
    nthreads = Threads.get()
    mode = 'a' if append else 'w'
    with open(filename, mode) as f:
        f.write(execTime+'\n')
        f.write(machine+'\n')
        f.write(nthreads+'\n')
        if gitInfo: f.write(gitInfo+'\n')
        if logTxt is not None: f.write(logTxt+'\n')

#==============================================================================
# Read un fichier star
#==============================================================================
def readStar(fileStar):
    star = ' '
    try:
        with open(fileStar, 'r') as f:
            star = f.readline().rstrip('\n')
    except: pass
    return star

#==============================================================================
# Lance un seul test unitaire ou un cas de la base de validation
#==============================================================================
def runSingleTest(no, module, test, update=False):
    if module == 'CFDBase': return runSingleCFDTest(no, module, test, update)
    return runSingleUnitaryTest(no, module, test, update)

#==============================================================================
# extrait le temps CPU d'un chaine output (utile sur windows)
# retourne le temps CPU comme une chaine
# moyenne si plusieurs repetitions d'un meme cas unitaire
#==============================================================================
def extractCPUTimeWindows(output1, output2):
    CPUtime = 'Unknown'
    try:
        split1 = output1.split(':')
        h1 = int(split1[0])
        m1 = int(split1[1])
        s1 = split1[2]; s1 = s1.split(',')
        ms1 = int(s1[1])
        s1 = int(s1[0])
        t1 = h1*3600. + m1*60. + s1 + 0.01*ms1
        split2 = output2.split(':')
        h2 = int(split2[0])
        m2 = int(split2[1])
        s2 = split2[2]; s2 = s2.split(',')
        ms2 = int(s2[1])
        s2 = int(s2[0])
        t2 = h2*3600. + m2*60. + s2 + 0.01*ms2
        tf = (t2-t1)
        hf = int(tf/3600.)
        tf = tf - 3600*hf
        mf = int(tf/60.)
        tf = tf - 60*mf
        sf = int(tf*100)/100.
        if hf > 0: CPUtime = '%dh%dm%gs'%(hf,mf,sf)
        else: CPUtime = '%dm%gs'%(mf,sf)
    except: pass
    return CPUtime

#=============================================================================
# Extrait le temps CPU d'une sortie time -p (unix)
# Moyenne si plusieurs repetitions d'un meme cas unitaire
#=============================================================================
def extractCPUTimeUnix(output):
    CPUtime = 'Unknown'
    try:
        i1 = output.find('real')
        if i1 == -1: # external time command
            i1 = output.find('elapsed')
            output = output[:i1+7].strip().replace(',', '.')
            output = re.split('\n| ', output)
            output = [s for s in output if s.endswith('elapsed')][0]
            output = re.sub('[a-z]', "", output)
        else: # shell built-in time command (shell keyword)
            output = output[i1+4:].strip()
            output = output.replace(',', '.').replace('s', '')
            output = re.sub('[h,m]', ':', output)
            output = re.split('\n|\t| ', output)[0]
        output = output.split(':')
        output = output[::-1]
        dt = [1., 60., 3600.]
        tf = 0.
        for i in range(len(output)):
            tf += dt[i]*float(output[i])
        hf = tf//3600
        tf = tf%3600
        output = [int(hf), int(tf//60), float(tf%60)]
        if hf > 0.: CPUtime = '{:d}h{:d}m{:.2f}'.format(*output)
        elif output[1] > 0.: CPUtime = '{:d}m{:.2f}'.format(*output[1:])
        else: CPUtime = '0m{:.2f}'.format(output[-1])
        CPUtime = CPUtime.rstrip('0').rstrip('.') + 's'
    except: pass
    return CPUtime

#==============================================================================
# Lance un seul test unitaire de module
#==============================================================================
def runSingleUnitaryTest(no, module, test, update=False):
    global TESTS
    testr = os.path.splitext(test)
    modulesDir = MODULESDIR[BASE4COMPARE][module]
    path = os.path.join(modulesDir, module, 'test')
    testName = module+'/'+test

    m1 = expTest1.search(test) # seq (True) ou distribue (False)

    nthreads = KCore.kcore.getOmpMaxThreads()
    nthreads = int(Threads.get())
    bktest = "bk_{0}".format(test) # backup

    if mySystem == 'mingw' or mySystem == 'windows':
        # Commande Dos (sans time)
        path = path.replace('/', '\\')
        pythonExec = os.getenv('PYTHONEXE', 'python')
        if m1 is not None: cmd = 'cd %s && %s %s'%(path, pythonExec, test)
        else: cmd = 'cd %s && set OMP_NUM_THREADS=%d && mpiexec -np 2 %s %s'%(path, nthreads//2, pythonExec, test)
        cmd2 = 'echo %time%'
    else:
        # Unix - le shell doit avoir l'environnement cassiopee
        #sformat = r'"real\t%E\nuser\t%U\nsys\t%S"'
        if m1 is None:
            sanitizerFlag = '' # LSAN always return a seg fault in parallel
            cmd = 'cd %s; time kpython %s -n 2 -t %d %s'%(
                path, sanitizerFlag, nthreads//2, test)
        else:
            sanitizerFlag = '-s' if any(USE_ASAN) else ''
            if module.startswith('FS'):
                cmd = 'cd %s; time kpython %s -n 1 -t %d %s'%(
                    path, sanitizerFlag, nthreads, test)
            else:
                cmd = 'cd %s; time kpython %s -t %d %s'%(
                    path, sanitizerFlag, nthreads, test)

    try:
        if mySystem == 'mingw' or mySystem == 'windows':
            output1 = check_output(cmd2, shell=True, stderr=subprocess.STDOUT)
            if sys.version_info[0] == 3: output1 = output1.decode()
        output = check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        if sys.version_info[0] == 3: output = output.decode()

        if mySystem == 'mingw' or mySystem == 'windows':
            output2 = check_output(cmd2, shell=True, stderr=subprocess.STDOUT)
            if sys.version_info[0] == 3: output2 = output2.decode()

        print(output)

        # Recupere success/failed
        success = 0
        if regLeakError.search(output) is not None: success = 2 # always check first
        if regDiff.search(output) is not None: success = 1
        if regFailed.search(output) is not None: success = 1
        if regError.search(output) is not None: success = 1

        # Recupere le CPU time
        if mySystem == 'mingw' or mySystem == 'windows':
            CPUtime = extractCPUTimeWindows(output1, output2)
        else:
            CPUtime = extractCPUTimeUnix(output)

        # Recupere le coverage
        i1 = output.find('coverage=')
        if i1 == -1: coverage = '0%'
        else:
            sub = output[i1+9:i1+13]
            i1 = sub.find('%')
            coverage = sub[:i1+1]
            coverage = coverage.strip()

    except subprocess.TimeoutExpired:
        # killed here because timeout of communicate doesnt kill child processes
        if mySystem == 'mingw' or mySystem == 'windows':
            subprocess.call(['taskkill', '/F', '/T', '/PID', str(PROCESS.pid)])
        else: # unix
            # try soft, then hard
            os.killpg(os.getpgid(PROCESS.pid), signal.SIGTERM)
            os.kill(PROCESS.pid, signal.SIGTERM)
            os.killpg(os.getpgid(PROCESS.pid), signal.SIGKILL)
            os.kill(PROCESS.pid, signal.SIGKILL)
        print('\nError: process TIMED OUT (killed).')
        success = 3; CPUtime = 'Unknown'; coverage='0%'

    except Exception as e:
        print(e)
        success = 1; CPUtime = 'Unknown'; coverage='0%' # Core dump/error

    # update le fichier .time (si non present)
    fileTime = '%s/%s/%s/%s/%s.time'%(MODULESDIR['LOCAL'][module], module, 'test', DATA, testr[0])
    if not os.access(fileTime, os.F_OK):
        writeTime(fileTime, CPUtime, coverage)

    if update or (testName not in TESTMETA) or (not TESTMETA[testName]['ref']):  # Update test metadata
        updateTestMetadata(module, test, CPUtime)

    # Recupere le tag local
    tag = TESTMETA[testName]['tag']

    # update status
    if success == 0: status = 'OK'
    elif success == 2: status = 'MEMLEAK'
    elif success == 3: status = 'TIMEOUT'
    else: status = 'FAILED'
    s = buildString(module, test, CPUtime, coverage, status, tag)
    regTest = re.compile(' '+test+' ')
    regModule = re.compile(module+' ')
    for c, tt in enumerate(TESTS):
        if regModule.search(tt) is not None:
            if regTest.search(tt) is not None: TESTS[c] = s; break
    Listbox.delete(no, no)
    Listbox.insert(no, s)
    Listbox.update()
    CPUtime = string2Time(CPUtime)
    return CPUtime

#==============================================================================
# Lance un seul test de la base CFD (CFDBase)
# module = 'CFDBase'
# test = nom du repertoire du cas CFD
#==============================================================================
def runSingleCFDTest(no, module, test, update=False):
    global TESTS
    print('Info: Running CFD test %s.'%test)
    path = os.path.join(MODULESDIR[BASE4COMPARE]['CFDBase'], CFDBASEPATH, test)
    testName = module+'/'+test

    m1 = None # si False=seq
    # force mpi test pour certains cas
    if test == 'RAE2822_IBC': m1 = True

    if m1 is not None:
        try: import mpi4py
        except: m1 = None

    nthreads = KCore.kcore.getOmpMaxThreads()
    nthreads = int(Threads.get())

    if mySystem == 'mingw' or mySystem == 'windows':
        # Commande Dos (sans time)
        path = path.replace('/', '\\')
        if m1 is None: cmd = 'cd %s && ./valid check'%(path)
        else: cmd = 'cd %s && ./valid check 0 0 0 2 %d'%(path, nthreads//2)
        cmd2 = 'echo %time%'
    else:
        # Unix - le shell doit avoir l'environnement cassiopee
        if m1 is None: cmd = 'cd %s; ./valid check'%(path)
        else: cmd = 'cd %s; ./valid check 0 0 0 2 %d'%(path, nthreads//2)
    try:
        if mySystem == 'mingw' or mySystem == 'windows':
            output1 = check_output(cmd2, shell=True, stderr=subprocess.STDOUT)
            if sys.version_info[0] == 3: output1 = output1.decode()
        output = check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        if sys.version_info[0] == 3: output = output.decode()
        if mySystem == 'mingw' or mySystem == 'windows':
            output2 = check_output(cmd2, shell=True, stderr=subprocess.STDOUT)
            if sys.version_info[0] == 3: output2 = output2.decode()

        print(output)

        # Recupere success/failed
        success = 0
        if regLeakError.search(output) is not None: success = 2 # always check first
        if regDiff.search(output) is not None: success = 1
        if regFailed.search(output) is not None: success = 1
        if regError.search(output) is not None: success = 1

        # Recupere le CPU time
        if mySystem == 'mingw' or mySystem == 'windows':
            CPUtime = extractCPUTimeWindows(output1, output2)
        else:
            CPUtime = extractCPUTimeUnix(output)
        # Recupere le coverage
        coverage = '100%'

    except subprocess.TimeoutExpired:
        # killed here because timeout of communicate doesnt kill child processes
        if mySystem == 'mingw' or mySystem == 'windows':
            subprocess.call(['taskkill', '/F', '/T', '/PID', str(PROCESS.pid)])
        else: # unix
            # try soft, then hard
            os.killpg(os.getpgid(PROCESS.pid), signal.SIGTERM)
            os.kill(PROCESS.pid, signal.SIGTERM)
            os.killpg(os.getpgid(PROCESS.pid), signal.SIGKILL)
            os.kill(PROCESS.pid, signal.SIGKILL)
        print('\nError: process TIMED OUT (killed).')
        success = 3; CPUtime = 'Unknown'; coverage='0%'

    except Exception as e:
        print(e)
        success = 1; CPUtime = 'Unknown'; coverage='0%' # Core dump/error

    # update le fichier .time (si non present)
    fileTime = '%s/%s/%s.time'%(path, DATA, test)
    if not os.access(fileTime, os.F_OK):
        writeTime(fileTime, CPUtime, coverage)

    if update or (testName not in TESTMETA) or (not TESTMETA[testName]['ref']):  # Update test metadata
        updateTestMetadata(module, test, CPUtime)

    # Recupere le tag local
    tag = TESTMETA[testName]['tag']

    # update status
    if success == 0: status = 'OK'
    elif success == 2: status = 'MEMLEAK'
    elif success == 3: status = 'TIMEOUT'
    else: status = 'FAILED'
    s = buildString(module, test, CPUtime, coverage, status, tag)
    regTest = re.compile(' '+test+' ')
    regModule = re.compile(module+' ')
    for c, tt in enumerate(TESTS):
        if regModule.search(tt) is not None:
            if regTest.search(tt) is not None: TESTS[c] = s; break
    Listbox.delete(no, no)
    Listbox.insert(no, s)
    Listbox.update()
    CPUtime = string2Time(CPUtime)
    return CPUtime

#==============================================================================
# Recupere le nbre de tests selectionnes et le temps total correspondant
#==============================================================================
def getTestsTime():
    selection = Listbox.curselection()
    total = len(selection)
    remaining = 0.
    for s in selection:
        t = Listbox.get(s)
        splits = t.split(separator)
        remaining += string2Time(splits[3])
    return (total, remaining)

#==============================================================================
# Run selected tests
# Update TESTS, update listbox, update progression
#==============================================================================
def runTests(update=False):
    global STOP, THREAD
    selection = Listbox.curselection()
    displayStatus(1)
    current = 0
    (total, remaining) = getTestsTime()
    elapsed = 0.

    for s in selection:
        no = int(s)
        t = Listbox.get(s)
        splits = t.split(separator)
        module = splits[0]
        test = splits[1]
        module = module.strip()
        test = test.strip()
        if update:  # Delete reference
            modulesDir = MODULESDIR[BASE4COMPARE][module]
            if module == 'CFDBase':
                pathl = os.path.join(modulesDir, CFDBASEPATH, test)
                testref = 'post.ref*'
            else:
                pathl = os.path.join(modulesDir, module, 'test')
                testref = os.path.splitext(test)[0] + '.ref*'
            rmFile(pathl, testref)
        current += 1; displayProgress(current, total, remaining, elapsed)
        remaining -= string2Time(splits[3])
        CPUtime = runSingleTest(no, module, test, update)
        elapsed += CPUtime # real elapsed time
        if STOP == 1: STOP = 0; displayStatus(0); return
    displayStatus(0)
    THREAD=None
    writeSessionLog()

def runTestsInThread(update=False):
    global THREAD, STOP
    if THREAD is not None: return
    STOP = 0
    THREAD = threading.Thread(target=runTests, args=(update,))
    THREAD.start()

#==============================================================================
# Load, fill, update, and write test metadata
#==============================================================================
def loadTestMetadata(loc=None):
    global TESTMETA
    if loc is None: loc = BASE4COMPARE
    metaFilename = os.path.join(VALIDDIR[loc], "testMetadata.json")
    if os.access(metaFilename, os.R_OK):
        with open(metaFilename, 'r') as f: TESTMETA = json.load(f)
        if loc == 'GLOBAL':  # overwrite global tags with local tags
            metaFilename = os.path.join(VALIDDIR['LOCAL'], "testMetadata.json")
            if os.access(metaFilename, os.R_OK):
                with open(metaFilename, 'r') as f: tmpMeta = json.load(f)
                for testName in TESTMETA:
                    if testName in tmpMeta:
                        TESTMETA[testName]['tag'] = tmpMeta[testName].get('tag', ' ')
    else: fillTestMetadata()

def newMetadata(coverage='...%', tag=' ', time='...'):
    return {'coverage': coverage, 'tag': tag, 'ref': [], 'time': time}

def newRefMetadata(date='...', hashCassiopee='Unknown', hashFast='Unknown'):
    return {'date': date, 'hashCassiopee': hashCassiopee, 'hashFast': hashFast}

def fillTestMetadata():
    global TESTMETA, TESTMETA_UPDATE
    TESTMETA = {}
    TESTMETA_UPDATE = True
    # Loop over all modules
    for module in getModules():
        tests = getTests(module)
        # Loop over all tests of that module
        for test in tests:
            # Get local tag and reference time
            modulesDirCmp = MODULESDIR[BASE4COMPARE][module]
            modulesDirLoc = MODULESDIR['LOCAL'][module]
            if module == 'CFDBase':
                fileTime = os.path.join(modulesDirCmp, CFDBASEPATH, test, DATA, test+'.time')
                fileStar = os.path.join(modulesDirLoc, CFDBASEPATH, test, DATA, test+'.star')
            else:
                testr = os.path.splitext(test)
                fileTime = os.path.join(modulesDirCmp, module, 'test', DATA, testr[0]+'.time')
                fileStar = os.path.join(modulesDirLoc, module, 'test', DATA, testr[0]+'.star')

            if os.access(fileTime, os.R_OK):
                with open(fileTime, 'r') as f: fileData = f.read().split('\n')
                if len(fileData) > 0: refDate = fileData[0]
                else: refDate = '...'
                if len(fileData) > 1: refCPUtime = fileData[1]
                else: refCPUtime = '...'
                if len(fileData) > 3 and fileData[2] != '': refCoverage = fileData[2]
                else: refCoverage = '...%'
            else:
                refDate = '...'
                refCPUtime = '...'
                refCoverage = '...%'
            if os.access(fileStar, os.R_OK): refTag = readStar(fileStar)
            else: refTag = ' '

            testName = module+'/'+test
            TESTMETA[testName] = newMetadata(coverage=refCoverage, tag=refTag,
                                             time=refCPUtime)
            if refDate != '...':
                TESTMETA[testName]['ref'].append(newRefMetadata(date=refDate))

def updateTestMetadata(module, test, CPUtime='...'):
    global TESTMETA, TESTMETA_UPDATE
    TESTMETA_UPDATE = True
    testName = module+'/'+test
    cassiopeeIncDir, fastIncDir = getInstallPaths()[:2]
    if fastIncDir is not None: hashFast = Dist.getGitHash(fastIncDir)[:7]
    else: hashFast = None
    refEntry = newRefMetadata(
        date=time.strftime('%d/%m/%y %Hh%M', time.localtime()),
        hashCassiopee=Dist.getGitHash(cassiopeeIncDir)[:7],
        hashFast=hashFast
    )
    if not isinstance(CPUtime, str):
        mins = int(CPUtime//60)
        secs = CPUtime%60
        fullsec = int(secs)
        csecs = int((secs - fullsec)*100)
        CPUtime = "{}m{}.{:02d}s".format(mins, fullsec, csecs)

    if testName in TESTMETA:
        TESTMETA[testName]['time'] = CPUtime
        if 'ref' in TESTMETA[testName]:
            TESTMETA[testName]['ref'].insert(0, refEntry)
            TESTMETA[testName]['ref'] = TESTMETA[testName]['ref'][:5]
        else: TESTMETA[testName]['ref'] = [refEntry]
    else:
        TESTMETA[testName] = newMetadata(time=CPUtime)
        TESTMETA[testName]['ref'] = [refEntry]

def writeTestMetadata():
    if TESTMETA_UPDATE:
        if BASE4COMPARE == 'GLOBAL':  # update tags in local test metadata only
            from copy import deepcopy
            tmpMeta = deepcopy(TESTMETA)
            loadTestMetadata(loc='LOCAL')  # reload original local metadata
            for testName in TESTMETA:
                if testName in tmpMeta:
                    TESTMETA[testName]['tag'] = tmpMeta[testName].get('tag', ' ')

        metaFilename = os.path.join(VALIDDIR['LOCAL'], "testMetadata.json")
        with open(metaFilename, 'w') as f:
            json.dump(TESTMETA, f, indent=4, sort_keys=True)

#==============================================================================
# Supprime un fichier
# IN: path: file path
# IN: file: file name
#==============================================================================
def rmFile(path, fileName):
    if mySystem == 'mingw' or mySystem == 'windows':
        path = path.replace('/', '\\')
        cmd = 'cd '+path+' && del '+DATA+'\\'+fileName
    else:
        cmd = 'cd '+path+'; rm -f '+DATA+'/'+fileName
    try:
        subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)
    except: pass

#==============================================================================
# Construit la liste des tests
# Update TESTS et la listBox
#==============================================================================
def buildTestList(sessionName=None, modules=[]):
    global TESTS
    TESTS = []
    Listbox.delete(0, 'end')
    if not modules:
        modules = getModules()

    if sessionName is not None:
        # Read last sessionLog conditionally
        ncolumns = 8
        logname = sorted(glob.glob(os.path.join(VALIDDIR['LOCAL'], "{}-*.log".format(sessionName))))
        if len(logname): logname = logname[-1]
        else: sessionName = None

    if (sessionName is not None and os.access(logname, os.R_OK) and
            os.path.getsize(logname) > 0):
        print("Loading last session: {}".format(logname))
        with open(logname, "r") as g:
            sessionLog = [line.rstrip().split(':') for line in g.readlines()]
        # Remove header from logfile
        sessionLog = [testLog for testLog in sessionLog
                      if (isinstance(testLog, list) and len(testLog) == ncolumns)]
        if not sessionLog:
            ncolumns = 7
            sessionLog = [testLog for testLog in sessionLog
                          if (isinstance(testLog, list) and len(testLog) == ncolumns)]
        # Create array and remove leading and trailing white spaces
        arr = np.array([entry.strip() for testLog in sessionLog for entry in testLog],
                       dtype=object)
        arr = arr.reshape(-1, ncolumns)

        # Read sessionLog and combine with lastSession. Priority given to
        # data from current session
        ncolumns = 8
        logname = os.path.join(VALIDDIR['LOCAL'], "session.log")
        if os.path.getsize(logname) > 0:
            with open(logname, "r") as g:
                sessionLog = [line.rstrip().split(':') for line in g.readlines()]
            sessionLog = [testLog for testLog in sessionLog
                          if (isinstance(testLog, list) and len(testLog) == ncolumns)]
            if not sessionLog:
                ncolumns = 7
                sessionLog = [testLog for testLog in sessionLog
                              if (isinstance(testLog, list) and len(testLog) == ncolumns)]
            arr2 = np.array([entry.strip() for testLog in sessionLog for entry in testLog],
                            dtype=object)
            arr2 = arr2.reshape(-1, ncolumns)

            testDict = {}
            for t in arr2: testDict[tuple(t[:2])] = t[2:]

            for t in arr:
                key = tuple(t[:2])
                if (key not in testDict) or ('...' in testDict[key]):
                    testDict[key] = t[2:]
            arr = np.array([list(key) + list(data) for key, data in testDict.items()])
    else:
        # Build an empty array
        arr = np.array([], dtype=object)

    for m in modules:
        tests = getTests(m)
        for t in tests:
            if sessionName is not None and arr.size:
                testArr = arr[np.logical_and(arr[:,0] == m, arr[:,1] == t)]
                if testArr.size:
                    # Args are CPU time, Coverage, Status, and Tag if present
                    if ncolumns == 8:
                        if testArr[0][6].strip() in ['OK', 'FAILED', 'MEMLEAK', 'TIMEOUT', '...']:
                            args = testArr[0][[2,5,6,7]]
                        else: args = testArr[0][[2,5,7,6]]
                    else: args = testArr[0][[2,5,6]]
                    s = buildString(m, t, *args)
                else:
                    s = buildString(m, t)
            else:
                s = buildString(m, t)
            TESTS.append(s)
            Listbox.insert('end', s)
    if sessionName is not None and arr.size: writeSessionLog()
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)

#=============================================================================
# Parse some entries of a test string
#=============================================================================
def parseCPUTimeStr(s):
    try:
        s = s.strip().replace('m', ':').replace('s', '')
        minutes, seconds = map(float, s.split(':'))
        return minutes * 60 + seconds  # total time in seconds
    except ValueError:  # undefined as in "..."
        return -1.

def parseCPUTimeRelDiffStr(s1, s2):
    t1 = parseCPUTimeStr(s1)
    t2 = parseCPUTimeStr(s2)
    if t1 <= 0. or t2 <= 0.: return 1e9
    else: return (t2 - t1)/t2

def parseDateStr(s):
    try:
        return time.strptime(s.strip(), "%d/%m/%y %Hh%M")
    except ValueError:  # undefined as in "..."
        return time.strptime(f"01/01/1970 00h00")

def parseCoverageStr(s):
    try:
        return float(s.strip('%'))
    except ValueError:  # undefined as in "..."
        return -1.

#==============================================================================
# Sort tests
#==============================================================================
def splitTestString(test):
    return [entry.strip() for entry in test.split(":")]

def sortTests(tests):
    global SORT_BY
    if SORT_BY not in SORT_CATEGORIES: SORT_BY = 'Name'
    if SORT_BY == 'Name': return sorted(tests, reverse=REV_SORT)
    sortIndex = SORT_CATEGORIES.index(SORT_BY) + 1
    if SORT_BY == 'Coverage':
        return sorted(
            tests,
            key=lambda t: parseCoverageStr(splitTestString(t)[sortIndex]),
            reverse=REV_SORT)
    elif SORT_BY in ['CPU time', 'Ref. CPU time']:
        return sorted(
            tests,
            key=lambda t: parseCPUTimeStr(splitTestString(t)[sortIndex]),
            reverse=REV_SORT)
    elif SORT_BY == 'CPU time relDiff.':
        sortIndex1 = SORT_CATEGORIES.index('CPU time') + 1
        sortIndex2 = SORT_CATEGORIES.index('Ref. CPU time') + 1
        return sorted(
            tests,
            key=lambda t: parseCPUTimeRelDiffStr(
                splitTestString(t)[sortIndex1],
                splitTestString(t)[sortIndex2]),
            reverse=REV_SORT)
    elif SORT_BY == 'Date':
        return sorted(
            tests,
            key=lambda t: parseDateStr(splitTestString(t)[sortIndex]),
            reverse=REV_SORT)
    return sorted(
        tests,
        key=lambda t: splitTestString(t)[sortIndex],
        reverse=REV_SORT)

def sortTestList(event=None, entry='Name', force=False):
    global SORT_BY
    if entry in SORT_CATEGORIES and (force or entry != SORT_BY):
        SORT_BY = entry
        tests = Listbox.get(0, 'end')
        Listbox.delete(0, 'end')
        for s in sortTests(tests): Listbox.insert('end', s)
    updateSortLabels()

def updateSortLabels():
    if not INTERACTIVE: return
    for i, c in enumerate(SORT_CATEGORIES):
        if c == SORT_BY: label = f'(*) By {c}'
        else: label = f'    By {c}'
        sortTab.entryconfig(i, label=label)

def updateSortingOrderLabel(event=None):
    global REV_SORT
    REV_SORT = not REV_SORT
    sortTestList(event=event, entry=SORT_BY, force=True)

#==============================================================================
# Filtre la liste des tests avec la chaine de filter
# Update la listbox
#==============================================================================
def filterTestList(event=None):
    def _rmSubstrings(filters):
        """Remove filters that are a substring of another filter and that as a
           first measure to prevent tests from appearing multiple times in the
           listbox"""
        outFilters = set()
        for str1 in filters:
            # Append if string is not part of a longer string in filters
            if not any(str1 != str2 and str1 in str2 for str2 in filters):
                outFilters.add(str1)
        return outFilters

    def _substituteCustomFilters(filters):
        """Substitute custom keyworded filters comprised between angle brackets
        by their regular expression"""
        outFilters = set()
        for filtr in filters:
            if not filtr: continue
            pos1 = filtr.find('<')
            pos2 = filtr.find('>')
            if pos1 != -1 and pos2 != -1 and pos1 < pos2:
                tmpFiltr = filtr[pos1+1:pos2]
                if filtr[0] == '!':
                    if tmpFiltr == 'SEQ': outFilters.add('&m.$')
                    elif tmpFiltr == 'DIST': outFilters.add('&t.$')
                    elif tmpFiltr == 'RUN': outFilters.update(['&/!FAILED', '&/!MEMLEAK', '&/!TIMEOUT', '&/!OK'])
                    elif tmpFiltr == 'UNRUN': outFilters.update(['/FAILED', '/MEMLEAK', '/TIMEOUT', '/OK'])
                    elif tmpFiltr == 'TAG': outFilters.add(r'@^(?![\*,r,g,b])')
                    elif tmpFiltr == 'UNTAG': outFilters.add(r'@[\*,r,g,b]')
                else:
                    if tmpFiltr == 'SEQ': outFilters.add('&t.$')
                    elif tmpFiltr == 'DIST': outFilters.add('&m.$')
                    elif tmpFiltr == 'RUN': outFilters.update(['/FAILED', '/MEMLEAK', '/TIMEOUT', '/OK'])
                    elif tmpFiltr == 'UNRUN': outFilters.update(['&/!FAILED', '&/!MEMLEAK', '&/!TIMEOUT', '&/!OK'])
                    elif tmpFiltr == 'TAG': outFilters.add(r'@[\*,r,g,b]')
                    elif tmpFiltr == 'UNTAG': outFilters.add(r'@^(?![\*,r,g,b])')
            else: outFilters.add(filtr)
        return outFilters

    filters = Filter.get()
    filters = _rmSubstrings(filters.split(' '))
    filters = _substituteCustomFilters(filters)
    if filters and all(filtr[0] == '&' for filtr in filters):
        filtr0 = filters.pop()
        if len(filtr0) > 1: filters.add(filtr0[1:])

    # Apply filters with an OR gate and append strings to set
    filteredTests = set()
    for filtr in filters:
        if (not filtr) or (filtr in ['#', '/', '!', '@', '%']) or (filtr[0] in ['&', '*']):
            continue
        shift = 1; endidx = 0
        if filtr[0] == '#': pos = 0 # filter modules
        elif filtr[0] == '/': pos = 7 # filter statuses
        elif filtr[0] == '@': pos = 6 # filter tags
        elif filtr[0] == '%': pos = 5 # filter coverage
        else: pos = 1; shift = 0; endidx = -3 # filter test names
        for s in TESTS:
            strg = s.split(separator)[pos].strip()
            if endidx != 0: strg = strg[:endidx]
            try:
                if filtr[shift] == '!':
                    if re.search(filtr[1+shift:], strg) is None:
                        filteredTests.add(s)
                elif re.search(filtr[shift:], strg) is not None:
                    filteredTests.add(s)
            except re.error: pass

    # Apply filters with an AND gate to remove strings from set
    insertedTests = filteredTests.copy()
    for filtr in filters:
        if len(filtr) > 1 and filtr[0] == '&':
            shift = 1; endidx = 0
            if filtr[1] == '#': pos = 0 # filter modules
            elif filtr[1] == '/': pos = 7 # filter statuses
            elif filtr[1] == '@': pos = 6 # filter tags
            elif filtr[1] == '%': pos = 5 # filter coverage
            else: pos = 1; shift = 0; endidx = -3 # filter test names
            for s in filteredTests:
                strg = s.split(separator)[pos].strip()
                if endidx != 0: strg = strg[:endidx]
                if len(filtr) < 3: continue
                try:
                    if filtr[1+shift] == '!':
                        if len(filtr) > 3 and re.search(filtr[2+shift:], strg) is not None:
                            insertedTests.discard(s)
                    elif re.search(filtr[1+shift:], strg) is None:
                        insertedTests.discard(s)
                except re.error: pass

    Listbox.delete(0, 'end')
    if filters:
        for s in sortTests(insertedTests): Listbox.insert('end', s)
    else:
        for s in sortTests(TESTS): Listbox.insert('end', s)
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)
    return True

#==============================================================================
# Ouvre un editeur sur le test
#==============================================================================
def viewTest(event=None):
    editor = PREFS.get("editor", "emacs")
    selection = Listbox.curselection()
    for s in selection:
        t = Listbox.get(s)
        splits = t.split(separator)
        module = splits[0]
        test = splits[1]
        module = module.strip()
        test = test.strip()
        modulesDir = MODULESDIR[BASE4COMPARE][module]
        if module == 'CFDBase':
            pathl = os.path.join(modulesDir, CFDBASEPATH, test)
            test = 'compute.py'
        else:
            pathl = os.path.join(modulesDir, module, 'test')
        if mySystem == 'mingw' or mySystem == 'windows':
            pathl = pathl.replace('/', '\\')
            cmd = f"cd {pathl} && {editor} {test}"
        else:
            cmd = f"cd {pathl}; {editor} {test}"
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

#==============================================================================
# Selectionne les tests affiche
# Met a jour les infos de progression
#==============================================================================
def selectAll(event=None):
    Listbox.selection_set(0, 'end')
    (total, remaining) = getTestsTime()
    displayProgress(0, total, remaining, 0.)

#==============================================================================
# Affiche les test FAILED ou MEMLEAK dans la listbox
#==============================================================================
def showFilter(filter='FAILED'):
    Listbox.delete(0, 'end')
    for s in TESTS:
        if re.search(filter, s) is not None:
            Listbox.insert('end', s)
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)
    Filter.set(''); TextFilter.update()
    return True

#==============================================================================
# Affiche les test qui ont deja tournes dans la listbox
#==============================================================================
def showRunCases():
    filter = r'\.\.\.'
    Listbox.delete(0, 'end')
    for s in TESTS:
        if re.search(filter, s) is None:
            Listbox.insert('end', s)
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)
    Filter.set(''); TextFilter.update()
    return True

#==============================================================================
# Affiche les test qui n'ont deja tournes dans la listbox
#==============================================================================
def showUnrunCases():
    filter = r'\.\.\.'
    Listbox.delete(0, 'end')
    for s in TESTS:
        if re.search(filter, s) is not None:
            Listbox.insert('end', s)
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)
    Filter.set(''); TextFilter.update()
    return True

#=============================================================================
# Affiche les tests couverts
#==============================================================================
def showCovered():
    filter = '100%'
    Listbox.delete(0, 'end')
    for s in TESTS:
        if re.search(filter, s) is not None:
            Listbox.insert('end', s)
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)
    Filter.set(''); TextFilter.update()
    return True

#=============================================================================
# Affiche les tests non couverts (0%)
#==============================================================================
def showUncovered():
    filter = ' 0%'
    Listbox.delete(0, 'end')
    for s in TESTS:
        if re.search(filter, s) is not None:
            Listbox.insert('end', s)
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)
    Filter.set(''); TextFilter.update()
    return True

#==============================================================================
# Affiche les tests partiellement couverts
#==============================================================================
def showPartialCovered():
    filter1 = '100%'
    filter2 = ' 0%'
    filter3 = r'\.%'
    Listbox.delete(0, 'end')
    for s in TESTS:
        if (re.search(filter1, s) is None and re.search(filter2, s) is None
                and re.search(filter3, s) is None):
            Listbox.insert('end', s)
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)
    Filter.set(''); TextFilter.update()
    return True

#==============================================================================
# Converti une chaine de type 0m0.12s (CPU time) en secondes
#==============================================================================
def string2Time(s):
    s = s.split('m')
    if len(s) != 2: return 0. # fail
    m1 = s[0]; s1 = s[1]; s1 = s1.replace('s', '')
    try: ret = float(m1)*60.+float(s1)
    except: ret = 0.
    return ret

#==============================================================================
# Converti un temps en secondes en une chaine 0h00m00s
#==============================================================================
def time2String(time):
    secs = time
    hours = int(secs / 3600)
    secs = secs - hours*3600
    mins = int(secs / 60)
    secs = secs - mins*60
    secs = int(secs)
    return "%1dh%2dm%2ds"%(hours,mins,secs)

#==============================================================================
# Affiche les tests plus rapide que ref CPUtime dans la listbox
#==============================================================================
def showFaster():
    Listbox.delete(0, 'end')
    for s in TESTS:
        s1 = s.split(separator)
        t1 = s1[2]; t2 = s1[3]
        t1 = string2Time(t1) # new time
        t2 = string2Time(t2)
        if t1 > 0 and t2 > 0:
            if t1 < t2-0.15*t2: Listbox.insert('end', s)
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)
    Filter.set(''); TextFilter.update()
    return True
def showFasterP():
    Listbox.delete(0, 'end')
    for s in TESTS:
        s1 = s.split(separator)
        t1 = s1[2]; t2 = s1[3]
        t1 = string2Time(t1) # new time
        t2 = string2Time(t2)
        if t1 > 0 and t2 > 0:
            if t1 < t2-0.5*t2: Listbox.insert('end', s)
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)
    Filter.set(''); TextFilter.update()
    return True

#==============================================================================
# Affiche les tests plus lent que la reference de 15%
#==============================================================================
def showSlower():
    Listbox.delete(0, 'end')
    for s in TESTS:
        s1 = s.split(separator)
        t1 = s1[2]; t2 = s1[3]
        t1 = string2Time(t1) # new time
        t2 = string2Time(t2)
        if t1 > 0 and t2 > 0:
            if (t1 > t2+0.15*t2): Listbox.insert('end', s)
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)
    Filter.set(''); TextFilter.update()
    return True
def showSlowerP():
    Listbox.delete(0, 'end')
    for s in TESTS:
        s1 = s.split(separator)
        t1 = s1[2]; t2 = s1[3]
        t1 = string2Time(t1) # new time
        t2 = string2Time(t2)
        if t1 > 0 and t2 > 0:
            if t1 > t2+0.5*t2: Listbox.insert('end', s)
    Listbox.config(yscrollcommand=Scrollbar.set)
    Scrollbar.config(command=Listbox.yview)
    Filter.set(''); TextFilter.update()
    return True


#==============================================================================
# Affiche tous les tests
#==============================================================================
def showAll():
    Filter.set(''); TextFilter.update()
    filterTestList()

#==============================================================================
# Stop l'execution des tests
#==============================================================================
def stopTests():
    global STOP, THREAD, PROCESS
    STOP = 1

    if PROCESS is not None:
        if mySystem == 'mingw' or mySystem == 'windows':
            subprocess.call(['taskkill', '/F', '/T', '/PID', str(PROCESS.pid)])
        else: # unix
            # try soft, then hard
            os.killpg(os.getpgid(PROCESS.pid), signal.SIGTERM)
            os.kill(PROCESS.pid, signal.SIGTERM)
            os.killpg(os.getpgid(PROCESS.pid), signal.SIGKILL)
            os.kill(PROCESS.pid, signal.SIGKILL)

        PROCESS = None
        displayStatus(0)

    if THREAD is not None:
        print("Info: stopping thread...")
        #THREAD._stop() # kill?
        #THREAD.join() # wait
        #THREAD.terminate()
        THREAD = None
        displayStatus(0)

#==============================================================================
# Affiche le status: running/stopped
#==============================================================================
def displayStatus(status):
    if status == 0: Status.set('Stopped'); StatusLabel.config(bg='red')
    else: Status.set('Running'); StatusLabel.config(bg='green')
    StatusLabel.update()

#==============================================================================
# Affiche la progression
# IN: current: le nbre de tests restants
# IN: total: nbre total de tests
# IN: remaining: temps restant avant la fin
# IN: elapsed: temps passe
#==============================================================================
def displayProgress(current, total, remaining, elapsed):
    Progress.set("%3d/%3d [%s/%s]"%
                 (current,total,time2String(remaining),time2String(elapsed)))
    ProgressLabel.update()

#==============================================================================
# Modifie le nbre de threads utilises pour la valid
#==============================================================================
def setThreads(event=None):
    nt = Threads.get()
    try:
        nti = int(nt)
        KCore.kcore.setOmpMaxThreads(nti)
        print('Info: Num threads set to %d.\n'%nti)
    except:
        print('Info: Bad thread number.\n')

#==============================================================================
# Recupere le nbre de threads (OMP_NUM_THREADS)
#==============================================================================
def getThreads():
    nt = KCore.kcore.getOmpMaxThreads()
    Threads.set(str(nt))
    TextThreads.update()

#==============================================================================
# Exporte les resultats de la valid dans un fichier texte
#==============================================================================
def export2Text():
    import tkinter.filedialog as tkFileDialog
    ret = tkFileDialog.asksaveasfilename()
    if ret == '' or ret is None or ret == (): # user cancel
        return
    with open(ret, 'w') as f:
        for t in TESTS: f.write(t); f.write('\n')

#==============================================================================
# writeSessionLog: write log and baseTime
#==============================================================================
def createEmptySessionLog():
    # Create an empty session log
    with open(os.path.join(VALIDDIR['LOCAL'], "session.log"), "w") as f:
        f.write("")

def writeSessionLog():
    cassiopeeIncDir = getInstallPaths()[0]
    gitOrigin = Dist.getGitOrigin(cassiopeeIncDir)
    gitBranch = Dist.getGitBranch(cassiopeeIncDir)
    gitHash = Dist.getGitHash(cassiopeeIncDir)[:7]
    gitInfo = "Git origin: {}\nGit branch: {}, commit hash: {}\n".format(
        gitOrigin, gitBranch, gitHash)

    messageText = "Base from {}\n{}\n".format(cassiopeeIncDir, gitInfo)
    for t in TESTS: messageText += t+'\n'

    # Write time stamp dans ValidData/base.time et
    # log dans ValidData/session.log
    writeFinal(os.path.join(VALIDDIR['LOCAL'], 'base.time'), gitInfo=gitInfo)
    writeFinal(os.path.join(VALIDDIR['LOCAL'], 'session.log'), logTxt=messageText)

#==============================================================================
# Notify "Commit ready"
#==============================================================================
def notifyValidOK():
    cassiopeeIncDir = getInstallPaths()[0]
    gitOrigin = Dist.getGitOrigin(cassiopeeIncDir)
    gitBranch = Dist.getGitBranch(cassiopeeIncDir)
    gitHash = Dist.getGitHash(cassiopeeIncDir)[:7]
    gitInfo = "Git origin: {}\nGit branch: {}, commit hash: {}\n".format(
        gitOrigin, gitBranch, gitHash)

    messageText = "Base from {}\n{}\n".format(cassiopeeIncDir, gitInfo)
    for t in TESTS: messageText += t+'\n'

    try:
        from KCore.notify import notify
        notify(recipients=['christophe.benoit@onera.fr'],
               messageSubject='[Cassiopee] Ready to commit',
               messageText=messageText)
    except ImportError:
        print("Error: KCore is required to import notify.")
        sys.exit()

#==============================================================================
def Quit(event=None, sessionName="session"):
    import os
    import shutil
    logname = os.path.join(VALIDDIR['LOCAL'], "session.log")
    # The session log is copied if it is not empty and if we have write
    # permissions
    if os.access(VALIDDIR['LOCAL'], os.W_OK) and (not os.path.getsize(logname) == 0):
        now = time.strftime("%y%m%d_%H%M%S", time.localtime())
        dst = os.path.join(VALIDDIR['LOCAL'], "{}-{}.log".format(sessionName, now))
        print("Saving session to: {}".format(dst))
        shutil.copyfile(logname, dst)
    # Write test metadata
    writeTestMetadata()
    os._exit(0)

#==============================================================================
# Ajoute une etoile a la selection. Tagger plusieurs fois une selection permet
# de changer de symbole: *, r, g, b
# Les tags sont locaux, cad, propres a l'utilisateur meme quand la base choisie
# est globale
#==============================================================================
def tagSelection(event=None):
    global TESTS, TESTMETA, TESTMETA_UPDATE
    TESTMETA_UPDATE = True
    tagSymbols = '* r g b'.split()
    ntags = len(tagSymbols)
    selection = Listbox.curselection()
    for s in selection:
        no = int(s)
        t = Listbox.get(s)
        splits = t.split(separator)
        module = splits[0].strip()
        test = splits[1].strip()
        testName = module+'/'+test
        tag = splits[6].strip()
        if not tag: tag = '*'
        else: tag = tagSymbols[(tagSymbols.index(tag)+1)%ntags]
        if testName in TESTMETA:
            TESTMETA[testName]['tag'] = tag
        splits[6] = ' {} '.format(tag)
        s = separator.join(i for i in splits)
        regTest = re.compile(' '+test+' ')
        regModule = re.compile(module+' ')
        for c, tt in enumerate(TESTS):
            if regModule.search(tt) is not None:
                if regTest.search(tt) is not None: TESTS[c] = s; break
        Listbox.delete(no, no)
        Listbox.insert(no, s)
        Listbox.selection_set(no)
    return

def untagSelection(event=None):
    global TESTS, TESTMETA, TESTMETA_UPDATE
    TESTMETA_UPDATE = True
    selection = Listbox.curselection()
    for s in selection:
        no = int(s)
        t = Listbox.get(s)
        splits = t.split(separator)
        module = splits[0].strip()
        test = splits[1].strip()
        testName = module+'/'+test
        if testName in TESTMETA:
            TESTMETA[testName]['tag'] = ' '
        splits[6] = ' '*3
        s = separator.join(i for i in splits)
        regTest = re.compile(' '+test+' ')
        regModule = re.compile(module+' ')
        for c, tt in enumerate(TESTS):
            if regModule.search(tt) is not None:
                if regTest.search(tt) is not None: TESTS[c] = s; break
        Listbox.delete(no, no)
        Listbox.insert(no, s)
        Listbox.selection_set(no)
    return

#==============================================================================
# Setups to either use the local or global databases
#==============================================================================
def toggleDB(**kwargs):
    if BASE4COMPARE == 'LOCAL': return setupGlobal(**kwargs)
    return setupLocal(**kwargs)

def setupLocal(**kwargs):
    global BASE4COMPARE
    # Change to local ref
    print('Info: comparing to local database.')
    BASE4COMPARE = 'LOCAL'
    os.environ['VALIDLOCAL'] = '.'
    casFolder = os.path.join(os.getenv('CASSIOPEE'), "Cassiopee")
    casValidFolder = os.path.join(casFolder, "Valid{}".format(DATA))
    # Create ValidData folder if not present and permissions are OK
    if os.access(casFolder, os.W_OK):
        os.makedirs(casValidFolder, exist_ok=True)
    else:
        os.environ['VALIDLOCAL'] = os.path.join(os.getcwd(), "Valid{}".format(DATA))

    if INTERACTIVE: WIDGETS['UpdateButton'].configure(state=TK.NORMAL)
    createEmptySessionLog()
    loadTestMetadata()
    buildTestList(**kwargs)
    updateDBLabel()
    setGUITitleBar(loc='LOCAL')
    return 0

def setupGlobal(**kwargs):
    global BASE4COMPARE
    if VALIDDIR['GLOBAL'] is None: return 1
    # Change to global ref
    print('Info: comparing to global database.')
    BASE4COMPARE = 'GLOBAL'
    os.environ['VALIDLOCAL'] = VALIDDIR['LOCAL']

    # No update on global ref
    if INTERACTIVE: WIDGETS['UpdateButton'].configure(state=TK.DISABLED)
    # Change to match the numthreads of global ref
    try:
        with open(os.path.join(VALIDDIR['GLOBAL'], 'base.time'), 'r') as f:
            db_info = f.read().split('\n')
            Threads.set(db_info[2])
            setThreads()
    except: pass
    createEmptySessionLog()
    loadTestMetadata()
    buildTestList(**kwargs)
    updateDBLabel()
    setGUITitleBar(loc='GLOBAL')
    return 0

def getDBInfo():
    dbInfo = ''
    if os.access('/stck/cassiope/git/Cassiopee/', os.R_OK):
        filename = '/stck/cassiope/git/Cassiopee/Cassiopee/Valid{}/base.time'
        try:
            with open(filename.format(DATA), 'r') as f:
                dbInfo = f.read().split('\n')
                dbInfo = '[{} - {} - {} threads]'.format(*dbInfo[0:3])
        except: pass
    return dbInfo

def updateDBLabel():
    if not INTERACTIVE: return
    dbInfo = getDBInfo()
    if not dbInfo: return
    if BASE4COMPARE == 'GLOBAL':
        label = 'Switch to local data base'
    else:
        label = 'Switch to global data base ' + dbInfo
    toolsTab.entryconfig(3, label=label)

#==============================================================================
# Parse command-line arguments
#==============================================================================
def parseArgs():
    import argparse
    def _checkInt(n):
        def _throwError():
            raise argparse.ArgumentTypeError("Number of remaining logs must be "
                                             "a positive integer")
            sys.exit()
        try: n = int(n)
        except: _throwError()
        if n > 0: return n
        else: _throwError()

    # Create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filters", type=str, default='',
                        help="Single-quoted test filters")
    parser.add_argument("-gdb", "--global-database", action="store_true",
                        dest="global_db",
                        help="Switch to global database. Default: local database")
    parser.add_argument("-l", "--load-session", dest='loadSession',
                        action="store_true",
                        help="Load last session. Default: False")
    parser.add_argument("--leak-sanitizer", action="store_true",
                        dest="leak_sanitizer",
                        help="Leak Sanitizer to detect memory leaks."
                             "Available in DEBUG mode only")
    parser.add_argument("--memory-sanitizer", action="store_true",
                        dest="memory_sanitizer",
                        help="Address Sanitizer to detect memory access errors."
                             "Available in DEBUG mode only")
    parser.add_argument("-p", "--purge", default=50, type=_checkInt,
                        help="Purge session logs down to the last X. Default: 50")
    parser.add_argument("-r", "--run", action="store_true",
                        help="Run selected tests")
    parser.add_argument("-s", "--session-name", type=str, default='session',
                        dest="sessionName",
                        help="Name of the session file. Default: session")
    parser.add_argument("--update", action="store_true",
                        help="Update local database")

    # Parse arguments
    return parser.parse_args()

# Purge session logs by date down to the last n most recent
def purgeSessionLogs(n, sessionName="session"):
    lognames = sorted(glob.glob(os.path.join(VALIDDIR['LOCAL'], '{}-*.log'.format(sessionName))))
    if len(lognames) > n:
        for log in lognames[:-n]: os.remove(log)
    return None

#==============================================================================
# Switches to control the state of USE_ASAN (Address/Leak Sanitizers)
#==============================================================================
def toggleASAN():
    global USE_ASAN
    USE_ASAN[0] = not USE_ASAN[0]
    updateASANLabel(3)

def toggleLSAN():
    global USE_ASAN
    USE_ASAN[1] = not USE_ASAN[1]
    updateASANLabel(4)
    updateASANOptions()

def updateASANOptions():
    # Update ASAN_OPTIONS according to USE_ASAN
    asan_opt = os.getenv("ASAN_OPTIONS", "")
    posArg = asan_opt.find('detect_leaks=')
    if posArg == -1:
        os.environ["ASAN_OPTIONS"] = "{}:detect_leaks={:d}".format(
            asan_opt, USE_ASAN[1])
    else:
        os.environ["ASAN_OPTIONS"] = "{}detect_leaks={:d}{}".format(
            asan_opt[:posArg], USE_ASAN[1], asan_opt[posArg+14:])
    print("Info: ASAN_OPTIONS = " + os.getenv("ASAN_OPTIONS", ""))
    print("      LSAN_OPTIONS = " + os.getenv("LSAN_OPTIONS", ""))

def updateASANLabel(entry_index):
    if not INTERACTIVE: return
    if getDBInfo(): entry_index += 2
    label = toolsTab.entrycget(entry_index, "label")
    if label[0].startswith('E'): label = 'Dis' + label[2:]
    else: label = 'En' + label[3:]
    toolsTab.entryconfig(entry_index, label=label)

#==============================================================================
# Set message in the title bar
#==============================================================================
def setGUITitleBar(loc='GLOBAL'):
    if not INTERACTIVE: return
    # Machine name
    machine = socket.gethostname()
    title = '*Cassiopee* valid {} : {} @ {}'.format(loc, os.getenv("ELSAPROD"),
                                                    machine)
    cassiopeeIncDir = getInstallPaths()[0]
    gitBranch = Dist.getGitBranch(cassiopeeIncDir)
    if gitBranch and gitBranch != "main":
        title += " (branch {})".format(gitBranch)
    Master.title(title)

#==============================================================================
# Load a dict storing user settings for validCassiopee and amend the PREFS dict
# Located in $HOME/.cassiopee/config_validCassiopee.json
#==============================================================================
def loadPrefFile():
    import json
    global PREFS
    homePath = os.path.expanduser('~')
    if homePath is None: homePath = os.getenv('HOME')
    if homePath is None: homePath = os.getenv('USERPROFILE')
    if homePath is None: homePath = ''
    kdir = os.path.join(homePath, '.cassiopee')
    kdirExist = os.path.exists(kdir) and os.path.isdir(kdir)
    try:
        if not kdirExist: os.makedir(kdir)
    except: return
    configFile = os.path.join(kdir, 'config_validCassiopee.json')
    if os.path.isfile(configFile) and os.access(configFile, os.R_OK):
        with open(configFile, 'r') as f:
            localPrefs = json.load(f)
        for k, v in localPrefs.items(): PREFS[k] = v
    else:
        savePrefFile()

#==============================================================================
# Save user settings to $HOME/.cassiopee/config_validCassiopee.json
#==============================================================================
def savePrefFile():
    import json
    homePath = os.path.expanduser('~')
    if homePath is None: homePath = os.getenv('HOME')
    if homePath is None: homePath = os.getenv('USERPROFILE')
    if homePath is None: homePath = ''
    kdir = os.path.join(homePath, '.cassiopee')
    try: os.makedirs(kdir, exist_ok=True)
    except OSError as e: return
    if not os.access(kdir, os.W_OK): return
    configFile = os.path.join(kdir, 'config_validCassiopee.json')
    with open(configFile, 'w') as f:
        json.dump(PREFS, f, indent=4)

#==============================================================================
# Main
#==============================================================================

if __name__ == '__main__':
    # Check that CASSIOPEE and the install paths are consistent
    checkEnvironment()
    # Get name of the ValidData folder
    DATA = Dist.getDataFolderName()
    # Set MODULESDIR and VALIDDIR paths once, both locally and globally
    setPaths()

    if INTERACTIVE:
        # --- Use GUI ---
        import tkinter as TK
        import CPlot.Tk as CTK
        import tkinter.font as Font
        from functools import partial
        # Load user settings
        loadPrefFile()

        # Main window
        Master = TK.Tk()
        setGUITitleBar()
        Master.columnconfigure(0, weight=1)
        Master.rowconfigure(0, weight=1)
        #GENERALFONT = ('Courier', 9)
        GENERALFONT = ('Andale Mono', 9)

        Master.option_add('*Font', GENERALFONT)
        generalFont = Font.Font(family=GENERALFONT[0], size=GENERALFONT[1])
        generalFontS = generalFont.measure(' ')*1.
        generalFontA = generalFont.measure('a')*1.
        generalFontFixed = generalFont.metrics('fixed')

        # Main menu
        menu = TK.Menu(Master)
        fileTab = TK.Menu(menu, tearoff=0)
        menu.add_cascade(label='File', menu=fileTab)
        toolsTab = TK.Menu(menu, tearoff=0)
        menu.add_cascade(label='Tools', menu=toolsTab)
        sortTab = TK.Menu(menu, tearoff=0)
        menu.add_cascade(label='Sort', menu=sortTab)
        viewTab = TK.Menu(menu, tearoff=0)
        menu.add_cascade(label='View', menu=viewTab)

        loadSessionWithArgs = partial(buildTestList, "session")
        fileTab.add_command(label='Load last session', command=loadSessionWithArgs)
        fileTab.add_command(label='Purge session', command=buildTestList)
        fileTab.add_command(label='Export to text file', command=export2Text)
        #fileTab.add_command(label='Notify Ready for commit', command=notifyValidOK)
        fileTab.add_command(label='Quit', command=Quit, accelerator='Ctrl+Q')
        viewTab.add_command(label='Show FAILED', command=showFilter)
        viewTab.add_command(label='Show MEMLEAK',
                            command=partial(showFilter, "MEMLEAK"))
        viewTab.add_command(label='Show TIMEOUT',
                            command=partial(showFilter, "TIMEOUT"))
        viewTab.add_command(label='Show ALL tests', command=showAll)
        viewTab.add_separator()
        viewTab.add_command(label='Show run cases', command=showRunCases)
        viewTab.add_command(label='Show unrun cases', command=showUnrunCases)
        viewTab.add_separator()
        viewTab.add_command(label='Show covered (100%)', command=showCovered)
        viewTab.add_command(label='Show partially covered (x%)',
                            command=showPartialCovered)
        viewTab.add_command(label='Show uncovered (0%)', command=showUncovered)
        viewTab.add_separator()
        viewTab.add_command(label='Show faster (-15%)', command=showFaster)
        viewTab.add_command(label='Show slower (+15%)', command=showSlower)
        viewTab.add_command(label='Show faster (-50%)', command=showFasterP)
        viewTab.add_command(label='Show slower (+50%)', command=showSlowerP)
        viewTab.add_separator()
        viewTab.add_command(label='Select all visible tests', command=selectAll,
                            accelerator='Ctrl+A')

        sortTab.add_command(label='(*) By Name', command=lambda: sortTestList(entry="Name"))
        sortTab.add_command(label='    By CPU time', command=lambda: sortTestList(entry="CPU time"))
        sortTab.add_command(label='    By Ref. CPU time', command=lambda: sortTestList(entry="Ref. CPU time"))
        sortTab.add_command(label='    By Date', command=lambda: sortTestList(entry="Date"))
        sortTab.add_command(label='    By Coverage', command=lambda: sortTestList(entry="Coverage"))
        sortTab.add_command(label='    By Tag', command=lambda: sortTestList(entry="Tag"))
        sortTab.add_command(label='    By Status', command=lambda: sortTestList(entry="Status"))
        sortTab.add_command(label='    By CPU time relDiff.', command=lambda: sortTestList(entry="CPU time relDiff."))
        sortTab.add_separator()
        sortTab.add_command(label='Reverse order', command=updateSortingOrderLabel)

        toolsTab.add_command(label='Tag selection', command=tagSelection)
        toolsTab.add_command(label='Untag selection', command=untagSelection)

        dbInfo = getDBInfo()
        if dbInfo:
            # Show this button if the global database can be interrogated
            toolsTab.add_separator()
            toolsTab.add_command(label='Switch to global data base ' + dbInfo,
                                 command=toggleDB)
        if Dist.DEBUG and os.getenv('ASAN_LIB') is not None:
            toolsTab.add_separator()
            toolsTab.add_command(label='Enable Address Sanitizer (ASan)',
                                 command=toggleASAN)
            toolsTab.add_command(label='Enable Leak Sanitizer (LSan)',
                                 command=toggleLSAN)
            updateASANOptions()

        Master.config(menu=menu)
        Master.bind_all("<Control-q>", Quit)
        Master.protocol("WM_DELETE_WINDOW", Quit)
        Master.bind_all("<Control-a>", selectAll)

        # Main frame
        Frame = TK.Frame(Master)
        Frame.columnconfigure(0, weight=1)
        Frame.rowconfigure(0, weight=1)
        Frame.columnconfigure(1, weight=1)
        Frame.grid(row=0, column=0, sticky=TK.EW)

        Listbox = TK.Listbox(Frame, selectmode=TK.EXTENDED, width=125,
                             height=39, background='White')
        Listbox.grid(row=0, column=0, columnspan=11, sticky=TK.NSEW)

        Scrollbar = TK.Scrollbar(Frame, orient=TK.VERTICAL)
        Scrollbar.grid(row=0, column=11, sticky=TK.NSEW)

        Status = TK.StringVar(Master)
        StatusLabel = TK.Label(Frame, textvariable=Status)
        Status.set('Stopped'); StatusLabel.config(bg='red')
        StatusLabel.grid(row=1, column=0, sticky=TK.EW)

        Progress = TK.StringVar(Master)
        ProgressLabel = TK.Label(Frame, textvariable=Progress)
        Progress.set('  0/  0 [0h 0m 0s/0h 0m 0s]')
        ProgressLabel.grid(row=1, column=1, sticky=TK.EW)

        Filter = TK.StringVar(Master)
        TextFilter = TK.Entry(Frame, textvariable=Filter, background='White',
                              width=50)
        TextFilter.bind('<KeyRelease>', filterTestList)
        TextFilter.grid(row=1, column=2, columnspan=3, sticky=TK.EW)

        filterInfoBulle = 'Filter test database using a regexp.\n'+'-'*70+'\n'\
            '1) White-spaced: ^cylinder ^sphere\n'\
            '2) Module filter using #: #Apps #Fast #FF   or simply   #[A,F] \n'\
            '3) Status filter using /: /FAILED /MEMLEAK   or simply   /F\n'\
            '4) Coverage filter using %: %100\n'\
            '5) Tag symbol filter using @: @r   to catch red-coloured cases\n'\
            '6) Keyworded filters: <SEQ>, <DIST>, <RUN>, <UNRUN>, <TAG>, <UNTAG>.\n'\
            '7) Logical OR ops unless prefixed with & (AND): #Converter &/FAILED\n'\
            '8) Negated using !: #Fast &#!FastC (innermost symbol)'

        RunButton = TK.Button(Frame, text='Run', command=runTestsInThread,
                              fg='blue')
        RunButton.grid(row=1, column=5, sticky=TK.EW)

        Button = TK.Button(Frame, text='Stop', command=stopTests, fg='red')
        Button.grid(row=1, column=6, sticky=TK.EW)
        updateTestsInThread = partial(runTestsInThread, True)
        UpdateButton = TK.Button(Frame, text='Update',
                                 command=updateTestsInThread, fg='blue')
        WIDGETS['UpdateButton'] = UpdateButton
        UpdateButton.grid(row=1, column=7, sticky=TK.EW)
        Button = TK.Button(Frame, text='Edit', command=viewTest)
        Button.grid(row=1, column=8, sticky=TK.EW)

        Threads = TK.StringVar(Master)
        TextThreads = TK.Entry(Frame, textvariable=Threads, background='White',
                               width=3)
        TextThreads.grid(row=1, column=9, columnspan=2, sticky=TK.EW)
        TextThreads.bind('<Return>', setThreads)
        TextThreads.bind('<KP_Enter>', setThreads)
        getThreads()

        Frame.grid(sticky=TK.NSEW)

        CTK.infoBulle(parent=TextFilter, text=filterInfoBulle)
        CTK.infoBulle(parent=RunButton, text='Run selected tests.')
        CTK.infoBulle(parent=UpdateButton,
                      text='Update tests (replace data base files).')
        CTK.infoBulle(parent=TextThreads, text='Number of threads.')
        if isDBAdmin(): setupLocal()  # Local valid for the DB admin 'cassiope'
        else:
            ierr = setupGlobal()  # Comparison is made against the global valid
            if ierr == 1: setupLocal()  # Global valid does not exist, default back to local
        TK.mainloop()
    else:
        # --- Command line execution ---
        vcargs = parseArgs()

        generalFontFixed = 1
        Listbox = NoDisplayListbox()
        Scrollbar = NoDisplayScrollbar()
        Status = NoDisplayStringVar()
        StatusLabel = NoDisplayLabel()
        Progress = NoDisplayStringVar()
        ProgressLabel = NoDisplayLabel()
        Filter = NoDisplayStringVar()
        TextFilter = NoDisplayEntry()
        WIDGETS['UpdateButton'] = NoDisplayButton()
        Threads = NoDisplayStringVar()
        TextThreads = NoDisplayEntry()
        getThreads()

        sessionName = vcargs.sessionName if vcargs.loadSession else None
        if (os.access('/stck/cassiope/git/Cassiopee/', os.R_OK) and
                vcargs.global_db and not (vcargs.update or isDBAdmin())):
            ierr = setupGlobal(sessionName=sessionName)
            if ierr == 1: setupLocal()  # Global valid does not exist, default back to local
        else: setupLocal(sessionName=sessionName)
        purgeSessionLogs(n=vcargs.purge, sessionName=vcargs.sessionName)
        if vcargs.filters:
            Filter.set(vcargs.filters)
            filterTestList()
        if vcargs.run:
            if Dist.DEBUG and os.getenv('ASAN_LIB') is not None:
                if vcargs.memory_sanitizer: USE_ASAN[0] = True
                if vcargs.leak_sanitizer: USE_ASAN[1] = True
                updateASANOptions()
            selectAll()
            runTests(update=vcargs.update)
            Quit(sessionName=vcargs.sessionName)
