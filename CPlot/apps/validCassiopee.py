# *Cassiopee* GUI for validation and tests
try: import Tkinter as TK
except: import tkinter as TK
try: import tkFont as Font
except: import tkinter.font as Font
import os, sys, re, signal, platform
from functools import partial
import numpy as np
import subprocess 
import threading
import time
import KCore
import KCore.Dist as Dist
import CPlot.Tk as CTK

# CASSIOPEE var
# doit etre le chemin des sources avec les tests unitaires
CASSIOPEE = os.getenv('CASSIOPEE_SOURCES')
if CASSIOPEE is None or CASSIOPEE == '':
    CASSIOPEE = os.getenv('CASSIOPEE')
    if CASSIOPEE is None or CASSIOPEE == '':
        print('Error: CASSIOPEE must be present in your environment.')
        sys.exit()

# CFD Base
CFDBASEPATH = '/Validation/Cases'

# System
mySystem = Dist.getSystem()[0]

# Machine name
import socket
machine = socket.gethostname()

# Support MPI?
try: 
    import mpi4py
    isMpi = True
except: isMpi = False

# Check svn version
CHECKSVNVERSION = True

# Regexprs
regDiff = re.compile('DIFF')
regFailed = re.compile('FAILED')
regError = re.compile('Error')
regErreur = re.compile('Erreur') # because of french system
regAbort = re.compile('Aborted')
regSegmentation = re.compile('Segmentation')
separator = ':'
separatorl = separator+' '
expTest1 = re.compile("_t[0-9]+") # normal tests
expTest2 = re.compile(".py")
expTest3 = re.compile("\~")
expTest4 = re.compile("_m[0-9]+") # distributed tests

# Liste des tous les tests obtenue en listant les repertoires
# Un element de la liste est une string comme affichee dans la listbox
TESTS = []
# Test filter: 0 (no filter), 1: sequential tests only, 2: distributed tests only
TESTS_FILTER = 0
# Repertoire 'module' des modules
MODULESDIR = {}

# Si THREAD est None, les test unitaires ne tournent pas
# Sinon, THREAD vaut le thread lance
THREAD = None
# Le process lance sinon None
PROCESS = None

# Est egal a 1 si on doit s'arreter
STOP = 0

# WIDGETS dict
WIDGETS = {}

# Cree sessionLog et le vide
if not os.path.exists(CASSIOPEE+'/Apps/Modules/ValidData'):
    os.mkdir(CASSIOPEE+'/Apps/Modules/ValidData')
    
f = open(CASSIOPEE+"/Apps/Modules/ValidData/session.log", "w")
f.write("")
f.close()

#==============================================================================
# Simulate check_output since it doesn't existe for early version of python
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
        PROCESS = subprocess.run(cmd, check=True, shell=shell, stderr=stderr, stdout=subprocess.PIPE)
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
                                   stderr=subprocess.PIPE, cwd=wdir, shell=shell, preexec_fn=ossid)
        
        # max accepted time is 2 minutes for one repetition of a test
        nreps = Repeats.get()
        stdout, stderr = PROCESS.communicate(None, timeout=60.*2.*nreps)
        
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
    if module == 'CFDBase':
        path = CASSIOPEE+CFDBASEPATH
        fileTime = '%s/%s/Data/%s.time'%(path, test, test)
        fileStar = '%s/%s/Data/%s.star'%(path, test, test)
    else:
        modulesDir = MODULESDIR[module]
        path = CASSIOPEE+'/Apps/'+modulesDir
        testr = os.path.splitext(test)
        fileTime = '%s/%s/test/Data/%s.time'%(path, module, testr[0])
        fileStar = '%s/%s/test/Data/%s.star'%(path, module, testr[0])
    a = os.access(fileTime, os.F_OK)
    if a:
        f = open(fileTime, 'r')
        list = f.read()
        f.close()
        list = list.split('\n')
        if len(list) > 0: refDate = list[0]
        else: refDate = '...'
        if len(list) > 1: refCPUtime = list[1]
        else: refCPUtime = '...'
        if len(list) > 3 and list[2] != '': refCoverage = list[2]
        else: refCoverage = '...%'
    else:
        refDate = '...'
        refCPUtime = '...'
        refCoverage = '...%'

    a = os.access(fileStar, os.F_OK)
    if a:
        f = open(fileStar, 'r')
        list = f.read()
        f.close()
        list = list.split('\n')
        if (len(list) > 0): refTag = list[0]
        else: refTag = ' '
    else: refTag = ' '

    execTime = '../../.. ..h..'
    if status != '...': # Not First call
        execTime = time.strftime('%d/%m/%y %Hh%M', time.localtime())

    if coverage == '...%': coverage = refCoverage
    if tag == ' ': tag = refTag

    s = ljust(module, 13)+separatorl+ljust(test, 40)+separatorl+\
        ljust(CPUtime, 10)+separatorl+ljust(refCPUtime, 10)+separatorl+\
        ljust(refDate, 16)+separatorl+ljust(coverage, 5)+separatorl+\
        tag.ljust(2)+separatorl+' '+ljust(status, 10)
    return s

#==============================================================================
# Retourne la liste des modules situes dans Apps/Modules et dans Apps/PModules
# Eventuellement peut ajouter "CFDBase", nom referencant les tests
# de validation des solveurs (CFDBase)
#==============================================================================
def getModules():
    # Tests unitaires des modules
    print('Info: Getting tests in:%s.'%CASSIOPEE)
    modules = []

    path = CASSIOPEE+'/Apps/PModules'
    try: mods = os.listdir(path)
    except: mods = []
    notTested = ['Upmost', 'FastP']
    for i in mods:
        if i not in notTested and i not in modules:
            a = os.access('%s/%s/test'%(path,i), os.F_OK)
            if a:
                modules.append(i)
                MODULESDIR[i] = 'PModules'

    path = CASSIOPEE+'/Apps/Modules'
    try: mods = os.listdir(path)
    except: mods = []
    for i in mods:
        if i not in modules:
            a = os.access('%s/%s/test'%(path,i), os.F_OK)
            if a: 
                modules.append(i)
                MODULESDIR[i] = 'Modules'
    
    # Validation CFD
    modules.append('CFDBase')
    MODULESDIR['CFDBase'] = ''
    return sorted(modules)

#==============================================================================
# Retourne la liste des tests unitaires d'un module
# si module == 'CFDBase', retourne la liste des cas de validation CFD
#==============================================================================
def getTests(module):
    a = []
    if module == 'CFDBase': a += getCFDBaseTests()
    else: a += getUnitaryTests(module)
    return a

#==============================================================================
# Retourne la liste des tests unitaires pour un module donne
# Les tests unitaires doivent etre dans module/test
# La variable globale TESTS_FILTER permet de filtrer tests sequentiels et tests
# distribues
#==============================================================================
def getUnitaryTests(module):
    modulesDir = MODULESDIR[module]
    path = '%s/Apps/%s/%s/test'%(CASSIOPEE, modulesDir, module)
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
    path = CASSIOPEE+CFDBASEPATH
    try: reps = os.listdir(path)
    except: reps = []
    out = []
    for r in reps: # a terme a supprimer
        if r == 'NACA': out.append(r) # MB 2D Euler
        elif r == 'NACA_IBC': out.append(r) # IBC 2D Euler
        elif r == 'DAUPHIN': out.append(r) # MB 3D Euler
        elif r == 'FLATPLATE': out.append(r) # MB 3D SA
        elif r == 'RAE2822': out.append(r) # MB 2D SA
        elif r == 'RAE2822_IBC': out.append(r) # IBC 2D SA
        elif r == 'CUBE_IBC': out.append(r) # IBC 3D SA
    return sorted(out)

#==============================================================================
# Ecrit un fichier contenant date, CPUtime, coverage
#==============================================================================
def writeTime(file, CPUtime, coverage):
    try:
        execTime = time.strftime('%d/%m/%y %Hh%M',time.localtime())
        f = open(file, 'w')
        f.write(execTime+'\n')
        f.write(CPUtime+'\n')
        f.write(coverage+'\n')
        f.close()
    except: pass

#==============================================================================
# Ecrit un fichier contenant date, machine, nbre de threads, svnVersion 
# et logTxt
#==============================================================================
def writeFinal(file, svnVersion=None, logTxt=None, append=False):
    execTime = time.strftime('%d/%m/%y %Hh%M', time.localtime())
    machine = platform.uname()
    if len(machine) > 1: machine = machine[1]
    else: machine = 'Unkwown'
    nthreads = Threads.get()
    mode = 'w'
    if append: mode = 'a'
    f = open(file, 'w')
    f.write(execTime+'\n')
    f.write(machine+'\n')
    f.write(nthreads+'\n')
    if svnVersion is not None: f.write(svnVersion+'\n')
    if logTxt is not None: f.write(logTxt+'\n')
    f.close()

#==============================================================================
# Read and update star dans un fichier star
#==============================================================================
def readStar(file):
    star = ' '
    try:
        f = open(file, 'r')
        star = f.readline().rstrip('\n')
        f.close()
    except: pass
    return star
    
def writeStar(file, star):
    try:
        f = open(file, 'w')
        f.write(star+'\n')
        f.close()
    except: pass

#==============================================================================
# Lance un seul test unitaire ou un cas de la base de validation
#==============================================================================
def runSingleTest(no, module, test):
    if module == 'CFDBase': return runSingleCFDTest(no, module, test)
    else: return runSingleUnitaryTest(no, module, test)

#==============================================================================
# extrait le temps CPU d'un chaine output (utile sur windows)
# retourne le temps CPU comme une chaine
# moyenne si plusieurs repetitions d'un meme cas unitaire
#==============================================================================
def extractCPUTime(output1, output2, nreps=1):
    CPUtime = 'Unknown'
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
    tf = (t2-t1)/float(nreps)
    hf = int(tf/3600.)
    tf = tf - 3600*hf
    mf = int(tf/60.)
    tf = tf - 60*mf
    sf = int(tf*100)/100.
    if hf > 0: CPUtime = '%dh%dm%gs'%(hf,mf,sf)
    else: CPUtime = '%dm%gs'%(mf,sf)
    return CPUtime

#=============================================================================
# Extrait le temps CPU d'une sortie time -p (unix)
# Moyenne si plusieurs repetitions d'un meme cas unitaire
#=============================================================================
def extractCPUTime2(output, nreps=1):
    i1 = output.find('real')
    output = output[i1+4:]
    output = output.replace(',', '.')
    output = output.lstrip()
    i2 = output.find(' ')
    if i2 != -1: output = output[:i2]
    i3 = output.find('\n')
    if i3 != -1: output = output[:i3]
    i1 = output.find('h')
    hf = 0; mf = 0; sf = 0
    if i1 != -1:
        hf = output[:i1]
        try: hf = float(hf)
        except: hf = 0.
        output = output[i1+1:]
    i1 = output.find('m')
    if i1 != -1:
        mf = output[:i1]
        try: mf = float(mf)
        except: mf = 0.
        output = output[i1+1:]
    sf = output.replace('s', '')
    try: sf = float(sf)
    except: sf = 0.
    tf = (hf*3600.+mf*60.+sf)/float(nreps)
    hf = int(tf/3600.)
    tf = tf - 3600*hf
    mf = int(tf/60.)
    tf = tf - 60*mf
    sf = int(tf*100)/100.
    if hf > 0: CPUtime = '%dh%dm%gs'%(hf,mf,sf)
    else: CPUtime = '%dm%gs'%(mf,sf)
    return CPUtime

#==============================================================================
# Lance un seul test unitaire de module
#==============================================================================
def runSingleUnitaryTest(no, module, test):
    global TESTS
    testr = os.path.splitext(test)
    modulesDir = MODULESDIR[module]
    path = '%s/Apps/%s/%s/test'%(CASSIOPEE, modulesDir, module)

    m1 = expTest1.search(test) # seq ou distribue

    pythonExec = os.getenv('PYTHONEXE', 'python')
    nthreads = KCore.kcore.getOmpMaxThreads()
    nreps = Repeats.get()
    bktest = "bk_{0}".format(test) # backup

    if mySystem == 'mingw' or mySystem == 'windows':
        # Commande Dos (sans time)
        path = path.replace('/', '\\')
        if m1 is not None: cmd = 'cd %s && %s %s'%(path, pythonExec, test)
        else: cmd = 'cd %s && set OMP_NUM_THREADS=%d && mpiexec -np 2 %s %s'%(path, nthreads//2, pythonExec, test)
        cmd2 = 'echo %time%'
    else:
        # Unix - le shell doit avoir l'environnement cassiopee
        #sformat = r'"real\t%E\nuser\t%U\nsys\t%S"'
        cmdReps = ""; test2 = test
        if nreps > 1:
          test2 = "tmp_{0}".format(test)
          cmpiImport = '\nimport Converter.Mpi as Cmpi\n'
          cmpiTrace = '\nCmpi.trace(">>> Iteration: ", cpu=False, stdout=False, filename="proc_{0}.txt")\n\n'.format(test[:-3])
          cmdReps = "rm -f proc_{5}???.txt tmp_{0} tmp2_{0}; cp {0} {4}; cp {0} tmp2_{0}; "\
                    "sed -i '/^test\.test/d' tmp2_{0}; sed -i 's/test\.test.*/pass/g' tmp2_{0};"\
                    "echo '{1}' > tmp_{0}; cat {0} >> tmp_{0}; echo '{2}' >> tmp_{0}; "\
                    "for i in {{2..{3}}}; do echo '{1}'; cat tmp2_{0}; echo '{2}'; done >> tmp_{0};"\
                    "mv tmp_{0} {0};".format(test, cmpiImport, cmpiTrace, nreps, bktest, test[:-3])
                
        if m1 is not None: cmd = 'cd %s; %s time %s %s'%(path, cmdReps, pythonExec, test)
        else: cmd = 'cd %s; %s time kpython -n 2 -t %d %s'%(path, cmdReps, nthreads//2, test)

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
        success = 1
        if regDiff.search(output) is not None: success = 0
        if regFailed.search(output) is not None: success = 0
        if regError.search(output) is not None: success = 0
        if regErreur.search(output) is not None: success = 0
        if regAbort.search(output) is not None: success = 0
        if regSegmentation.search(output) is not None: success = 0
        
        # Recupere l'utilisation memoire lors des runs successifs d'un meme cas
        # unitaire si celui-ci est OK. Pas de check memoire sur les cas FAILED
        if nreps > 1:
            if success == 1:
                tolMem = 0.1
                getMemCmd = "cd {0}; ls; cut -f 2 -d '[' proc_{1}000.txt | cut -f 1 -d ' ' > tmp2_{2};".format(path, test[:-3], test)
                _ = check_output(getMemCmd, shell=True, stderr=subprocess.STDOUT)
                memData = np.loadtxt("{0}/tmp2_{1}".format(path, test))
                if memData.size > 0:
                    # MEMFAILED if the delta mem if greater than a tolerance
                    # for each successive rerun
                    relDMem = np.diff(memData)/memData[:-1]*100.
                    print("Successive relative MEM increments (%)", relDMem)
                    if np.all(relDMem > tolMem): success = -1

        # Recupere le CPU time
        if mySystem == 'mingw' or mySystem == 'windows':
            try: CPUtime = extractCPUTime(output1, output2, nreps=nreps)
            except: CPUtime = 'Unknown'
        else:
            i1 = output.find('\nreal')
            if i1 == -1: CPUtime = 'Unknown'
            else:
                try: CPUtime = extractCPUTime2(output, nreps=nreps)
                except: CPUtime = 'Unknown'
                #CPUtime = output[i1+5:i1+15]; CPUtime = CPUtime.strip()

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
        success = 0; CPUtime = 'Unknown'; coverage='0%' # Core dump/error

    except Exception as e:
        print(e)
        success = 0; CPUtime = 'Unknown'; coverage='0%' # Core dump/error
        
    finally:
        if nreps > 1 and not (mySystem == 'mingw' or mySystem == 'windows'):
            cleanCmd = "cd {0}; mv {1} {2}; rm -f tmp_{2} tmp2_{2} proc_{3}*.txt;".format(
                path, bktest, test, test[:-3])
            _ = check_output(cleanCmd, shell=True, stderr=subprocess.STDOUT)

    # update le fichier .time (si non present)
    fileTime = '%s/Data/%s.time'%(path, testr[0])
    if not os.access(fileTime, os.F_OK):
        writeTime(fileTime, CPUtime, coverage)
        
    # Recupere le tag
    fileStar = '%s/Data/%s.star'%(path, testr[0])
    tag = ' '
    if os.access(fileStar, os.R_OK):
        tag = readStar(fileStar)

    # update status
    if success == 1: status = 'OK'
    elif success == -1: status = 'MEMFAILED'
    else: status = 'FAILED'
    s = buildString(module, test, CPUtime, coverage, status, tag)
    regTest = re.compile(' '+test+' ')
    regModule = re.compile(module+' ')
    for c, tt in enumerate(TESTS):
        if regModule.search(tt) is not None:
            if regTest.search(tt) is not None: TESTS[c] = s; break
    listbox.delete(no, no)
    listbox.insert(no, s)
    listbox.update()
    CPUtime = string2Time(CPUtime)
    return CPUtime

#==============================================================================
# Lance un seul test de la base CFD (CFDBase)
# module = 'CFDBase'
# test = nom du repertoire du cas CFD
#==============================================================================
def runSingleCFDTest(no, module, test):
    global TESTS
    print('Info: Running CFD test %s.'%test)
    path = CASSIOPEE+CFDBASEPATH+'/'+test

    m1 = None # si None=seq
    if test == 'RAE2822_IBC': m1 = True
    try: import mpi4py
    except: m1 = None

    if mySystem == 'mingw' or mySystem == 'windows':
        # Commande Dos (sans time)
        path = path.replace('/', '\\')
        if m1 is None: cmd = 'cd %s && ./valid check'%(path)
        else: cmd = 'cd %s && ./valid check 0 0 0 2 4'%(path)
        cmd2 = 'echo %time%'
    else:
        # Unix - le shell doit avoir l'environnement cassiopee
        if m1 is None: cmd = 'cd %s; ./valid check'%(path)
        else: cmd = 'cd %s; ./valid check 0 0 0 2 4'%(path)
    try:
        if mySystem == 'mingw' or mySystem == 'windows':
            output1 = check_output(cmd2, shell=True, stderr=subprocess.STDOUT)
            output1 = output1.decode()
        output = check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        output = output.decode()
        if mySystem == 'mingw' or mySystem == 'windows':
            output2 = check_output(cmd2, shell=True, stderr=subprocess.STDOUT)
            output2 = output2.decode()
        
        print(output)

        # Recupere success/failed
        success = True
        if regDiff.search(output) is not None: success = False
        if regFailed.search(output) is not None: success = False
        if regError.search(output) is not None: success = False
        if regErreur.search(output) is not None: success = False
        if regAbort.search(output) is not None: success = False
        if regSegmentation.search(output) is not None: success = False

        # Recupere le CPU time
        if mySystem == 'mingw' or mySystem == 'windows':
            CPUtime = extractCPUTime(output1, output2)
        else:
            i1 = output.find('real')
            if i1 == -1: CPUtime = 'Unknown'
            else:
                try: CPUtime = extractCPUTime2(output)
                except: CPUtime = 'Unknown'
                #CPUtime = output[i1+4:i1+14]; CPUtime = CPUtime.strip()
        # Recupere le coverage
        coverage = '100%'
    except Exception as e:
        print(e)
        success = False; CPUtime = 'Unknown'; coverage='0%' # Core dump/error

    # update le fichier .time (si non present)
    fileTime = '%s/Data/%s.time'%(path, test)
    if not os.access(fileTime, os.F_OK):
        writeTime(fileTime, CPUtime, coverage)
        
    # Recupere le tag
    fileStar = '%s/Data/%s.star'%(path, test)
    tag = ' '
    if os.access(fileStar, os.R_OK):
        tag = readStar(fileStar)

    # update status
    if success: status = 'OK'
    else: status = 'FAILED'
    s = buildString(module, test, CPUtime, coverage, status, tag)
    regTest = re.compile(' '+test+' ')
    regModule = re.compile(module+' ')
    for c, tt in enumerate(TESTS):
        if regModule.search(tt) is not None:
            if regTest.search(tt) is not None: TESTS[c] = s; break
    listbox.delete(no, no)
    listbox.insert(no, s)
    listbox.update()
    CPUtime = string2Time(CPUtime)
    return CPUtime

#==============================================================================
# Recupere le nbre de tests selectionnes et le temps total correspondant
#==============================================================================
def getTestsTime():
    selection = listbox.curselection()
    total = len(selection)
    remaining = 0.
    for s in selection:
        t = listbox.get(s)
        splits = t.split(separator)
        remaining += string2Time(splits[3])
    return (total, remaining)

#==============================================================================
# Run selected tests
# Update TESTS, update listbox, update progression
#==============================================================================
def runTests():
    global STOP, THREAD
    selection = listbox.curselection()
    displayStatus(1)
    current = 0
    (total, remaining) = getTestsTime()
    elapsed = 0.

    for s in selection:
        no = int(s)
        t = listbox.get(s)
        splits = t.split(separator)
        module = splits[0]
        test = splits[1]
        module = module.strip()
        test = test.strip()
        current += 1; displayProgress(current, total, remaining, elapsed)
        remaining -= string2Time(splits[3])
        CPUtime = runSingleTest(no, module, test)
        elapsed += CPUtime # real elapsed time
        if STOP == 1: STOP = 0; displayStatus(0); return
    displayStatus(0)
    THREAD=None
    if len(selection) == len(TESTS): notifyValidOK()
    writeSessionLog()
    
def runTestsInThread():
    global THREAD, STOP
    if THREAD is not None: return
    STOP = 0
    THREAD = threading.Thread(target=runTests)
    THREAD.start()

#==============================================================================
# Update the data base for selected tests
#==============================================================================
def updateTests():
    # Supprime les references
    selection = listbox.curselection()
    for s in selection:
        t = listbox.get(s)
        splits = t.split(separator)
        module = splits[0]
        test = splits[1]
        module = module.strip()
        test = test.strip()
        if module == 'CFDBase':
            pathl = CASSIOPEE+CFDBASEPATH+'/'+test
            test2 = test+'.time'
            test = 'post.ref*'
        else:
            modulesDir = MODULESDIR[module]
            path = CASSIOPEE+'/Apps/'+modulesDir
            d = os.path.splitext(test)
            test = d[0]+'.ref*'
            test2 = d[0]+'.time'
            pathl = '%s/%s/test'%(path,module)
        rmFile(pathl, test)
        rmFile(pathl, test2)
    # Set le nombre de fois qu'un cas unitaire doit etre execute a 1
    nreps = Repeats.get()
    Repeats.set(1)
    # Run les tests
    runTests()
    # Reset le nombre de fois qu'un cas unitaire doit etre execute a sa valeur initiale
    Repeats.set(nreps)

def updateTestsInThread():
    global THREAD
    if THREAD is not None: return
    THREAD = threading.Thread(target=updateTests)
    THREAD.start()

#==============================================================================
# Supprime un fichier
# IN: path: file path
# IN: file: file name
#==============================================================================
def rmFile(path, fileName):
    if mySystem == 'mingw' or mySystem == 'windows':
        path = path.replace('/', '\\')
        cmd = 'cd '+path+' && del Data\\'+fileName
    else:
        cmd = 'cd '+path+'; rm -f Data/'+fileName
    try:
        subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)
    except: pass

#==============================================================================
# Construit la liste des tests
# Update TESTS et la listBox
#==============================================================================
def buildTestList(loadSession=False, modules=[]):
    global TESTS
    TESTS = []
    listbox.delete(0, TK.END)
    if not modules:
        modules = getModules()
    # Read last sessionLog conditionally
    ncolumns = 8
    logname = CASSIOPEE+'/Apps/Modules/ValidData/lastSession.log'
    if loadSession and os.access(logname, os.R_OK) and os.path.getsize(logname) > 0:
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
        logname = CASSIOPEE+'/Apps/Modules/ValidData/session.log'
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
            if loadSession and arr.size:
                testArr = arr[np.logical_and(arr[:,0] == m, arr[:,1] == t)]
                if testArr.size:
                    # Args are CPU time, Coverage, Status, and Tag if present
                    if ncolumns == 8:
                        if testArr[0][6].strip() in ['OK', 'FAILED', 'FAILEDMEM', '...']:
                            args = testArr[0][[2,5,6,7]]
                        else: args = testArr[0][[2,5,7,6]]
                    else: args = testArr[0][[2,5,6]]
                    s = buildString(m, t, *args)
                else:
                    s = buildString(m, t)
            else:
                s = buildString(m, t)
            TESTS.append(s)
            listbox.insert(TK.END, s)
    if loadSession and arr.size: writeSessionLog()
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)

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
                    elif tmpFiltr == 'RUN': outFilters.update(['&/!FAILED', '&/!FAILEDMEM', '&/!OK'])
                    elif tmpFiltr == 'UNRUN': outFilters.update(['/FAILED', '/FAILEDMEM', '/OK'])
                    elif tmpFiltr == 'TAG': outFilters.add('£^(?!\*)')
                    elif tmpFiltr == 'UNTAG': outFilters.add('£\*')
                else:
                    if tmpFiltr == 'SEQ': outFilters.add('&t.$')
                    elif tmpFiltr == 'DIST': outFilters.add('&m.$')
                    elif tmpFiltr == 'RUN': outFilters.update(['/FAILED', '/FAILEDMEM', '/OK'])
                    elif tmpFiltr == 'UNRUN': outFilters.update(['&/!FAILED', '&/!FAILEDMEM', '&/!OK'])
                    elif tmpFiltr == 'TAG': outFilters.add('£\*')
                    elif tmpFiltr == 'UNTAG': outFilters.add('£^(?!\*)')
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
        if (not filtr) or (filtr in ['#', '/', '!', '£', '%']) or (filtr[0] in ['&', '*']):
            continue
        shift = 1; endidx = 0
        if filtr[0] == '#': pos = 0 # filter modules
        elif filtr[0] == '/': pos = 7 # filter statuses
        elif filtr[0] == '£': pos = 6 # filter tags
        elif filtr[0] == '%': pos = 5 # filter coverage
        else: pos = 1; shift = 0; endidx = -3 # filter test names
        for s in TESTS:
            strg = s.split(separator)[pos].strip()
            if endidx != 0: strg = strg[:endidx]
            if filtr[shift] == '!':
                if re.search(filtr[1+shift:], strg) is None:
                    filteredTests.add(s)
            elif re.search(filtr[shift:], strg) is not None:
                filteredTests.add(s)
    
    # Apply filters with an AND gate to remove strings from set
    insertedTests = filteredTests.copy()
    for filtr in filters:
        if len(filtr) > 1 and filtr[0] == '&':
            shift = 1; endidx = 0
            if filtr[1] == '#': pos = 0 # filter modules
            elif filtr[1] == '/': pos = 7 # filter statuses
            elif filtr[1] == '£': pos = 6 # filter tags
            elif filtr[1] == '%': pos = 5 # filter coverage
            else: pos = 1; shift = 0; endidx = -3 # filter test names
            for s in filteredTests:
                strg = s.split(separator)[pos].strip()
                if endidx != 0: strg = strg[:endidx]
                if len(filtr) < 3: continue
                if filtr[1+shift] == '!':
                    if len(filtr) > 3 and re.search(filtr[2+shift:], strg) is not None:
                        insertedTests.discard(s)
                elif re.search(filtr[1+shift:], strg) is None:
                    insertedTests.discard(s)
        
    listbox.delete(0, TK.END)
    if filters:
        for s in sorted(insertedTests): listbox.insert(TK.END, s)
    else:
        for s in TESTS: listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    return True

#==============================================================================
# Ouvre un editeur sur le test (emacs)
#==============================================================================
def viewTest(event=None):
    selection = listbox.curselection()
    for s in selection:
        t = listbox.get(s)
        splits = t.split(separator)
        module = splits[0]
        test = splits[1]
        module = module.strip()
        test = test.strip()
        if module == 'CFDBase':
            pathl = CASSIOPEE+CFDBASEPATH+'/'+test
            test = 'compute.py'
        else:
            modulesDir = MODULESDIR[module]
            path = CASSIOPEE+'/Apps/'+modulesDir
            pathl = '%s/%s/test'%(path,module)
        if mySystem == 'mingw' or mySystem == 'windows':
            pathl = pathl.replace('/', '\\')
            cmd = 'cd '+pathl+' && emacs '+test
        else:
            cmd = 'cd '+pathl+'; emacs '+test
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

#==============================================================================
# Selectionne les tests affiche
# Met a jour les infos de progression
#==============================================================================
def selectAll(event=None):
    listbox.selection_set(0, TK.END)
    (total, remaining) = getTestsTime()
    displayProgress(0, total, remaining, 0.)

#==============================================================================
# Affiche les test FAILED ou MEMFAILED dans la listbox
#==============================================================================
def showFilter(filter='FAILED'):
    listbox.delete(0, TK.END)
    for s in TESTS:
        if re.search(filter, s) is not None:
            listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    Filter.set(''); text.update()
    return True

#==============================================================================
# Affiche les test qui ont deja tournes dans la listbox
#==============================================================================
def showRunCases():
    filter = '\.\.\.'
    listbox.delete(0, TK.END)
    for s in TESTS:
        if re.search(filter, s) is None:
            listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    Filter.set(''); text.update()
    return True

#==============================================================================
# Affiche les test qui n'ont deja tournes dans la listbox
#==============================================================================
def showUnrunCases():
    filter = '\.\.\.'
    listbox.delete(0, TK.END)
    for s in TESTS:
        if re.search(filter, s) is not None:
            listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    Filter.set(''); text.update()
    return True

#=============================================================================
# Affiche les tests couverts
#==============================================================================
def showCovered():
    filter = '100%'
    listbox.delete(0, TK.END)
    for s in TESTS:
        if re.search(filter, s) is not None:
            listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    Filter.set(''); text.update()
    return True

#=============================================================================
# Affiche les tests non couverts (0%)
#==============================================================================
def showUncovered():
    filter = ' 0%'
    listbox.delete(0, TK.END)
    for s in TESTS:
        if re.search(filter, s) is not None:
            listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    Filter.set(''); text.update()
    return True

#==============================================================================
# Affiche les tests partiellement couverts
#==============================================================================
def showPartialCovered():
    filter1 = '100%'
    filter2 = ' 0%'
    filter3 = '\.%'
    listbox.delete(0, TK.END)
    for s in TESTS:
        if (re.search(filter1, s) is None and re.search(filter2, s) is None
            and re.search(filter3, s) is None):
            listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    Filter.set(''); text.update()
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
    listbox.delete(0, TK.END)
    for s in TESTS:
        s1 = s.split(separator)
        t1 = s1[2]; t2 = s1[3]
        t1 = string2Time(t1) # new time
        t2 = string2Time(t2)
        if t1 > 0 and t2 > 0:
            if t1 < t2-0.15*t2: listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    Filter.set(''); text.update()
    return True
def showFasterP():
    listbox.delete(0, TK.END)
    for s in TESTS:
        s1 = s.split(separator)
        t1 = s1[2]; t2 = s1[3]
        t1 = string2Time(t1) # new time
        t2 = string2Time(t2)
        if t1 > 0 and t2 > 0:
            if t1 < t2-0.5*t2: listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    Filter.set(''); text.update()
    return True

#==============================================================================
# Affiche les tests plus lent que la reference de 15%
#==============================================================================
def showSlower():
    listbox.delete(0, TK.END)
    for s in TESTS:
        s1 = s.split(separator)
        t1 = s1[2]; t2 = s1[3]
        t1 = string2Time(t1) # new time
        t2 = string2Time(t2)
        if t1 > 0 and t2 > 0:
            if (t1 > t2+0.15*t2): listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    Filter.set(''); text.update()
    return True
def showSlowerP():
    listbox.delete(0, TK.END)
    for s in TESTS:
        s1 = s.split(separator)
        t1 = s1[2]; t2 = s1[3]
        t1 = string2Time(t1) # new time
        t2 = string2Time(t2)
        if t1 > 0 and t2 > 0:
            if t1 > t2+0.5*t2: listbox.insert(TK.END, s)
    listbox.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    Filter.set(''); text.update()
    return True


#==============================================================================
# Affiche tous les tests
#==============================================================================
def showAll():
    Filter.set(''); text.update()
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
    if status == 0: Status.set('Stopped'); label.config(bg='red')
    else: Status.set('Running'); label.config(bg='green')
    label.update()

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
    progressLabel.update()
        
#==============================================================================
# Modifie le nbre de fois qu'un test unitaire doit etre execute.
# De multiples repetitions permettent de detecter d'eventuelles fuites memoire
#==============================================================================
def setNRepeats(event=None):
    if mySystem == 'mingw' or mySystem == 'windows':
      # Feature not implemented for non-Unix OS
      Repeats.set(1)
    else:
      nr = Repeats.get()
      try:
          nri = int(nr)
          print('Info: Number of times each test gets executed set to %d.\n'%nri)
      except:
          Repeats.set(1)
          print('Info: Bad unit test repetition number.\n')
        
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
    text.update()

#==============================================================================
# Exporte les resultats de la valid dans un fichier texte
#==============================================================================
def export2Text():
    try: import tkFileDialog
    except: import tkinter.filedialog as tkFileDialog 
    ret = tkFileDialog.asksaveasfilename()
    if ret == '' or ret is None or ret == (): # user cancel
        return

    file = open(ret, 'w')
    for t in TESTS: file.write(t); file.write('\n')
    file.close()

#=======================================
# writeSessionLog: write log and baseTime
#=======================================
def writeSessionLog():
    svnVersion = 'Unknown'
    if CHECKSVNVERSION:
        try:
            CASSIOPEEL = CASSIOPEE.replace('D:', '/d/') # patch pour msys2/CB
            svnInfo = subprocess.check_output("svn info %s/Apps/Modules"%CASSIOPEEL, shell=True)
            svnInfo = svnInfo.decode('utf-8', 'ignore')
            ss = svnInfo.split('\n')
            for s in ss:
                t = s.split(':')
                if 'vision' in t[0]: svnVersion = t[1]
        except: pass

    messageText = 'Base from'+CASSIOPEE+'\n'
    messageText += 'Based on version %s (can be locally modified).\n'%svnVersion
    for t in TESTS:
        messageText += t+'\n'

    # Write time stamp dans ValidData/base.time et
    # log dans ValidData/session.log
    cassiopee = os.getenv("CASSIOPEE")
    if not os.path.exists(cassiopee+'/Apps/Modules/ValidData'):
        os.mkdir(cassiopee+'/Apps/Modules/ValidData')
    writeFinal(cassiopee+'/Apps/Modules/ValidData/base.time', svnVersion)
    writeFinal(cassiopee+'/Apps/Modules/ValidData/session.log', svnVersion, 
               messageText, append=True)

#=======================================
# Notify "Commit ready" 
#=======================================
def notifyValidOK():
    try:
        import smtplib
        from email.mime.text import MIMEText
        from email.mime.multipart import MIMEMultipart
    except: return
    svnVersion = 'Unknown'
    if CHECKSVNVERSION:
        try:
            CASSIOPEEL = CASSIOPEE.replace('D:', '/d/') # patch pour msys2/CB
            svnInfo = subprocess.check_output("svn info %s/Apps/Modules"%CASSIOPEEL, shell=True)
            svnInfo = svnInfo.decode('utf-8', 'ignore')
            ss = svnInfo.split('\n')
            for s in ss:
                t = s.split(':')
                if 'vision' in t[0]: svnVersion = t[1]
        except: pass

    messageText = 'Base from'+CASSIOPEE+'\n'
    messageText += 'Based on version %s (can be locally modified).\n'%svnVersion
    for t in TESTS:
        messageText += t+'\n'

    try:
        me = os.getenv('USER')+'@onera.fr'
        if me is None: me = ''
        to = 'Christophe.Benoit@onera.fr'
        msg = MIMEMultipart()
        msg['Subject'] = '[Cassiopee] Ready to commit'
        msg['From'] = me
        msg['To'] = to
        msg.preamble = 'Send by Cassiopee.'
        if messageText != '':
            msg.attach(MIMEText(messageText, 'plain'))
        s = smtplib.SMTP('localhost')
        s.sendmail(me, to, msg.as_string())
        s.quit()
    except: pass
    
#==============================================================================
def Quit(event=None):
    import os
    import shutil
    dirname = CASSIOPEE+"/Apps/Modules/ValidData/"
    logname = dirname+"session.log"
    # The session log is copied if it is not empty and if we have write
    # permissions
    if os.access(dirname, os.W_OK) and (not os.path.getsize(logname) == 0):
        dst = CASSIOPEE+"/Apps/Modules/ValidData/lastSession.log"
        shutil.copyfile(logname, dst)
    os._exit(0)

#==============================================================================
# Ajoute une etoile a la selection
#==============================================================================
def tagSelection(event=None):
    global TESTS
    selection = listbox.curselection()
    for s in selection:
        no = int(s)
        t = listbox.get(s)
        splits = t.split(separator)
        module = splits[0].strip()
        test = splits[1].strip()
        modulesDir = MODULESDIR[module]
        path = CASSIOPEE+'/Apps/'+modulesDir+'/'+module+'/test'
        testr = os.path.splitext(test)
        fileStar = path+'/Data/'+testr[0]+'.star'
        writeStar(fileStar, '*')
        splits[7] = ' *'.ljust(4)
        s = separator.join(i for i in splits)
        regTest = re.compile(' '+test+' ')
        regModule = re.compile(module+' ')
        for c, tt in enumerate(TESTS):
            if regModule.search(tt) is not None:
                if regTest.search(tt) is not None: TESTS[c] = s; break
        listbox.delete(no, no)
        listbox.insert(no, s)
        listbox.update()
    return

def untagSelection(event=None):
    global TESTS
    selection = listbox.curselection()
    for s in selection:
        no = int(s)
        t = listbox.get(s)
        splits = t.split(separator)
        module = splits[0].strip()
        test = splits[1].strip()
        modulesDir = MODULESDIR[module]
        path = CASSIOPEE+'/Apps/'+modulesDir+'/'+module+'/test'
        testr = os.path.splitext(test)
        fileStar = path+'/Data/'+testr[0]+'.star'
        rmFile(path, testr[0]+'.star')
        splits[7] = ' '*5
        s = separator.join(i for i in splits)
        regTest = re.compile(' '+test+' ')
        regModule = re.compile(module+' ')
        for c, tt in enumerate(TESTS):
            if regModule.search(tt) is not None:
                if regTest.search(tt) is not None: TESTS[c] = s; break
        listbox.delete(no, no)
        listbox.insert(no, s)
        listbox.update()
    return

#===================================
# Setup for use of global data base
#===================================
def setupGlobal():
    global CASSIOPEE
    CASSIOPEE = os.getenv('CASSIOPEE')
    # Create Local directory for valid products
    if not os.path.exists(CASSIOPEE+'/Apps/Modules/ValidData'):
        os.mkdir(CASSIOPEE+'/Apps/Modules/ValidData')
    os.environ['VALIDLOCAL'] = CASSIOPEE+'/Apps/Modules/ValidData'
    # Change to global ref
    CASSIOPEE = '/stck/benoit/Cassiopee'
    # No update on global ref!
    WIDGETS['updateButton'].configure(state=TK.DISABLED)
    # Change also to match the numthreads of global
    try:
        file = open('/stck/benoit/Cassiopee/Apps/Modules/ValidData/base.time')
        d = file.read(); d = d.split('\n')
        Threads.set(d[2])
        setThreads()
    except: pass

    buildTestList()

def setupLocal():
    global CASSIOPEE
    os.environ['VALIDLOCAL'] = '.'
    CASSIOPEE = os.getenv('CASSIOPEE')
    WIDGETS['updateButton'].configure(state=TK.NORMAL)
    buildTestList()

#===================================================================================
# Filter modules and tests to display. Tests can either be Sequential or Distributed
#===================================================================================
def filterModulesTests(master, event=None):
    def _onDeselectAllModules(modSwitches):
        [sw.set(0) for sw in modSwitches]
    def _onDeselectPModules(modules, modSwitches):
        [sw.set(0) for m, sw in zip(modules, modSwitches)
            if (m.startswith("Fast") or m.startswith("Apps") or m.startswith("FF"))]
    def _onSelectAllModules(modSwitches):
        [sw.set(1) for sw in modSwitches]
    def _onUpdating(modules, modSwitches, testSwitches):
        global TESTS_FILTER
        values = [testSwitches[i].get() for i in range(2)]
        if sum(values) == 2:
            TESTS_FILTER = 0
        elif values[1] == 1:
            TESTS_FILTER = 2
        else:
            TESTS_FILTER = 1
        modules = [m for i, m in enumerate(modules) if modSwitches[i].get() == 1]
        buildTestList(modules=modules)
        newWin.destroy()
    
    newWin = TK.Toplevel(master)
    newWin.title('Filter Modules and Tests')
    newWin.geometry("460x775")

    # Define new frame with lists of modules and tests to load as sub-frames
    newFrame = TK.LabelFrame(newWin)
    newFrame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
    
    leftFrame = TK.LabelFrame(newFrame, text="Modules to load:")
    leftFrame.grid(row=0, column=0, padx=5, pady=10, sticky="nsew")
    rightFrame = TK.LabelFrame(newFrame, text="Tests to load:")
    rightFrame.grid(row=0, column=1, padx=5, pady=10, sticky="nsew")
    
    # > Modules
    modules = getModules()
    nmodules = len(modules)
    modSwitches = [TK.IntVar(value=1) for _ in range(nmodules)]
    
    for i in range(nmodules):
        label = TK.Label(leftFrame, text=modules[i], anchor="w")
        label.grid(row=i, column=1, sticky="w")
        button = TK.Checkbutton(leftFrame, variable=modSwitches[i],
                                onvalue=1, offvalue=0, height=1, width=5)
        button.grid(row=i, column=0, sticky="w")
        
    for i in range(nmodules):
        leftFrame.grid_rowconfigure(i, weight=1)
    leftFrame.grid_columnconfigure(0, weight=1)
    leftFrame.grid_columnconfigure(1, weight=1)
        
    # > Tests
    if not isMpi or TESTS_FILTER == 1: values = [1, 0]
    elif TESTS_FILTER == 2: values = [0, 1]
    else: values = [1, 1]
    testSwitches = [TK.IntVar(value=values[i]) for i in range(2)]
    labels = ["Sequential", "Distributed"]
    
    for i in range(2):
        label = TK.Label(rightFrame, text=labels[i], anchor="w")
        label.grid(row=i, column=1, sticky="nsew")
        button = TK.Checkbutton(rightFrame, variable=testSwitches[i],
                                onvalue=1, offvalue=0, height=1, width=5)
        if i == 1 and not isMpi: button.configure(state='disabled')
        button.grid(row=i, column=0, sticky="nsew")
    
    rightFrame.grid_columnconfigure(0, weight=1)
    rightFrame.grid_columnconfigure(1, weight=1)
        
    # > Buttons at the bottom of the new frame
    deselectAllModulesWithArgs = partial(_onDeselectAllModules, modSwitches)
    deselectPModulesWithArgs = partial(_onDeselectPModules, modules, modSwitches)
    selectAllModulesWithArgs = partial(_onSelectAllModules, modSwitches)
    
    updateWithArgs = partial(_onUpdating, modules, modSwitches, testSwitches)
    
    button = TK.Button(newFrame, text='Deselect All',
                       command=deselectAllModulesWithArgs,
                       height=1, fg='black', bg='white')
    button.grid(row=1, column=0, padx=5, pady=1, sticky="nsew")
    button = TK.Button(newFrame, text='Deselect PModules',
                       command=deselectPModulesWithArgs,
                       height=1, fg='black', bg='white')
    button.grid(row=1, column=1, padx=5, pady=1, sticky="nsew")
    button = TK.Button(newFrame, text='Select All',
                       command=selectAllModulesWithArgs,
                       height=1, fg='black', bg='white')
    button.grid(row=2, column=0, padx=5, pady=1, sticky="nsew")
    button = TK.Button(newFrame, text='Update', command=updateWithArgs,
                       height=1, fg='blue', bg='white')
    button.grid(row=2, column=1, padx=5, pady=1, sticky="nsew")
    
    for i in range(2):
        newFrame.grid_columnconfigure(i, weight=1, minsize=220)
                             
    newWin.protocol("WM_DELETE_WINDOW", newWin.destroy)
    newWin.mainloop()
    
#==============================================================================
# Main
#==============================================================================

# Main window
master = TK.Tk()
master.title('*Cassiopee* valid @ '+machine)
master.columnconfigure(0, weight=1)
master.rowconfigure(0, weight=1)
#GENERALFONT = ('Courier', 9)
GENERALFONT = ('Andale Mono', 9)

master.option_add('*Font', GENERALFONT)
generalFont = Font.Font(family=GENERALFONT[0], size=GENERALFONT[1])
generalFontS = generalFont.measure(' ')*1.
generalFontA = generalFont.measure('a')*1.
generalFontFixed = generalFont.metrics('fixed')

# Main menu
menu = TK.Menu(master)
file = TK.Menu(menu, tearoff=0)
menu.add_cascade(label='File', menu=file)
tools = TK.Menu(menu, tearoff=0)
menu.add_cascade(label='Tools', menu=tools)
view = TK.Menu(menu, tearoff=0)
menu.add_cascade(label='View', menu=view)

loadSessionWithArgs = partial(buildTestList, True)
file.add_command(label='Load last session', command=loadSessionWithArgs)
file.add_command(label='Purge session', command=buildTestList)
file.add_command(label='Export to text file', command=export2Text)
file.add_command(label='Notify Ready for commit', command=notifyValidOK)
file.add_command(label='Quit', command=Quit, accelerator='Ctrl+Q')
view.add_command(label='Show FAILED', command=showFilter)
showFilterWithArgs = partial(showFilter, "MEMFAILED")
view.add_command(label='Show MEMFAILED', command=showFilterWithArgs)
view.add_command(label='Show ALL tests', command=showAll)
view.add_separator()
view.add_command(label='Show run cases', command=showRunCases)
view.add_command(label='Show unrun cases', command=showUnrunCases)
view.add_separator()
view.add_command(label='Show covered (100%)', command=showCovered)
view.add_command(label='Show partially covered (x%)',
                 command=showPartialCovered)
view.add_command(label='Show uncovered (0%)', command=showUncovered)
view.add_separator()
view.add_command(label='Show faster (-15%)', command=showFaster)
view.add_command(label='Show slower (+15%)', command=showSlower)
view.add_command(label='Show faster (-50%)', command=showFasterP)
view.add_command(label='Show slower (+50%)', command=showSlowerP)
view.add_separator()
view.add_command(label='Select all visible tests', command=selectAll,
                 accelerator='Ctrl+A')

filterModulesTestsWithArgs = partial(filterModulesTests, master)
#filterTestsWithArgs = partial(filterSeqParTests, master)
#tools.add_command(label='Filter tests', command=filterTestsWithArgs)
tools.add_command(label='Filter modules and tests', command=filterModulesTestsWithArgs)
tools.add_separator()
tools.add_command(label='Tag selection', command=tagSelection)
tools.add_command(label='Untag selection', command=untagSelection)
tools.add_separator()

try:
    file = open('/stck/benoit/Cassiopee/Apps/Modules/ValidData/base.time')
    d = file.read(); d = d.split('\n')
    d = ' ['+d[0]+'/'+d[1]+'/'+d[2]+' threads]'
except: d = ''

tools.add_command(label='Switch to global data base'+d, command=setupGlobal)
tools.add_command(label='Switch to local data base', command=setupLocal)

master.config(menu=menu)
master.bind_all("<Control-q>", Quit)
master.protocol("WM_DELETE_WINDOW", Quit)
master.bind_all("<Control-a>", selectAll)

# Main frame
frame = TK.Frame(master)
frame.columnconfigure(0, weight=1)
frame.rowconfigure(0, weight=1)
frame.columnconfigure(1, weight=1)
frame.grid(row=0, column=0, sticky=TK.EW)

listbox = TK.Listbox(frame, selectmode=TK.EXTENDED, width=120, height=39,
                     background='White')
listbox.grid(row=0, column=0, columnspan=11, sticky=TK.NSEW)

scrollbar = TK.Scrollbar(frame, orient=TK.VERTICAL)
scrollbar.grid(row=0, column=11, sticky=TK.NSEW)

Status = TK.StringVar(master)
label = TK.Label(frame, textvariable=Status)
Status.set('Stopped'); label.config(bg='red')
label.grid(row=1, column=0, sticky=TK.EW)

Progress = TK.StringVar(master)
progressLabel = TK.Label(frame, textvariable=Progress)
Progress.set('  0/  0 [0h 0m 0s/0h 0m 0s]')
progressLabel.grid(row=1, column=1, sticky=TK.EW)

Filter = TK.StringVar(master)
text = TK.Entry(frame, textvariable=Filter, background='White', width=50)
text.bind('<KeyRelease>', filterTestList)
text.grid(row=1, column=2, columnspan=3, sticky=TK.EW)
"""Filter tests by this regexp.\n'+'-'*70+'\n'\
  '1) Filters are separated by a white space, for ex.\n'\
  '    ^cylinder ^sphere\nselects tests whose names start with either cylinder or sphere.\n'\
  '2) Filters can be applied on modules by prefixing their names\nwith the symbol #, for ex.\n'\
  '    #Fast #FF #Apps\nwill load all PModules.\n'\
  '3) Filters can be applied on statuses by prefixing their names\nwith the symbol /, for ex.\n'\
  '    /FAILED /FAILEDMEM\nwill select all tests that failed because of bugs or memory leaks.\n'\
  '4) Logical OR operations are performed unless a filter is \n'\
  'prefixed with the symbol & to force an AND operation.\n'\
  '    #Converter &/FAILED\nwill select all buggy tests in the Converter module.\n'\
  '5) Negation of an expression with symbol ! should be the innermost\noperator, for ex.\n'\
  '    #Fast &#!FastC\nwill load Fast modules but FastC.\n'\
  '6) A few keyworded filters exist and can be called using angle\nbrackets, for ex.\n'\
  '    <SEQ> &<UNRUN>\nselects all sequential tests that have not been run yet.\n'\
  'These keyworded filters are:\n    <SEQ>, <DIST>, <RUN>, <UNRUN>, <TAG>, <UNTAG>.\n'\
  '7) Tests can be filtered by coverage using the % symbol, for ex.\n'\
  '    <TAG> &/OK &%100\nselects all bug-free tagged tests that have 100% coverage."""
  
filterInfoBulle = 'Filter test database using a regexp.\n'+'-'*70+'\n'\
  '1) White-spaced: ^cylinder ^sphere\n'\
  '2) Module filter using #: #Fast #FF #Apps\n'\
  '3) Status filter using /: /FAILED /FAILEDMEM\n'\
  '4) Coverage filter using %: %100\n'\
  '5) Keyworded filters: <SEQ>, <DIST>, <RUN>, <UNRUN>, <TAG>, <UNTAG>.\n'\
  '6) Logical OR ops unless prefixed with & (AND): #Converter &/FAILED\n'\
  '7) Negation of a filter using !: #Fast &#!FastC (innermost symbol)'
BB = CTK.infoBulle(parent=text, text=filterInfoBulle)

button = TK.Button(frame, text='Run', command=runTestsInThread, fg='blue')
BB = CTK.infoBulle(parent=button, text='Run selected tests.')
button.grid(row=1, column=5, sticky=TK.EW)
button = TK.Button(frame, text='Stop', command=stopTests, fg='red')
button.grid(row=1, column=6, sticky=TK.EW)
button = TK.Button(frame, text='Update', command=updateTestsInThread, fg='blue')
BB = CTK.infoBulle(parent=button, text='Update tests (replace data base files).')
WIDGETS['updateButton'] = button
button.grid(row=1, column=7, sticky=TK.EW)
button = TK.Button(frame, text='Edit', command=viewTest)
button.grid(row=1, column=8, sticky=TK.EW)

Threads = TK.StringVar(master)
text = TK.Entry(frame, textvariable=Threads, background='White', width=3)
text.grid(row=1, column=9, sticky=TK.EW)
text.bind('<Return>', setThreads)
text.bind('<KP_Enter>', setThreads)
getThreads()
BB = CTK.infoBulle(parent=text, text='Number of threads.')

Repeats = TK.IntVar(master, value=1)
repeatsEntry = TK.Entry(frame, textvariable=Repeats, background='White',
                        width=3)
repeatsEntry.grid(row=1, column=10, sticky=TK.EW)
repeatsEntry.bind('<Return>', setNRepeats)
repeatsEntry.bind('<KP_Enter>', setNRepeats)
BB = CTK.infoBulle(parent=repeatsEntry,
                   text='Number of times each unit test gets executed.')
frame.grid(sticky=TK.NSEW)

buildTestList()

TK.mainloop()
