#
# Support pour le log
#
import Tk as CTK
import Converter.Internal as Internal
import inspect
import re
# Balise de fin de log
endExp = re.compile("Log.end()")
# Log file handler
LOGFILE = None
# Saved log file handler (temporary stop recording)
SAVELOG = None

# LOGGING
LOGGING = True
# LOG string
LOG = []

#==============================================================================
# Open a log file
#==============================================================================
def openLogFile(file):
    global LOGFILE
    if LOGFILE is not None: # deja ouvert
        LOGFILE.close()
        
    LOGFILE = open(file, 'w')
    # Entete
    l = ['import Converter.PyTree as C',
         'import Converter.Internal as Internal',
         'import Generator.PyTree as G',
         'import Transform.PyTree as T']
    w(l)

#==============================================================================
# Close a log file
#==============================================================================
def closeLogFile():
    global LOGFILE
    if LOGFILE is not None: LOGFILE.close()
    LOGFILE = None

#==============================================================================
# Pause log
#==============================================================================
def pause():
    global SAVELOG, LOGFILE
    SAVELOG = LOGFILE; LOGFILE = None

#==============================================================================
# balise start
#==============================================================================
def start():
    if LOGFILE is None: return
    curr = inspect.currentframe()
    f = curr.f_back
    lineStart = f.f_lineno
    file = f.f_code.co_filename
    ptf = open(file, 'r')
    gotoLine(ptf, lineStart)
    l = []
    readToEnd(ptf, l)
    ptf.close()
    w(l)
    return

#==============================================================================
# balise end
#==============================================================================
def end(): return

#==============================================================================
# Write a list of lines to Log file
#==============================================================================
def w(lines):
    if LOGFILE is None: return
    for l in lines:
        LOGFILE.write(l)
        LOGFILE.write('\n')
    LOGFILE.flush()

#==============================================================================
# Passe no lignes dans le fichier en entree
#==============================================================================
def gotoLine(ptf, no):
    c = 0
    while (c < no): res = ptf.readline() ; c += 1

#==============================================================================
# Lit jusqu'a Log.end()
#==============================================================================
def readToEnd(ptf, ret):
    cont = True
    while cont:
        res = ptf.readline()
        rend = endExp.search(res)
        if rend is not None: cont = False
        else: ret.append(res) 

#===============================================================================
# Retourne une liste de chaines contenant le code python de la liste des zones 
# selectionnees (zones)
# Retourne None si pas de zone a traiter
#===============================================================================
def getSelectedZones():
  import CPlot
  nzs = CPlot.getSelectedZones()
  if nzs == []: return None
  if len(nzs) == len(Internal.getZones(CTK.t)): return ['zones = Internal.getZones(t)']
  nobu = -1
  for nz in nzs:
        nob = CTK.Nb[nz]+1
        if nobu == -1: nobu = nob
        elif nobu != nob: nobu = -2; break
  
  if nobu > 0:
    base = Internal.getBases(CTK.t)[nobu-1]
    if len(nzs) == len(Internal.getZones(base)):
        return ['base = Internal.getBases(t)[%d]'%(nobu-1), 
        'zones = Internal.getZones(base)']
  ret = ['zones = []']
  for nz in nzs:
        nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
        ret.append('zones.append(t[2][%d][2][%d])'%(nob,noz))
  return ret

# Idem mais retourne le chemin des zones selectionnees (zonesPath)
def getSelectedZones2():
  import CPlot
  nzs = CPlot.getSelectedZones()
  if nzs == []: return None
  if len(nzs) == len(Internal.getZones(CTK.t)):
    bases = Internal.getBases(CTK.t)
    ret = ['zonesPath = []']
    for b in bases:
        path = Internal.getPath(b)
        ret.append('zonesPath.append(%s)'%path)
    return ret
  nobu = -1
  for nz in nzs:
        nob = CTK.Nb[nz]+1
        if nobu == -1: nobu = nob
        elif nobu != nob: nobu = -2; break
  if nobu > 0:
    path = Internal.getPath(t[2][nobu])     
    return ['zonesPath = [%s]'%path]

  ret = ['zonesPath = []']
  for nz in nzs:
        nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
        path = Internal.getPath(CTK.t, CTK.t[2][nob][2][noz])
        ret.append('zonesPath.append(%s)'%(path))
  return ret

#================================================================================
# Init LOG
#================================================================================
def initLog():
    global LOG
    LOG = ['import Converter.PyTree as C', 
           'import Converter.Internal as Internal']
    if CTK.FILE != '':
        LOG.append('C.convertFile2PyTree(%s)\n'%CTK.FILE)

#=================================================================================
# Ajoute le import module au LOG
#=================================================================================
def addImportModule(moduleName):
    global LOG
    c = 0; modules = []
    for i in LOG:
        sp = i.split(' ')
        if sp[0] == 'import': modules.append(sp[1])
        else: break 
        c += 1
    if moduleName not in modules:
        LOG.insert(c, 'import %s'%moduleName)

#=================================================================================
def displayLog():
    print(LOG)
