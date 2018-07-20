# - Filters -
import Internal
import PyTree
import Converter
from Distributed import convert2PartialTree, _convert2PartialTree, convertFile2SkeletonTree, _readPyTreeFromPaths, _readZones, \
_convert2SkeletonTree, readNodesFromPaths, writeNodesFromPaths, writePyTreeFromPaths, deletePaths

# Prend un fileName, si c'est toto/*, rend la liste des fichiers
def expand(fileName):
  s = fileName.rsplit('/', 1)
  if len(s) > 1 and (s[1] == '*' or s[1] == '*.cgns' or s[1] == '*.hdf'): # multiple
    import glob
    if s[1] == '*': files = glob.glob(s[0]+'/*.cgns')
    else: files = glob.glob(fileName)
  return files

#============================================================================
# Lecture des noms Base/Zones + dims
# Retourne un squelette (depth=3) + la liste des zones path names (znp)
# si baseNames: retourne les znp pour la base donnee
# si familyZoneNames: retourne les znp des zones de famille donnee
# si BCType: retourne les znp des zones contenant une BC de type BCType
# ou de famille 'familySpecified:WALL'
#============================================================================
def readZoneHeaders(fileName, format=None, baseNames=None, familyZoneNames=None, BCType=None):
    a = convertFile2SkeletonTree(fileName, format, maxDepth=3, maxFloatSize=6)
    # filter by base names
    if baseNames is not None:
        if not isinstance(baseNames, list): baseNames = [baseNames]
        znp = []
        for b in baseNames:
            l = Internal.getNodeFromName1(a, b)
            if l is not None:
               zs = Internal.getZones(l)
               for z in zs: znp.append('/'+b+'/'+z[0])
    # filter by zone family names
    elif familyZoneNames is not None:
        if not isinstance(familyZoneNames, list): familyZoneNames = [familyZoneNames]
        znp = []
        bases = Internal.getBases(a)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                f = Internal.getNodeFromType1(z, 'FamilyName_t')
                if f is not None:
                    for fz in familyZoneNames:
                        if Internal.getValue(f) == fz:
                            znp.append('/'+b[0]+'/'+z[0])
    elif BCType is not None:
        families = PyTree.getFamilyBCNamesOfType(a, BCType)
        s = BCType.split(':')
        if len(s) == 2: families.append(s[1])
        print families
        if BCType == 'BCWall':
          families1 = PyTree.getFamilyBCNamesOfType(a, 'BCWallInviscid')
          families2 = PyTree.getFamilyBCNamesOfType(a, 'BCWallViscous*')
        znp = []
        bases = Internal.getBases(a)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                path = '/'+b[0]+'/'+z[0]
                _readPyTreeFromPaths(a, fileName, path+'/ZoneBC')
                nodes = Internal.getNodesFromValue(z, BCType)
                nodes += PyTree.getFamilyBCs(z, families)
                #_convert2SkeletonTree(z)
                if nodes != []: znp.append(path)
    else:
        znp = Internal.getZonePaths(a, pyCGNSLike=True)
    return a, znp

#========================================================================
# Load par containeurs
# Load les containeurs "cont" dans a pour les znp donnes
# cont='GridCoordinates', 'FlowSolution'
#========================================================================
def _loadContainers(a, fileName, znp, cont, format=None):
    if isinstance(cont, list): conts = cont
    else: conts = [cont]
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    for p in znps:
        paths = []
        for c in conts: paths.append('%s/%s'%(p,c))
        _readPyTreeFromPaths(a, fileName, paths, format)
    return None

#=========================================================================
# Load connectivity
# Load les connectivites dans a pour les znps donnes
#=========================================================================
def _loadConnectivity(a, fileName, znp, format=None):
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    for p in znps:
        z = Internal.getNodeFromPath(a, p)
        elts = Internal.getNodesFromType1(z, 'Elements_t')
        paths = []
        for e in elts: paths.append(p+'/'+e[0])
        _readPyTreeFromPaths(a, fileName, paths, format)    
    return None

#==========================================================================
# Load par variables
# Load les variables "var" pour les znp donnes
# var='Density', 'centers:Density', 'FlowSolution/Density' or list of vars
#==========================================================================
def _loadVariables(a, fileName, znp, var, format):
    if isinstance(var, list): vars = var
    else: vars = [var]
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    fvars = []; cont = []
    for v in vars:
       s = v.split(':')
       if v[0:10] == 'Coordinate': 
        fvars.append(Internal.__GridCoordinates__+'/'+v); cont = Internal.__GridCoordinates__
       elif len(s) == 2 and s[0] == 'centers': 
        fvars.append(Internal.__FlowSolutionCenters__+'/'+s[1]); cont = Internal.__FlowSolutionCenters__
       elif len(s) == 2 and s[0] == 'nodes': 
        fvars.append(Internal.__FlowSolutionNodes__+'/'+s[1]); cont = Internal.__FlowSolutionNodes__
       else:
        s = v.split('/')
        if len(s) == 2: 
          fvars.append(s[0]+'/'+s[1]); cont = Internal.__FlowSolutionNodes__
        else: 
          fvars.append(Internal.__FlowSolutionNodes__+'/'+v); cont = Internal.__FlowSolutionNodes__

    for p in znps:
        paths = []
        for v in fvars: paths.append('%s/%s'%(p,v))
        _readPyTreeFromPaths(a, fileName, paths, format)
        # Ensure location in containers
        zp = Internal.getNodeFromPath(a, p)
        fp = Internal.getNodeFromPath(a, paths[0])
        if zp is not None and fp is not None:
            c = Internal.getNodeFromName1(zp, cont)
            if c is not None:
              if PyTree.getNPts(zp) != fp[1].size:
                Internal.createUniqueChild(c, 'GridLocation', 'GridLocation_t', 'CellCenter')
    return None

#==================================================================================
# Fully load all path in a
#==================================================================================
def _loadZones(a, fileName, znp, format=None):
  if isinstance(znp, list): znps = znp
  else: znps = [znp]
  _readPyTreeFromPaths(a, fileName, znps, format)

# Fully load zoneBC_t and GridConnectivity_t of znp
def _loadZoneBCs(a, fileName, znp, format=None):
  if isinstance(znp, list): znps = znp
  else: znps = [znp]
  paths = []
  for p in znps:
    n = Internal.getNodeFromPath(a, p)
    if n is not None:
      children = n[2]
      for i in children:
        if i[3] == 'ZoneBC_t': paths.append(p+'/'+i[0])
        if i[3] == 'GridConnectivity_t': paths.append(p+'/'+i[0])
  if paths != []: _readPyTreeFromPaths(a, fileName, paths, format)
  return None

# Load zone BC but without loading BCDataSet fields
def _loadZoneBCsWoData(a, fileName, znp, format=None):
  if isinstance(znp, list): znps = znp
  else: znps = [znp]
  paths = []
  for p in znps:
    n = Internal.getNodeFromPath(a, p)
    if n is not None:
      children = n[2]
      for i in children:
        if i[3] == 'GridConnectivity_t': paths.append(p+'/'+i[0])
  if paths != []: _readPyTreeFromPaths(a, fileName, paths, format)
  paths = []
  for p in znps:
    n = Internal.getNodeFromPath(a, p)
    if n is not None:
      children = n[2]
      for i in children:
        if i[3] == 'ZoneBC_t': paths.append(p+'/'+i[0])
  print paths
  if paths != []: _readPyTreeFromPaths(a, fileName, paths, format, maxDepth=2, setOnlyValue=False)
  return None

# Load fully extras node (not Zone or Base)
def _loadExtras(a, fileName, format=None):
  # Fully load nodes of level 1 except CGNSBase_t
  children = a[2]
  paths = []
  for i in children:
    if i[3] != 'CGNSBase_t': paths.append(i[0])
  
  # Fully load nodes of level 2 except Zone_t
  bases = Internal.getBases(a)
  for b in bases:
    children = b[2]
    for i in children:
      if i[3] != 'Zone_t': paths.append(b[0]+'/'+i[0])
  if paths != []: _readPyTreeFromPaths(a, fileName, paths, format)
  return None

# get variables: return a list
# this doesnt mean that all vars are defined in all zones!
def getVariables(fileName, znp, cont=None, format=None):
  if isinstance(znp, list): znps = znp
  else: znps = [znp]
  if cont is None: # check containers in files
    conts = set()
    nodes = readNodesFromPaths(fileName, znp, format, maxDepth=1)
    for n in nodes:
      for j in n[2]:
        if j[3] == 'FlowSolution_t': conts.add(j[0])
    cont = list(conts)
  elif cont == 'centers': cont = [Internal.__FlowSolutionCenters__]
  elif cont == 'nodes': cont = [Internal.__FlowSolutionNodes__]
  elif isinstance(cont, str): cont = [cont]
  paths = []
  for c in cont:
    for p in znps:
      paths.append(p+'/'+c)
  vars = set()
  nodes = readNodesFromPaths(fileName, paths, format, maxDepth=1)
  for n in nodes:
    for j in n[2]:
      if j[3] == 'DataArray_t': vars.add(n[0]+'/'+j[0])
  return list(vars)

# Return the list of variables in BCDataSet
def getBCVariables(a, fileName, znp, cont=None, format=None):
  if isinstance(znp, list): znps = znp
  else: znps = [znp]
  paths = []
  for p in znps:
      zp = Internal.getNodeFromPath(a, p)
      if zp is not None:
        zoneBC = Internal.getNodesFromType1(zp, 'ZoneBC_t')
        for zbc in zoneBC:
          paths.append(p+'/'+zbc[0])
  nodes = readNodesFromPaths(fileName, paths, format, maxDepth=-1, maxFloatSize=0)
  vars = set()
  for n in nodes:
    p = Internal.getNodesFromType(n, 'BCData_t')
    for j in p:
      for k in j[2]:
        if k[3] == 'DataArray_t': vars.add(k[0])
  return list(vars)

# Load un bboxTree
def bboxTree(fileName):
    # Load only bboxes
    return None

# Load only zones that match a bbox
def filterBBox(a, fileName, znp, bbox, F):
    xmin = bb[0]; ymin = bb[1]; zmin = bb[2]
    xmax = bb[3]; ymax = bb[4]; zmax = bb[5]
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    out = []
    for p in znps:
       z = Internal.getNodeFromPath(a, p)
       path = ['%s/CoordinateX']
       _readPyTreeFromPaths(a, fileName, path)
       vmin = C.getMinValue(z, 'CoordinateX')
       vmax = C.getMaxValue(z, 'CoordinateX')
       if vmax < xmin or vmin > xmax: continue
       path = ['%s/CoordinateY']
       _readPyTreeFromPaths(a, fileName, path)
       vmin = C.getMinValue(z, 'CoordinateY')
       vmax = C.getMaxValue(z, 'CoordinateY') 
       if vmax < ymin or vmin > ymax: continue
       path = ['%s/CoordinateZ']
       _readPyTreeFromPaths(a, fileName, path)
       vmin = C.getMinValue(z, 'CoordinateZ')
       vmax = C.getMaxValue(z, 'CoordinateZ') 
       if vmax < zmin or zmin > ymax: continue
       # Perform the action, finish loading
       # ...
       # slice, selectCells, isoSurfMC
       ztype = Internal.getZoneType(z)
       if ztype == 2: _loadConnectivity(a, fileName, p)
       out = F(z)
       # Release eventuellement
       _convert2SkeletonTree(z)
    return out


#==========================================================
class Handle:
  """Handle for filter."""
  def __init__(self, fileName, where=None):
    self._fileName = fileName
    self._fileVars = None # vars existing in file
    self._fileBCVars = None # BCDataSet vars existing in file
    self._znp = None # zone paths existing in file
    self._format = Converter.checkFileType(fileName) # Real format of file

  # Retourne les zones de a
  def getZonePaths(self, a):
    out = []
    bases = Internal.getBases(a)
    for b in bases:
      zones = Internal.getZones(b)
      for z in zones: out.append('/'+b[0]+'/'+z[0])
    return out

  def setZnp(a):
    self._znp = []
    bases = Internal.getBase()
    for b in bases:
      zones = Internal.getZones(b)
      for z in zones: self._znp.append('/'+b[0]+'/'+z[0])

  # Retourne les variables du fichier
  def getVariables(self, a=None, cont=None):
    if a is not None: p = self.getZonePaths(a)
    else: p = self._znp
    vars = getVariables(self._fileName, p, cont, self._format)
    self._fileVars = vars
    return vars

  # Retourne les variables de BCDatSet du fichier
  def getBCVariables(self, a):
    if a is not None: p = self.getZonePaths(a)
    else: p = self._znp
    vars = getBCVariables(a, self._fileName, p, None, self._format)
    self._fileBCVars = vars
    return vars

  # Charge un squelette, stocke les chemins des zones du fichier (znp)
  def loadSkeleton(self):
    a, self._znp = readZoneHeaders(self._fileName, self._format)
    return a

  # Charge les grid coordinates + grid connectivity + BC
  # pour toutes les zones de a
  def _loadGCBC(self, a):
    self._loadContainer(a, 'GridCoordinates')
    self._loadConnectivity(a)
    self._loadZoneBCs(a)
    Internal._fixNGon(a)
    return None

  # Charge tous les noeuds extra de a
  def _loadExtras(self, a):
    _loadExtras(a, self._fileName, self._format)
    return None

  # Charge completement toutes les zones de a ou de chemin fourni
  def _loadZones(self, a, znp=None):
    if znp is None: znp = self.getZonePaths(a)
    _loadZones(a, self._fileName, znp, self._format)
    return None

  # Charge les grid coordinates + grid connectivity + BC
  # pour les zones specifiees
  def _loadZonesWoVars(self, a, znp=None):
    if znp is None: znp = self._znp
    _loadContainers(a, self._fileName, znp, 'GridCoordinates', self._format)
    _loadConnectivity(a, self._fileName, znp, self._format)
    _loadZoneBCs(a, self._fileName, znp, self._format)
    Internal._fixNGon(a)
    return None

  # Charge toutes les BCs (avec BCDataSet) des zones de a  
  def _loadZoneBCs(self, a):
    p = self.getZonePaths(a)
    _loadZoneBCs(a, self._fileName, self._znp, self._format)
    return None

  # Charge toutes les BCs (sans BCDataSet) des zones de a
  def _loadZoneBCsWoData(self, a):
    p = self.getZonePaths(a)
    _loadZoneBCsWoData(a, self._fileName, p, self._format)
    return None

  # Charge la connectivite pour toutes les zones de a
  def _loadConnectivity(self, a):
    p = self.getZonePaths(a)
    _loadConnectivity(a, self._fileName, p, self._format)
    return None

  # Charge le container "cont" pour toutes les zones de a
  def _loadContainer(self, a, cont):
    p = self.getZonePaths(a)
    _loadContainers(a, self._fileName, p, cont, self._format)
    return None

  # Charge la ou les variables "var" pour toutes les zones de a
  def _loadVariables(self, a, var):
    p = self.getZonePaths(a)
    _loadVariables(a, self._fileName, p, var, self._format)
    return None
