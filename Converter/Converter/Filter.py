# - Filters -
from . import Internal
from . import PyTree
from . import Converter
from . import Mpi as Cmpi
from .Distributed import convert2PartialTree, _convert2PartialTree, convert2SkeletonTree, _convert2SkeletonTree, convertFile2SkeletonTree, _readPyTreeFromPaths, readPyTreeFromPaths, _readZones, \
_convert2SkeletonTree, readNodesFromPaths, writeNodesFromPaths, writePyTreeFromPaths, deletePaths, fixPaths__
import numpy

# Prend un fileName, si c'est toto/*, rend la liste des fichiers
def expand(fileName):
  s = fileName.rsplit('/', 1)
  if len(s) > 1 and (s[1] == '*' or s[1] == '*.cgns' or s[1] == '*.hdf'): # multiple
    import glob
    if s[1] == '*': files = glob.glob(s[0]+'/*.cgns')
    else: files = glob.glob(fileName)
  return files

#==============================================================================
# Lit seulement une partie des tableaux d'un fichier a partir de la definition d'un filtre.
# filter est un dictionnaire pour chaque path
# pour les grilles structurees : [[imin,jmin,kmin], [1,1,1], [imax,jmax,kmax], [1,1,1]]
# pour les grilles non structurees : [[istart], [1], [iend], [1]]
# Uniquement HDF
# Retourne un dictionnaire du numpy loade pour chaque path
#==============================================================================
def readNodesFromFilter(fileName, filter, format='bin_hdf', com=None):
  """Read nodes from file given a filter."""
  for i in filter:
    b = fixPaths__([i])[0]
    val = filter.pop(i)
    filter[b] = val
  ret = Converter.converter.convertFile2PartialPyTree(fileName, format, None, com, filter)
  return ret

# Ecrit des tableaux ou des morceaux de tableau a certains endroits du fichier
# definis par filter
# t: pyTree avec les memes chemins
def writePyTreeFromFilter(t, fileName, filter, format='bin_hdf', com=None, skelData=None):
  """Write nodes to file given a filter."""
  for i in filter:
    b = fixPaths__([i])[0]
    val = filter.pop(i)
    filter[b] = val
  Converter.converter.convertPyTree2FilePartial(t, fileName, format, skelData, com, filter)
  return None

#============================================================================
# Lecture des noms Base/Zones + dims
# Retourne un squelette (depth=3) + la liste des zones path names (znp)
# si baseNames: retourne les znp pour les bases donnees
# si familyZoneNames: retourne les znp des zones de famille donnee
# si BCType: retourne les znp des zones contenant une BC de type BCType
# ou de famille 'familySpecified:WALL'
# maxDepth: peut-etre mis a 2 ou 3 suivant l'utilisant que l'on veut faire du squelette
# readProcNode: si True, ajoute le procNode (pas lu si depth < 4)
#============================================================================
def readZoneHeaders(fileName, format=None, baseNames=None, familyZoneNames=None, BCType=None,
                    maxDepth=3, readProcNode=False):
    a = convertFile2SkeletonTree(fileName, format, maxDepth=maxDepth, maxFloatSize=6)
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
        if BCType == 'BCWall':
          families1 = PyTree.getFamilyBCNamesOfType(a, 'BCWallInviscid')
          families2 = PyTree.getFamilyBCNamesOfType(a, 'BCWallViscous*')
        znp = []
        bases = Internal.getBases(a)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                path = '/'+b[0]+'/'+z[0]
                _readPyTreeFromPaths(a, fileName, path+'/ZoneBC', format)
                nodes = Internal.getNodesFromValue(z, BCType)
                nodes += PyTree.getFamilyBCs(z, families)
                #_convert2SkeletonTree(z)
                if nodes != []: znp.append(path)
    else:
        znp = Internal.getZonePaths(a, pyCGNSLike=True)
    if readProcNode:
      if maxDepth < 2: raise ValueError('loadSkeleton: maxDepth must be >= 2 for use with readProcNode.')
      paths = []
      bases = Internal.getBases(a)
      for b in bases:
        zones = Internal.getZones(b)
        for z in zones:
          path = '/'+b[0]+'/'+z[0]+'/.Solver#Param/proc'
          paths.append(path)
        for z in zones:
          Internal._createUniqueChild(z, '.Solver#Param', 'UserDefinedData_t')
      _readPyTreeFromPaths(a, fileName, paths, format)
    return a, znp

#========================================================================
# Load par containeurs
# Load les containeurs "cont" dans a pour les znp donnes
# cont='GridCoordinates', 'FlowSolution'
#========================================================================
def _loadContainer(a, fileName, znp, cont, format=None):
    if isinstance(cont, list): conts = cont
    else: conts = [cont]
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    for p in znps:
        paths = []
        for c in conts: paths.append('%s/%s'%(p,c))
        _readPyTreeFromPaths(a, fileName, paths, format)
    return None

# variablesN = ['GridCoordinates/CoordinateX',...]
# variablesC = ['FlowSolution#Centers/Density',...]
def _loadContainerPartial(a, fileName, znp, variablesN=[], variablesC=[], format=None):
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    for p in znps:
      f = {}
      spl = p.rsplit('/',1)
      zname = spl[1]
      bname = spl[0].rsplit('/',1)[1]
      # Get loc2Glob
      zone = Internal.getNodeFromPath(a, p)
      src = None
      if zone is not None: 
        src, loc2glob = Internal.getLoc2Glob(zone)
      if src is None: continue

      j = loc2glob
      pname = src

      # Variables aux noeuds
      paths = []
      for v in variablesN: paths.append('/%s/%s/%s'%(bname, pname, v))

      DataSpaceMMRY = [[0,0,0], [1,1,1], [j[1]-j[0]+1,j[3]-j[2]+1,j[5]-j[4]+1], [1,1,1]]
      DataSpaceFILE = [[j[0]-1,j[2]-1,j[4]-1], [1,1,1], [j[1]-j[0]+1,j[3]-j[2]+1,j[5]-j[4]+1], [1,1,1]]
      DataSpaceGLOB = [[0]]
    
      for p in paths: f[p] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB

      # Variables aux centres
      paths = []
      for v in variablesC: paths.append('/%s/%s/%s'%(bname, pname, v))

      DataSpaceMMRYC = [[0,0,0], [1,1,1], [max(j[1]-j[0],1),max(j[3]-j[2],1),max(j[5]-j[4],1)], [1,1,1]]
      DataSpaceFILEC = [[j[0]-1,j[2]-1,j[4]-1], [1,1,1], [max(j[1]-j[0],1),max(j[3]-j[2],1),max(j[5]-j[4],1)], [1,1,1]]
      DataSpaceGLOBC = [[0]]

      for p in paths: f[p] = DataSpaceMMRYC+DataSpaceFILEC+DataSpaceGLOBC

      r = readNodesFromFilter(fileName, f, format)

      # Repositionne les chemins dans la zone
      for k in r:
        k2 = k.replace(pname, zname)
        spl = k2.rsplit('/',1)
        varName = spl[1]
        contName = spl[0].rsplit('/',1)[1]
        n = Internal.getNodeFromPath(a, k2)
        if n is None:
          parent = k2.rsplit('/',1)[0]
          parent = Internal.getNodeFromPath(a, parent)
          if parent is not None:
            Internal.createChild(parent, varName, 'DataArray_t', value=r[k])
            loc = Internal.getNodeFromType1(parent, 'GridLocation_t')
            if loc is None and contName+'/'+varName in variablesC:
              Internal.newGridLocation(value='CellCenter', parent=parent)
          else: print('Cannot set %s.'%k2)
        else: n[1] = r[k]
    
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

# force le proc node des zones au processeur courant
def _enforceProcNode(a):
  zones = Internal.getZones(a)
  for z in zones:
    p = Internal.createUniqueChild(z, '.Solver#Param', 'UserDefinedData_t')
    Internal.createUniqueChild(p, 'proc', 'DataArray_t', value=Cmpi.rank)
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
        elif i[3] == 'ZoneGridConnectivity_t': paths.append(p+'/'+i[0])
        if i[0] == '.Solver#define': paths.append(p+'/'+i[0])
        elif i[0] == '.Solver#ownData': paths.append(p+'/'+i[0])
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
  if paths != []: _readPyTreeFromPaths(a, fileName, paths, format, maxDepth=2, setOnlyValue=False)
  return None

# Load fully extra nodes of a tree (but no Zone or Base)
def _loadTreeExtras(a, fileName, format=None):
  # Fully load nodes of level 1 except CGNSBase_t
  paths = []
  children = a[2]
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

# Load extra data of zones (cartesianData or solver#define)
def _loadZoneExtras(a, fileName, znp, format=None, uncompress=True):
  if isinstance(znp, list): znps = znp
  else: znps = [znp]
  paths = []
  for p in znps:
    n = Internal.getNodeFromPath(a, p)
    # Level1
    for i in n[2]:
      if i[3] == 'DataArray_t': paths.append(p+'/'+i[0])
  
    # Level2
    for i in n[2]:
      if i[3] == 'UserDefinedData_t':
        for j in i[2]:
          if i[3] == 'DataArray_t': paths.append(p+'/'+i[0]+'/'+j[0])
  if paths != []: _readPyTreeFromPaths(a, fileName, paths, format)
  # Decompression cartesien evetuellement
  if uncompress:
    for p in znps:
      n = Internal.getNodeFromPath(a, p)
      for i in n[2]:
        if i[0] == 'CartesianData':
          try:
            import Compressor.PyTree as Compressor
            Compressor._uncompressCartesian(n)
          except: pass
  return None

# get variables: return a list
# this doesnt mean that all vars are defined in all zones!
def getVariables(fileName, znp, cont=None, format=None):
  if isinstance(znp, list): znps = znp
  else: znps = [znp]
  if cont is None: # check containers in files
    conts = set()
    nodes = readNodesFromPaths(fileName, znp, format, maxDepth=1, maxFloatSize=0)
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
  nodes = readNodesFromPaths(fileName, paths, format, maxDepth=1, maxFloatSize=0)
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
def isInBBox(a, fileName, format, bbox, znp):
    xmin = bbox[0]; ymin = bbox[1]; zmin = bbox[2]
    xmax = bbox[3]; ymax = bbox[4]; zmax = bbox[5]
    if isinstance(znp, list): znps = znp
    else: znps = [znp]
    out = []
    for p in znps:
       z = Internal.getNodeFromPath(a, p)
       path = ['%s/GridCoordinates/CoordinateX'%p]
       _readPyTreeFromPaths(a, fileName, path, format)
       pt = Internal.getNodeFromName2(z, 'CoordinateX')
       vmin = numpy.min(pt[1])
       vmax = numpy.max(pt[1])
       pt[1] = None
       if vmax < xmin or vmin > xmax: out.append(False); continue
       path = ['%s/GridCoordinates/CoordinateY'%p]
       _readPyTreeFromPaths(a, fileName, path, format)
       pt = Internal.getNodeFromName2(z, 'CoordinateY')
       vmin = numpy.min(pt[1])
       vmax = numpy.max(pt[1])
       pt[1] = None
       if vmax < ymin or vmin > ymax: out.append(False); continue
       path = ['%s/GridCoordinates/CoordinateZ'%p]
       _readPyTreeFromPaths(a, fileName, path, format)
       pt = Internal.getNodeFromName2(z, 'CoordinateZ')
       vmin = numpy.min(pt[1])
       vmax = numpy.max(pt[1])
       pt[1] = None
       if vmax < zmin or vmin > zmax: out.append(False); continue
       out.append(True)
    if len(out) == 1: out = out[0]
    return out

# Load only a split of zones of file
def loadAndSplit(fileName, NParts=None, noz=None, NProc=None, rank=None, variables=[]):
  if NParts is None: NParts = NProc
  if NProc is None: NProc = NParts
  import Transform.PyTree as T
  a = convertFile2SkeletonTree(fileName)

  # split on skeleton
  T._splitNParts(a, N=NParts)
  # distribute skeleton
  import Distributor2.PyTree as D2   
  D2._distribute(a, NProc)

  zones = Internal.getZones(a)
  if rank is not None:
    import Distributor2.PyTree as D2   
    D2._distribute(b, NProc)
    noz = []
    for i, z in enumerate(zones):
      p = D2.getProc(z)
      if p == rank: noz.append(i)

  f = {}
  for i in noz:
    z = zones[i]
    j = splitDict[z[0]]
    # Correction de dimension en attendant que subzone fonctionne sur un skel
    ni = j[2]-j[1]+1; nj = j[4]-j[3]+1; nk = j[6]-j[5]+1
    ni1 = max(ni-1,1); nj1 = max(nj-1,1); nk1 = max(nk-1,1)
    d = numpy.empty((3,3), numpy.int32, order='Fortran')
    d[0,0] = ni;  d[1,0] = nj;  d[2,0] = nk
    d[0,1] = ni1; d[1,1] = nj1; d[2,1] = nk1
    d[0,2] = 0;   d[1,2] = 0;   d[2,2] = 0
    z[1] = d
    # END correction
    
    zname = z[0] # zone name in b
    pname = zname.rsplit('.',1)[0] # parent name (guessed)
    (base, c) = Internal.getParentOfNode(b, z)
    bname = base[0]
    
    f = {}
    # Node fields
    path = []
    cont = Internal.getNodeFromName2(z, Internal.__GridCoordinates__)
    if cont is not None:
      for c in cont[2]:
        if c[3] == 'DataArray_t':
          path += ['/%s/%s/%s/%s'%(bname, pname,cont[0],c[0])]
    cont = Internal.getNodeFromName2(z, Internal.__FlowSolutionNodes__)
    if cont is not None:
      for c in cont[2]:
        if c[3] == 'DataArray_t' and c[0] in variables:
          path += ['/%s/%s/%s/%s'%(bname, pname,cont[0],c[0])]
    
    DataSpaceMMRY = [[0,0,0], [1,1,1], [ni,nj,nk], [1,1,1]]
    DataSpaceFILE = [[j[1]-1,j[3]-1,j[5]-1], [1,1,1], [ni,nj,nk], [1,1,1]]
    DataSpaceGLOB = [[ni,nj,nk]]
    
    for p in path:
      f[p] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB

    # Center fields
    path = []
    cont = Internal.getNodeFromName2(z, Internal.__FlowSolutionCenters__)
    if cont is not None:
      for c in cont[2]:
        if c[3] == 'DataArray_t' and 'centers:'+c[0] in variables:
          path += ['/%s/%s/%s/%s'%(bname, pname,cont[0],c[0])]
    
    DataSpaceMMRY = [[0,0,0], [1,1,1], [ni1,nj1,nk1], [1,1,1]]
    DataSpaceFILE = [[j[1]-1,j[3]-1,j[5]-1], [1,1,1], [ni1,nj1,nk1], [1,1,1]]
    DataSpaceGLOB = [[ni1,nj1,nk1]]

    for p in path:
      f[p] = DataSpaceMMRY+DataSpaceFILE+DataSpaceGLOB

    r = readNodesFromFilter(fileName, f, format)

    # Repositionne les chemins
    for k in r:
      k2 = k.replace(pname, zname)
      n = Internal.getNodeFromPath(b, k2)
      n[1] = r[k]

  b = convert2PartialTree(b)
  return b

#==========================================================
# a : must be a tree or a zone list coherent with znp
# znp: is the full path from top
def writeZones(a, fileName, znp, format=None):
  if isinstance(znp, list): znps = znp
  else: znps = [znp]
  if len(a) < 4: zones = a # suppose single zone in a list
  else:
    if a[3] == 'CGNSTree_t':
      zones = []
      for p in znps: zones.append(Internal.getNodeFromPath(a, p))
    elif a[3] == 'Zone_t': zones = [a]
    else: zones = a

  if len(zones) != len(znps): raise ValueError('writeZones: data and znp are incompatible.')
  
  paths = []
  for p in znps:
    pp = p.rsplit('/', 1)
    paths.append(pp[0])
  writeNodesFromPaths(fileName, paths, zones, format, mode=0)
  
# Write all except flowfield containers
def writeZonesWoVars(a, fileName, znp, format=None):
  if isinstance(znp, list): znps = znp
  else: znps = [znp]
  if len(a) < 4: zones = a # suppose single zone in a list
  else:
    if a[3] == 'CGNSTree_t':
      zones = []
      for p in znps: zones.append(Internal.getNodeFromPath(a, p))
    elif a[3] == 'Zone_t': zones = [a]
    else: zones = a
  
  if len(zones) != len(znps): raise ValueError('writeZones: data and znp are incompatible.')
  
  ppaths = []
  for p in znps:
    pp = p.rsplit('/', 1)
    ppaths.append(pp[0])
  writeNodesFromPaths(fileName, ppaths, zones, format, maxDepth=0, mode=0)

  nodes = []; paths = []
  # skipType=FlowSolution_t
  for c, z in enumerate(zones):
    children = z[2]
    for i in children:
      if i[3] != 'FlowSolution_t': 
        nodes.append(i)
        paths.append(znps[c])
  writeNodesFromPaths(fileName, paths, nodes, format, mode=0)

# Write zone variables only
def writeVariables(a, fileName, var, znp, format=None):
  if isinstance(znp, list): znps = znp
  else: znps = [znp]
  if isinstance(var, list): vars = var
  else: vars = [var]
  if len(a) < 4: zones = a # suppose single zone in a list
  else:
    if a[3] == 'CGNSTree_t':
      zones = []
      for p in znps: zones.append(Internal.getNodeFromPath(a, p))
    elif a[3] == 'Zone_t': zones = [a]
    else: zones = a
  
  if len(zones) != len(znps): raise ValueError('writeZones: data and znp are incompatible.')
  
  loc = []
  for v in vars:
    vs = v.split(':')
    if len(vs) == 2: 
      if vs[0] == 'centers': loc.append(Internal.__FlowSolutionCenters__)
      else: loc.append(Internal.__FlowSolutionNodes__)
    else: loc.append(Internal.__FlowSolutionNodes__)

  conts = []; cpaths = []; nodes = []; npaths = []
  for c, z in enumerate(zones):
    for i in z[2]:
      if i[3] == 'FlowSolution_t': 
        conts.append(i); cpaths.append(znps[c])
      for j in i[2]:
        if j[3] == 'GridLocation_t': 
          conts.append(j); cpaths.append(znps[c]+'/'+i[0])
    for d, v in enumerate(vars):
      ns = PyTree.getStdNodesFromName(z, v)
      if ns != []:
        npaths.append(znps[c]+'/'+loc[d])
        nodes.append(ns[0])
        
  writeNodesFromPaths(fileName, cpaths, conts, format, maxDepth=0, mode=0)
  writeNodesFromPaths(fileName, npaths, nodes, format, mode=0)

#==========================================================
class Handle:
  """Handle for partial reading."""
  def __init__(self, fileName):
    self.fileName = fileName
    self.fileVars = None # vars existing in file
    self.fileBCVars = None # BCDataSet vars existing in file
    self.znp = None # zone paths existing in file
    self.size = None # size des zones du fichier
    self.bary = None # Barycentres zones du fichier
    self.bbox = None # BBox zones du fichier
    self.hmoy = None # pas moyen des zones du fichier
    self.format = Converter.checkFileType(fileName) # Real format of file
    
  # Retourne les chemins des zones de a
  def getZonePaths(self, a):
    out = []
    bases = Internal.getBases(a)
    if len(bases) > 0:
      for b in bases:
        zones = Internal.getZones(b)
        for z in zones: out.append('/'+b[0]+'/'+z[0])
    else:
      zones = Internal.getZones(a)
      for z in zones: 
        for p in self.znp:
          r = p.rsplit('/',1)[1]
          if r == z[0]: out.append(p); break
    return out

  def setZnp(a):
    self.znp = []
    bases = Internal.getBase()
    for b in bases:
      zones = Internal.getZones(b)
      for z in zones: self.znp.append('/'+b[0]+'/'+z[0])

  # Retourne les variables du fichier
  def getVariables(self, a=None, cont=None):
    """Return the variable names contained in file."""
    if a is not None: p = self.getZonePaths(a)
    else: p = [self.znp[0]] # only first zone
    vars = getVariables(self.fileName, p, cont, self.format)
    self.fileVars = vars
    return vars

  # Retourne les variables de BCDatSet du fichier
  def getBCVariables(self, a):
    if a is not None: p = self.getZonePaths(a)
    else: p = self.znp
    vars = getBCVariables(a, self.fileName, p, None, self.format)
    self.fileBCVars = vars
    return vars

  # Charge un squelette, stocke les chemins des zones du fichier (znp)
  # Stocke le nombre de pts de chaque zone
  def loadSkeleton(self, maxDepth=3, readProcNode=False):
    """Load a skeleton tree."""
    if Cmpi.rank == 0:
      a, self.znp = readZoneHeaders(self.fileName, self.format, maxDepth=maxDepth, readProcNode=readProcNode)
      # evaluation taille des zones
      self.size = {}
      for zn in self.znp:
        z = Internal.getNodeFromPath(a, zn)
        if z[1] is not None:
          pt = z[1].ravel('k')
          if pt[2] == 0: s = pt[0]
          else: s = pt[0]*pt[1]*pt[2]
        else: s = 0
        self.size[zn] = s
    else: a = None
    a = Cmpi.bcast(a)
    return a

  # Charge le squelette, le split et conserve les infos de split
  def loadAndSplitSkeleton(self, NParts=None, NProc=Cmpi.size):
    """Load and split skeleton."""
    if Cmpi.rank == 0:
      a, self.znp = readZoneHeaders(self.fileName, self.format)
      # force loc2glob to himself from himself
      for z in Internal.getZones(a):
        dim = Internal.getZoneDim(z)
        Internal._setLoc2Glob(z, z[0], win=[1,dim[1],1,dim[2],1,dim[3]], sourceDim=dim[1:])
      import Transform.PyTree as T
      # split on skeleton
      if NParts is not None: T._splitNParts(a, N=NParts)
      else: T._splitNParts(a, N=NProc)
      if NProc is not None:
        import Distributor2.PyTree as D2   
        D2._distribute(a, NProc)
    else: a = None
    a = Cmpi.bcast(a)
    return a

  def loadAndSplit(self, NParts=None, NProc=Cmpi.size):
    """Load and split a file."""
    a = self.loadAndSplitSkeleton(NParts, NProc)
    _convert2PartialTree(a, rank=Cmpi.rank)
    varsN = ['%s/CoordinateX'%Internal.__GridCoordinates__,
             '%s/CoordinateY'%Internal.__GridCoordinates__,
             '%s/CoordinateZ'%Internal.__GridCoordinates__]
    #varsN += self.getVariables(cont=Internal.__FlowSolutionNodes__)
    #varsC = self.getVariables(cont=Internal.__FlowSolutionCenters__)
    varsC = []
    self._loadContainerPartial(a, variablesN=varsN, variablesC=varsC)
    return a
    
  def loadFromProc(self):
    """Load and distribute zones from proc node."""
    if Cmpi.rank == 0:
      t = convertFile2SkeletonTree(self.fileName, self.format, maxDepth=2, maxFloatSize=6)
      zones = Internal.getZones(t)
      paths = []
      for z in zones:
        p = Internal.getPath(t, z)+'/.Solver#Param/proc'
        paths.append(p)
      for z in zones:
        p = Internal._createUniqueChild(z, '.Solver#Param', 'UserDefinedData_t')
      _readPyTreeFromPaths(t, self.fileName, paths)
      self._loadTreeExtras(t)
    else: t = None
    t = Cmpi.bcast(t)
    paths = []
    zones = Internal.getZones(t)
    for z in zones:
      proc = Internal.getNodeFromName2(z, 'proc')
      if Internal.getValue(proc) == Cmpi.rank:
        paths.append(Internal.getPath(t, z))
      else:
        Internal._rmNode(t, z)
    if paths != []: _readPyTreeFromPaths(t, self.fileName, paths)
    return t

  def loadAndDistribute(self, useCom=None):
    """Load and distribute zones."""
    if Cmpi.rank == 0:
      
      if useCom == 'match':
        t = convertFile2SkeletonTree(self.fileName, self.format, maxDepth=3, maxFloatSize=6)
        paths = []
        zones = Internal.getZones(t)
        for z in zones:
          pz = Internal.getPath(t, z)
          gcs = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
          for gc in gcs:
            p = pz+'/%s'%gc[0]
            paths.append(p)
        if paths != []: _readPyTreeFromPaths(t, self.fileName, paths, self.format)
      else:
        t = convertFile2SkeletonTree(self.fileName, self.format, maxDepth=2, maxFloatSize=6)
        paths = []
        zones = Internal.getZones(t)
        for z in zones:
          p = Internal.getPath(t, z)+'/ZoneType'
          paths.append(p)
        _readPyTreeFromPaths(t, self.fileName, paths, self.format)  
      self._loadTreeExtras(t)
      import Distributor2.PyTree as D2
      D2._distribute(t, Cmpi.size, useCom=useCom)
    else: t = None
    t = Cmpi.bcast(t)
    paths = []
    zones = Internal.getZones(t)
    for z in zones:
      proc = Internal.getNodeFromName2(z, 'proc')
      if Internal.getValue(proc) == Cmpi.rank:
        paths.append(Internal.getPath(t, z))
      else:
        Internal._rmNode(t, z)
    if paths != []: _readPyTreeFromPaths(t, self.fileName, paths, self.format)
    _enforceProcNode(t)
    return t

  # Calcul et stocke des infos geometriques sur les zones
  def geomProp(self, znp=None):
    """Store some zone properties."""
    if znp is None: znp = self.znp

    if self.bbox is None: self.bbox = {}
    if self.bary is None: self.bary = {}
    if self.hmoy is None: self.hmoy = {}
    for p in znp:
       if p in self.bbox: continue
       xc = 0.; yc = 0.; zc = 0.
       bbox = [0.,0.,0.,0.,0.,0.]
       hmoy = 0.; hmax = 0.; hmin = 0.
       # Coordinate X
       path = '%s/GridCoordinates/CoordinateX'%p
       n = readNodesFromPaths(self.fileName, path, maxDepth=0)
       vmin = numpy.min(n[1])
       vmax = numpy.max(n[1])
       bbox[0] = vmin; bbox[3] = vmax
       xc = 0.5*(vmin+vmax)
       # Coordinate Y
       path = '%s/GridCoordinates/CoordinateY'%p
       n = readNodesFromPaths(self.fileName, path, maxDepth=0)
       vmin = numpy.min(n[1])
       vmax = numpy.max(n[1])
       bbox[1] = vmin; bbox[4] = vmax
       yc = 0.5*(vmin+vmax)
       # Coordinate Z
       path = '%s/GridCoordinates/CoordinateZ'%p
       n = readNodesFromPaths(self.fileName, path, maxDepth=0)
       vmin = numpy.min(n[1])
       vmax = numpy.max(n[1])
       bbox[2] = vmin; bbox[5] = vmax
       zc = 0.5*(vmin+vmax)
       self.bary[p] = (xc,yc,zc)
       self.bbox[p] = bbox
       dhx = bbox[3]-bbox[0]
       dhy = bbox[4]-bbox[1]
       dhz = bbox[5]-bbox[2]
       npts = self.size[p]**0.33
       self.hmoy[p] = (dhx+dhy+dhz)*0.33*npts

  # Charge les Grid coordinates + grid connectivity + BC
  # pour toutes les zones de a
  def _loadGCBC(self, a):
    self._loadContainer(a, 'GridCoordinates')
    self._loadConnectivity(a)
    self._loadZoneBCs(a)
    Internal._fixNGon(a)
    return None

  # Charge tous les noeuds extra d'un arbre a
  def _loadTreeExtras(self, a):
    _loadTreeExtras(a, self.fileName, self.format)
    return None

  # Charge tous les noeuds extra des zones
  def _loadZoneExtras(self, a):
    _loadZoneExtras(a, self.fileName, self.format)
    return None

  # Charge completement toutes les zones de a ou de chemin fourni
  def loadZones(self, a, znp=None):
    """Fully load zones."""
    b = Internal.copyRef(a)
    self._loadZones(b, znp)
    return b

  def _loadZones(self, a, znp=None):
    """Fully load zones."""
    if znp is None: znp = self.getZonePaths(a)
    _loadZones(a, self.fileName, znp, self.format)
    return None

  # Charge les grid coordinates + grid connectivity + BC
  # pour les zones specifiees
  def loadZonesWoVars(self, a, znp=None, bbox=None):
    """Load zones without loading variables."""
    b = Internal.copyRef(a)
    self._loadZonesWoVars(b, znp, bbox)
    return b

  def _loadZonesWoVars(self, a, znp=None, bbox=None):
    if znp is None: znp = self.getZonePaths(a)
    if bbox is None:
      # Read paths as skeletons
      _readPyTreeFromPaths(a, self.fileName, znp, self.format, maxFloatSize=0)
      _loadContainer(a, self.fileName, znp, 'GridCoordinates', self.format)
      _loadConnectivity(a, self.fileName, znp, self.format)
      _loadZoneBCs(a, self.fileName, znp, self.format)
      _loadZoneExtras(a, self.fileName, znp, self.format)
      for zp in znp:
        _convert2PartialTree(Internal.getNodeFromPath(a, zp))
    else:
      if self.bbox is None: self.geomProp()
      for zp in znp:
        zbb = self.bbox[zp]
        #print zbb[3],'<', bbox[0]
        #print zbb[0],'>', bbox[3]
        #print zbb[4],'<', bbox[1]
        #print zbb[1],'>', bbox[4]
        #print zbb[5],'<', bbox[2]
        #print zbb[2],'>', bbox[5]
        if zbb[3] >= bbox[0] and zbb[0] <= bbox[3] and zbb[4] >= bbox[1] and zbb[1] <= bbox[4] and zbb[5] >= bbox[2] and zbb[2] <= bbox[5]:
          print('loading: %s'%zp)
          _readPyTreeFromPaths(a, self.fileName, [zp], self.format, maxFloatSize=0)
          _loadContainer(a, self.fileName, [zp], 'GridCoordinates', self.format)
          _loadConnectivity(a, self.fileName, [zp], self.format)
          _loadZoneBCs(a, self.fileName, [zp], self.format)
          _convert2PartialTree(Internal.getNodeFromPath(a, zp))
    return None

  # Charge toutes les BCs (avec BCDataSet) des zones de a  
  def _loadZoneBCs(self, a, znp=None):
    if znp is None: znp = self.getZonePaths(a)
    _loadZoneBCs(a, self.fileName, znp, self.format)
    return None

  # Charge toutes les BCs (sans BCDataSet) des zones de a
  def _loadZoneBCsWoData(self, a, znp=None):
    if znp is None: znp = self.getZonePaths(a)
    _loadZoneBCsWoData(a, self.fileName, znp, self.format)
    return None

  # Charge la connectivite pour toutes les zones de a
  def _loadConnectivity(self, a, znp=None):
    if znp is None: znp = self.getZonePaths(a)
    _loadConnectivity(a, self.fileName, znp, self.format)
    return None

  # Charge le container "cont" pour toutes les zones de a
  def _loadContainer(self, a, cont, znp=None):
    if znp is None: znp = self.getZonePaths(a)
    _loadContainer(a, self.fileName, znp, cont, self.format)
    return None

  # Charge le container "cont" en partiel (si la zone a loc2glob)
  def _loadContainerPartial(self, a, variablesN=[], variablesC=[], znp=None):
    if znp is None: znp = self.getZonePaths(a)
    _loadContainerPartial(a, self.fileName, znp, variablesN, variablesC, self.format)
    return None
   
  def loadVariables(self, a, var, znp=None):
    """Load specified variables."""
    b = Internal.copyRef(a)
    self._loadVariables(b, var, znp)
    return b

  # Charge la ou les variables "var" pour toutes les zones de a
  def _loadVariables(self, a, var, znp=None):
    """Load specified variables."""
    if znp is None: znp = self.getZonePaths(a)
    _loadVariables(a, self.fileName, znp, var, self.format)
    return None
  
  def isInBBox(self, a, bbox, znp=None):
    if znp is None: znp = self.getZonePaths(a)
    return isInBBox(a, self.fileName, self.format, bbox, znp)

  # Ecrit des zones
  def writeZones(self, a, znp=None):
    """Write specified zones."""
    if znp is None: znp = self.getZonePaths(a)
    writeZones(a, self.fileName, znp, self.format)

  # Ecrit des zones sans les FlowSolution_t
  def writeZonesWoVars(self, a, znp=None):
    """Write specified zones without variables."""
    if znp is None: znp = self.getZonePaths(a)
    writeZonesWoVars(a, self.fileName, znp, self.format)
  
  # Ecrit des variables
  def writeVariables(self, a, var, znp=None):
    """Write specified variables."""
    if znp is None: znp = self.getZonePaths(a)  
    writeVariables(a, self.fileName, var, znp, self.format)

  # save zones + field
  # mode 0: parallele, 1: écriture chacun son tour
  def save(self, a, fileName=None, mode=0):
    if fileName is None: fileName = self.fileName
    if mode == 0: self.writeZones(a)
    else: Cmpi.convertPyTree2File(a, self.fileName, self.format, ignoreProcNode=True)
    
  def mergeAndSave(self, a, fileName=None):
    if fileName is None: fileName = self.fileName
    zones = Internal.getZones(a)
    # write skeleton data corresponding to a
    for z in zones:
      p = Internal.getPath(a, z)
      
      