"""OpenCascade definition module (pyTree).
"""
try:
    import OCC
    import OCC.occ as occ
    import Converter
    import Converter.PyTree as C
    import Generator.PyTree as G
    import Converter.Internal as Internal
    import Converter.Mpi as Cmpi
except ImportError:
  raise ImportError("OCC.PyTree: requires Converter module.")

__version__ = OCC.__version__
import numpy

#==============================================================================
# -- convertCAD2PyTree --
#==============================================================================
def convertCAD2PyTree(fileName, format=None, h=0., chordal_err=0., 
  growth_ratio=0., merge_tol=-1, algo=1, join=True):
  """Convert a CAD (IGES or STEP) file to pyTree.
  Usage: convertCAD2PyTree(fileName, options)"""
  a = OCC.convertCAD2Arrays(fileName, format, h, chordal_err, growth_ratio, merge_tol, algo, join)
  
  t = C.newPyTree([])
  base1 = False; base2 = False; base3 = False; base = 1
  
  for i in a:
    if len(i) == 5: # Structure
      if i[3] == 1 and i[4] == 1:
        if not base1:
          t = C.addBase2PyTree(t, 'Base1', 1); base1 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base1][2].append(z)
      elif i[4] == 1:
        if not base2:
          t = C.addBase2PyTree(t, 'Base2', 2); base2 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__) 
        t[2][base2][2].append(z)
      else:
        if not base3:
          t = C.addBase2PyTree(t, 'Base', 3); base3 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base3][2].append(z)
    else: # non structure
      if i[3] == 'BAR':
        if not base1:
          t = C.addBase2PyTree(t, 'Base1', 1); base1 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base1][2].append(z)
      elif i[3] == 'TRI' or i[3] == 'QUAD':
        if not base2:
          t = C.addBase2PyTree(t, 'Base2', 2); base2 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base2][2].append(z)
      else:
        if not base3:
          t = C.addBase2PyTree(t, 'Base', 3); base3 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base3][2].append(z)

  Internal._correctPyTree(t, level=2) # force unique name
  Internal._correctPyTree(t, level=7) # create familyNames
  return t

#================================================================================
def meshSTRUCT(fileName, format="fmt_step", N=11):
  """Return a STRUCT discretisation of CAD."""
  hook = OCC.occ.readCAD(fileName, format)
  return meshSTRUCT__(hook, N) 

def meshSTRUCT__(hook, N=11, faceSubset=None, linkFaceNo=None):
  """Return a STRUCT discretisation of CAD."""
  faceNoA = []
  a = OCC.meshSTRUCT__(hook, N, faceSubset, faceNoA)
  out = []
  for c, i in enumerate(a):
    z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    out.append(z)
    if linkFaceNo is not None: linkFaceNo[z[0]] = faceNoA[c]
  return out

#================================================================================
def meshTRI(fileName, format="fmt_step", N=11, hmax=-1., order=1):
  """Return a TRI discretisation of CAD."""
  hook = OCC.occ.readCAD(fileName, format)
  return meshTRI__(hook, N, hmax, order) 

def meshTRI__(hook, N=11, hmax=-1., order=1, faceSubset=None, linkFaceNo=None):
  """Return a TRI discretisation of CAD."""
  faceNoA = []
  a = OCC.meshTRI__(hook, N, hmax, order, faceSubset, faceNoA)
  out = []
  for c, i in enumerate(a):
    z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    out.append(z)
    if linkFaceNo is not None: linkFaceNo[z[0]] = faceNoA[c]
  return out

def meshTRIHO(fileName, format="fmt_step", N=11):
  """Return a TRI HO discretisation of CAD."""
  a = OCC.meshTRIHO(fileName, format, N)
  out = []
  for i in a:
    z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    out.append(z)
  return out

#================================================================================
def meshQUAD(fileName, format="fmt_step", N=11, order=1):
  """Return a QUAD discretisation of CAD."""
  hook = OCC.occ.readCAD(fileName, format)
  return meshQUAD__(hook, N, order)

def meshQUAD__(hook, N=11, order=1, faceSubset=None, linkFaceNo=None):
  """Return a QUAD discretisation of CAD."""
  faceNoA = []
  a = OCC.meshQUAD__(hook, N, order, faceSubset, faceNoA)
  out = []
  for c, i in enumerate(a):
    z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    out.append(z)
    if linkFaceNo is not None: linkFaceNo[z[0]] = faceNoA[c]
  return out

def meshQUADHO(fileName, format="fmt_step", N=11):
  """Return a QUAD HO discretisation of CAD."""
  hook = OCC.occ.readCAD(fileName, format)
  return meshQUADHO__(hook, N)

def meshQUADHO__(hook, N=11, faceSubset=None, linkFaceNo=None):
  """Return a QUAD HO discretisation of CAD."""
  faceNoA = []
  a = OCC.meshQUADHO__(hook, N, faceSubset, faceNoA)
  out = []
  for c, i in enumerate(a):
    z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    out.append(z)
    if linkFaceNo is not None: linkFaceNo[z[0]] = faceNoA[c]
  return out

#===========================================================================
class Edge:
  """CAD Edge."""
  def __init__(self, i, cad):
    self.number = i # no in CAD edge list
    self.name = 'XXX' # CAD edge name
    self.hook = None # hook on OCC TOPODS::edge
    self.cad = cad # master CAD object

  def valueAt(self, distribution):
    """Evaluate edge at given parameters."""
    return self.cad.evalEdge(self.number, distribution)

  def _projectOn(self, z):
    """Project z on edge."""
    a = C.getFields(Internal.__GridCoordinates__, z, api=2)
    for i in a:
      self.cad._projectOnEdges(i, [self.number])
    return None

class Face:
  """CAD Face."""
  def __init__(self, i, cad):
    self.number = i # no in CAD face list
    self.name = 'XXX' # CAD face name
    self.hook = None # hook on OCC TOPODS::face
    self.cad = cad # master CAD object

  def valueAt(self, distribution):
    """Evaluate face at given parameters."""
    return self.cad.evalFace(self.number, distribution)

  def _projectOn(self, z):
    """Project z on face."""
    a = C.getFields(Internal.__GridCoordinates__, z, api=2)
    for i in a:
      self.cad._projectOnFaces(i, [self.number])
    return None

class CAD:
  """CAD top tree."""
  def __init__(self, fileName, format='fmt_iges'):
    self.fileName = fileName
    self.format = format
    self.hook = None # hook on OCC tree
    self.faces = [] # list of CAD faces (class) 
    self.edges = [] # list of CAD edges (class)

    self.zones = [] # associated discretization (list of zones)
    self.linkFaceNo = {} # association zone Name -> CAD face no
    self.linkEdgeNo = {} # association zone Name -> CAD edge no

    # read CAD
    self.hook = OCC.occ.readCAD(fileName, format)
    nbfaces = OCC.occ.getNbFaces(self.hook)
    for i in range(nbfaces): self.faces.append(Face(i+1, self))
    nbedges = OCC.occ.getNbEdges(self.hook)
    for i in range(nbedges): self.edges.append(Edge(i+1, self))
  
  def evalFace(self, face, distribution):
    """Evaluate face at given parameters."""
    if isinstance(face, int): no = face
    else: no = face.number
    if isinstance(distribution, tuple): 
      d = Converter.array('x,y,z', 1,1,1)
      d[1][0,0] = distribution[0]
      d[1][1,0] = distribution[1]
      d[1][2,0] = 0.      
    else:
      d = C.getFields(Internal.__GridCoordinates__, distribution, api=2)[0]
    m = OCC.occ.evalFace(self.hook, d, no)
    z = Internal.createZoneNode(C.getZoneName('Face'), m, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    return z

  def evalEdge(self, edge, distribution):
    """Evaluate edge at given parameters."""
    if isinstance(edge, int): no = edge
    else: no = edge.number
    if isinstance(distribution, float):
      d = Converter.array('x,y,z', 1,1,1)
      d[1][0,0] = distribution
      d[1][1,0] = 0.
      d[1][2,0] = 0.
    else:
      d = C.getFields(Internal.__GridCoordinates__, distribution, api=2)[0]
    m = OCC.occ.evalEdge(self.hook, d, no)
    z = Internal.createZoneNode(C.getZoneName('Edge'), m, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    return z

  def _project(self, z, faceList=None):
    """Project z on CAD."""
    if faceList is not None:
        out = []
        for f in faceList:
          if isinstance(f, int): out.append(f)
          else: out.append(f.number)
    else: out = None
    a = C.getFields(Internal.__GridCoordinates__, z, api=2)
    for i in a: OCC.occ.projectOnFaces(self.hook, i, out)
    return None

  def project(self, z, faceList=None):
    """Project z on CAD."""
    zp = Internal.copyTree(z)
    self._project(zp, faceList)
    return zp

  # faceList ne marche pas encore
  def mesh(self, mtype='STRUCT', N=11, hmax=-1.):
    """Mesh CAD with given type."""
    if mtype == 'STRUCT':
      zones = meshSTRUCT__(self.hook, N, None, self.linkFaceNo)
    elif mtype == 'TRI':
      zones = meshTRI__(self.hook, N, hmax, None, self.linkFaceNo)
    elif mtype == 'QUAD':
      zones = meshQUAD__(self.hook, N, 1, None, self.linkFaceNo)
    elif mtype == 'TRIHO':
      zones = meshTRI__(self.hook, N, -1., 2, None, self.linkFaceNo)
    elif mtype == 'QUADHO':
      zones = meshQUAD__(self.hook, N, 2, None, self.linkFaceNo)
    else: raise ValueError("mesh: not a valid meshing type.")
    self.zones += zones
    return zones

  def getLinkFace(self, zone):
    """Return the faces linked to zone."""
    if isinstance(zone, str): name = zone
    else: name = zone[0]
    no = self.linkFaceNo[name]
    return self.faces[no-1]

#========================
#=== nouvelle vision ====
#========================
def readCAD(fileName, format='fmt_step'):
  """Read CAD and return a CAD hook."""
  return OCC.occ.readCAD(fileName, format)

def writeCAD(hook, fileName, format='fmt_step'):
  """Write CAD file from CAD hook."""
  OCC.occ.writeCAD(hook, fileName, format)
  return None

def freeHook(hook):
  """Free hook."""
  OCC.occ.freeHook(hook)
  return None

def _linkCAD2Tree(hook, t):
  """Put hook in CAD/hook for each zone."""
  zones = Internal.getZones(t)
  for z in zones:
    r = Internal.getNodeFromName1(z, "CAD")
    if r is not None:
      l = Internal.getNodeFromName1(r, "hook")
      if l is not None: l[1] = hook
  return None

def _unlinkCAD2Tree(t):
  """Suppress hook in CAD/hook for each zone."""
  zones = Internal.getZones(t)
  for z in zones:
    r = Internal.getNodeFromName1(z, "CAD")
    if r is not None:
      l = Internal.getNodeFromName1(r, "hook")
      if l is not None: l[1] = 0
  return None

# Get the first tree for structured CAD meshing
def getTree(hook, N=11, hmax=-1, hausd=-1.):
  """Get a first TRI meshed tree linked to CAD."""
  
  t = C.newPyTree(['EDGES', 'FACES'])

  # Add CAD top container containing the CAD file name
  fileName, fileFmt = OCC.occ.getFileAndFormat(hook)
  _setCADcontainer(t, fileName, fileFmt, hmax, hausd)

  # Edges
  if hmax > 0.: edges = OCC.occ.meshGlobalEdges1(hook, hmax)
  elif hausd > 0.: edges = OCC.occ.meshGlobalEdges4(hook, hausd)
  else: edges = OCC.occ.meshGlobalEdges2(hook, N)

  b = Internal.getNodeFromName1(t, 'EDGES')
  for c, e in enumerate(edges):
    z = Internal.createZoneNode('edge%03d'%(c+1), e, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    # Conserve name, type et no de l'edge dans la CAD
    r = Internal.createChild(z, "CAD", "UserDefinedData_t")
    Internal._createChild(r, "name", "DataArray_t", value="edge%03d"%(c+1))
    Internal._createChild(r, "type", "DataArray_t", value="edge")
    Internal._createChild(r, "no", "DataArray_t", value=(c+1))
    b[2].append(z)

  # Faces
  b = Internal.getNodeFromName1(t, 'FACES')
  faceNo = []
  #m = OCC.meshTRI__(hook, N=N, hmax=hmax, hausd=hausd, faceNo=faceNo)
  m = OCC.meshSTRUCT__(hook, N=N, faceNo=faceNo)

  for c, f in enumerate(m):
    noface = faceNo[c]
    z = Internal.createZoneNode(C.getZoneName('face%03d'%noface), f, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    edgeNo = OCC.occ.getEdgeNoByFace(hook, noface)
    # conserve name, type
    r = Internal.createChild(z, "CAD", "UserDefinedData_t")
    Internal._createChild(r, "name", "DataArray_t", value="face%03d"%noface)
    Internal._createChild(r, "type", "DataArray_t", value="face")
    Internal._createChild(r, "no", "DataArray_t", value=noface)
    Internal._createChild(r, "edgeList", "DataArray_t", value=edgeNo)
    Internal._createChild(r, "hsize", "DataArray_t", value=[0.])
    b[2].append(z)

  _updateEdgesFaceList__(t)
  _setLonelyEdgesColor(t)

  return t

# remesh tree from edges with new ue from EDGES
def remeshTreeFromEdges(hook, tp):

  t = C.newPyTree(['EDGES', 'FACES'])

  # Edges
  b = Internal.getNodeFromName1(tp, 'EDGES')
  prevEdges = Internal.getZones(b)
  arrays = C.getAllFields(prevEdges, 'nodes', api=3)
  edges = OCC.occ.meshGlobalEdges3(hook, arrays)

  b = Internal.getNodeFromName1(t, 'EDGES')
  for c, e in enumerate(edges):
    z = Internal.createZoneNode('edge%03d'%(c+1), e, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    # Conserve name, type et no de l'edge dans la CAD
    r = Internal.createChild(z, "CAD", "UserDefinedData_t")
    Internal._createChild(r, "name", "DataArray_t", value="edge%03d"%(c+1))
    Internal._createChild(r, "type", "DataArray_t", value="edge")
    Internal._createChild(r, "no", "DataArray_t", value=(c+1))
    b[2].append(z)

  # Faces
  b = Internal.getNodeFromName1(t, 'FACES')
  faceNo = []
  m = OCC.meshTRIU__(hook, arrays, faceNo=faceNo)
  for c, f in enumerate(m):
    noface = faceNo[c]
    z = Internal.createZoneNode(C.getZoneName('face%03d'%noface), f, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    edgeNo = OCC.occ.getEdgeNoByFace(hook, noface)
    # conserve name, type
    r = Internal.createChild(z, "CAD", "UserDefinedData_t")
    Internal._createChild(r, "name", "DataArray_t", value="face%03d"%noface)
    Internal._createChild(r, "type", "DataArray_t", value="face")
    Internal._createChild(r, "no", "DataArray_t", value=noface)
    Internal._createChild(r, "edgeList", "DataArray_t", value=edgeNo)
    b[2].append(z)

  _updateEdgesFaceList__(t)
  _setLonelyEdgesColor(t)

  return t

#====================================================================================
# ULTIMATE FUNCTIONS
#====================================================================================

# return the cad no of entity (edge or face)
def getNo(e):
    cad = Internal.getNodeFromName1(e, 'CAD')
    no = Internal.getNodeFromName1(cad, 'no')
    no = Internal.getValue(no)
    return no

# return the position of entities in base baseName by number
def getPos(t, baseName=None):
  pos = {}; posi = {}
  if baseName is not None:
    b = Internal.getNodeFromName1(t, baseName)
  else: b = t # suppose t is already chosen base
  for c, e in enumerate(b[2]):
    cad = Internal.getNodeFromName1(e, 'CAD')
    if cad is not None: # this is a CAD edge or face
      no = Internal.getNodeFromName1(cad, 'no')
      no = Internal.getValue(no)
      pos[no] = c
      posi[c] = no
  return pos, posi

# return the position of edges in EDGES base
def getPosEdges(t):
  return getPos(t, 'EDGES')

# return the position of faces in FACES base
def getPosFaces(t):
  return getPos(t, 'FACES')

# return the position of edges/faces in t
# return pose, posf, posei (reverse), posfi (reverse)
def getAllPos(t):
  pose, posei = getPos(t, 'EDGES')
  posf, posfi = getPos(t, 'FACES')
  return [pose, posf, posei, posfi]

#==================================================
# get the first tree from CAD - mesh TRI the CAD
# IN: hook: hook on CAD
# IN: hmax: hmax
# IN: hausd: hausd deflection
# IN: faceList: si fourni, ne maille que ces faces
# OUT: meshed CAD with CAD links
#=================================================
def getFirstTree(hook, hmax=-1., hausd=-1., faceList=None):
  """Get a first TRI meshed tree linked to CAD."""
  
  t = C.newPyTree(['EDGES', 'FACES'])

  # Add CAD top container containing the CAD file name
  fileName, fileFmt = OCC.occ.getFileAndFormat(hook)
  _setCADcontainer(t, fileName, fileFmt, hmax, hausd)

  # - Edges -
  edges = OCC.meshAllEdges(hook, hmax, hausd, -1)

  b = Internal.getNodeFromName1(t, 'EDGES')
  for c, e in enumerate(edges):
    z = Internal.createZoneNode('edge%03d'%(c+1), e, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    # Conserve hook, name, type et no de l'edge dans la CAD
    r = Internal.createChild(z, "CAD", "UserDefinedData_t")
    Internal._createChild(r, "name", "DataArray_t", value="edge%03d"%(c+1))
    Internal._createChild(r, "type", "DataArray_t", value="edge")
    Internal._createChild(r, "no", "DataArray_t", value=(c+1))
    #Internal._createChild(r, "hook", "UserDefinedData_t", value=hook)
    b[2].append(z)

  # - Faces -
  b = Internal.getNodeFromName1(t, 'FACES')
  nbFaces = occ.getNbFaces(hook)
  # distribution parallele (CAD already split)
  if faceList is None:
    N = nbFaces // Cmpi.size
    nstart = Cmpi.rank*N
    nend = nstart+N
    if Cmpi.rank == Cmpi.size-1: nend = nbFaces
    faceList = range(nstart+1, nend+1)

  if hausd < 0:
    hList = [(hmax,hmax,hausd)]*len(faceList)
  else:
    hList = [(hmax*0.8,hmax*1.2,hausd)]*len(faceList)

  faces = OCC.meshAllFacesTri(hook, edges, True, faceList, hList)
  
  for c, f in enumerate(faces):
    if f is None: continue # Failed face
    noface = faceList[c]
    z = Internal.createZoneNode('face%03d'%(noface), f, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    edgeNo = OCC.occ.getEdgeNoByFace(hook, noface)
    # conserve hook, name, type
    r = Internal.createChild(z, "CAD", "UserDefinedData_t")
    Internal._createChild(r, "name", "DataArray_t", value="face%03d"%(noface))
    Internal._createChild(r, "type", "DataArray_t", value="face")
    Internal._createChild(r, "no", "DataArray_t", value=noface)
    Internal._createChild(r, "edgeList", "DataArray_t", value=edgeNo)
    Internal._createChild(r, "hsize", "DataArray_t", value=hList[c])
    #Internal._createChild(r, "hook", "UserDefinedData_t", value=hook)
    b[2].append(z)

  _updateEdgesFaceList__(t)
  _setLonelyEdgesColor(t)

  return t

# the first version of parallel CAD split and TRI meshing
def getFirstTreePara(hook, area, hmax=-1., hausd=-1.):
  import Distributor2
  import Distributor2.PyTree as D2

  # split CAD with max area
  OCC.occ.splitFaces(hook, area)

  # write split CAD
  #OCC.occ.writeCAD(hook, "cube_split.step", "fmt_step")

  # distribute faces
  nfaces = OCC.occ.getNbFaces(hook)

  arrays = []; weights = []
  for i in range(nfaces):
    area = OCC.occ.getFaceArea(hook, i+1)
    arrays.append(['x',None,1,1,1])
    weights.append(area)

  out = Distributor2.distribute(arrays, weight=weights, NProc=Cmpi.size)
  dis = out['distrib']
  if Cmpi.rank == 0: print(out['varMax'])

  faceList = [] # list of face to mesh on this proc
  for i in range(nfaces):
    if dis[i] == Cmpi.rank: faceList.append(i+1)

  #print(Cmpi.rank, faceList)
  t = getFirstTree(hook, hmax, hausd, faceList=faceList)
  D2._addProcNode(t, Cmpi.rank)
  return t

#================================================
# remesh faces from an external modified edges
# edges is a list of externally remeshed edges
#================================================
def _remeshTreeFromEdges(hook, t, edges):
  
  # find impacted faces by edges
  faceList = set()
  for edge in edges:
    cad = Internal.getNodeFromName1(edge, 'CAD')
    facel = Internal.getNodeFromName1(cad, 'faceList')
    if facel is not None:
      facel = facel[1]
      facel = list(facel)
      faceList.update(facel)
  faceList = list(faceList)
  
  # build hList from CAD/hsize
  hList = []
  b = Internal.getNodeFromName1(t, 'FACES')
  if b is None: faceList = [] # forced when no FACES
  be = Internal.getNodeFromName1(t, 'EDGES')
  for f in faceList:
    z = b[2][f-1]
    CAD = Internal.getNodeFromName1(z, 'CAD')
    hsize = Internal.getNodeFromName1(CAD, 'hsize')
    hsize = hsize[1]

    # modify hmax/hmin from edge sizes
    edgeList = Internal.getNodeFromName1(CAD, 'edgeList')
    edgeList = edgeList[1]
    fedges = []
    for e in edgeList:
      ze = be[2][e-1]
      fedges.append(ze)
    a = G.getMaxLength(fedges)
    hmine = C.getMinValue(a, 'centers:MaxLength')
    hmaxe = C.getMaxValue(a, 'centers:MaxLength')
    hausde = hsize[2]
    #print("hsize=",hmine,hmaxe,hausde)
    hsize = ( min(hmine, hsize[0]), max(hmaxe, hsize[1]), min(hausde, hsize[2]) )
    hList.append(hsize)
    
  # get dedges (all CAD edges - suppose CAD order)
  b = Internal.getNodeFromName1(t, 'EDGES')
  dedges = []
  for e in Internal.getZones(b):
    dedges.append(C.getFields([Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__], e, api=2)[0])
  
  # set edge in dedges and in t
  for edge in edges:
    edgeno = getNo(edge)
    aedge = C.getFields([Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__], edge, api=2)[0]
    e = occ.meshOneEdge(hook, edgeno, -1, -1, -1, aedge)
    dedges[edgeno-1] = e
    cad = Internal.getNodeFromName1(edge, 'CAD')
    render = Internal.getNodeFromName1(edge, '.RenderInfo')
    z = Internal.createZoneNode('edge%03d'%(edgeno), e, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    z[2].append(cad)
    if render is not None: z[2].append(render)
    b[2][edgeno-1] = z

  # eval the impacted faces
  faces = OCC.meshAllFacesTri(hook, dedges, metric=True, faceList=faceList, hList=hList)
  
  # replace faces in t
  b = Internal.getNodeFromName1(t, 'FACES')
  if b is not None: pos, posi = getPosFaces(t)
  for c, f in enumerate(faceList):
    cd = pos[f]
    zp = b[2][cd]
    cad = Internal.getNodeFromName1(zp, 'CAD')
    noface = getNo(zp)
    z = Internal.createZoneNode('face%03d'%(noface), faces[c], [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    z[2].append(cad)
    b[2][cd] = z

  return None

# modify hsize for faces from hlist
# IN: faceList: liste d'entiers start 1
# IN: hList: liste de (hmin,hmax,hausd)
def _modifyHSizeForFaces(t, faceList, hList):
  b = Internal.getNodeFromName1(t, 'FACES')
  for c, i in enumerate(faceList):
    hsize = hList[c]
    face = b[2][i-1]
    CAD = Internal.getNodeFromName1(face, 'CAD')
    node = Internal.getNodeFromName1(CAD, 'hsize')
    Internal._setValue(node, hsize)
  return None

# remesh from face hmin/hmax/hausd
def _remeshTreeFromFaces(hook, t, faceList, hList):
  _modifyHSizeForFaces(t, faceList, hList)

  # get dedges (all CAD edges - suppose CAD order)
  b = Internal.getNodeFromName1(t, 'EDGES')
  dedges = []
  for e in Internal.getZones(b):
    dedges.append(C.getFields([Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__], e, api=2)[0])
  
  # eval the impacted faces
  faces = OCC.meshAllFacesTri(hook, dedges, metric=True, faceList=faceList, hList=hList)
  
  # replace faces in t
  pos, posi = getPosFaces(t)
  b = Internal.getNodeFromName1(t, 'FACES')
  for c, f in enumerate(faceList):
    cd = pos[f]
    zp = b[2][cd]
    cad = Internal.getNodeFromName1(zp, 'CAD')
    noface = getNo(zp)
    z = Internal.createZoneNode('face%03d'%(noface), faces[c], [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    z[2].append(cad)
    b[2][cd] = z

  return None

# add CAD container to tree (file and format from hook)
def _setCADcontainer(t, fileName, fileFmt, hmax, hausd):
  CAD = Internal.getNodeFromName1(t, 'CAD')
  if CAD is None:
    CAD = Internal.createChild(t, 'CAD', 'UserDefinedData_t')
  if fileName is not None: Internal._createUniqueChild(CAD, 'file', 'DataArray_t', value=fileName)
  if fileFmt is not None: Internal._createUniqueChild(CAD, 'format', 'DataArray_t', value=fileFmt)
  if hmax is not None: Internal._createUniqueChild(CAD, 'hsize', 'DataArray_t', value=hmax)
  if hausd is not None: Internal._createUniqueChild(CAD, 'hausd', 'DataArray_t', value=hausd)
  return None

# mesh all edges
def _meshAllEdges(hook, t, hmax=-1, hausd=-1, N=-1, edgeList=None):
  
  if edgeList is None: 
    nbEdges = occ.getNbEdges(hook)
    edgeList = range(1, nbEdges+1)

  edges = OCC.meshAllEdges(hook, hmax, hausd, N)
  b = Internal.getNodeFromName1(t, 'EDGES')
  if b is None: b = Internal.newCGNSBase('EDGES', parent=t)

  for c, e in enumerate(edges):
    z = Internal.createZoneNode('edge%03d'%(c+1), e, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    # Conserve hook, name, type et no de l'edge dans la CAD
    r = Internal.createChild(z, "CAD", "UserDefinedData_t")
    Internal._createChild(r, "name", "DataArray_t", value="edge%03d"%(c+1))
    Internal._createChild(r, "type", "DataArray_t", value="edge")
    Internal._createChild(r, "no", "DataArray_t", value=(c+1))
    b[2].append(z)

  _setCADcontainer(t, None, None, hmax, hausd)
  return None

# remesh all CAD edges with odd number of points
def _remeshAllEdgesOdd(hook, t):
  """Remesh all edges to have an odd number of points."""
  import Geom.PyTree as D
  b = Internal.getNodeFromName1(t, 'EDGES')
  for edge in Internal.getZones(b):
    npts = C.getNPts(edge)
    if npts//2-0.5*npts == 0:
      factor = (npts*1.)/(npts-1.)
      G._refine(edge, factor, 1)
      print('rem=',npts, C.getNPts(edge), factor)
      D._getCurvilinearAbscissa(edge)
      edgeno = getNo(edge)
      aedge = C.getFields([Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__], edge, api=2)[0]
      e = occ.meshOneEdge(hook, edgeno, -1, -1, -1, aedge)
      cad = Internal.getNodeFromName1(edge, 'CAD')
      render = Internal.getNodeFromName1(edge, '.RenderInfo')
      z = Internal.createZoneNode('edge%03d'%(edgeno), e, [],
                                  Internal.__GridCoordinates__,
                                  Internal.__FlowSolutionNodes__,
                                  Internal.__FlowSolutionCenters__)
      z[2].append(cad)
      if render is not None: z[2].append(render)
      b[2][edgeno-1] = z
  return None    

def getCADcontainer(t):
  hmax = None; hausd = None
  CAD = Internal.getNodeFromName1(t, 'CAD')
  if CAD is None: return [hmax, hausd]
  hmax = Internal.getNodeFromName1(CAD, 'hsize')
  if hmax is not None: hmax = Internal.getValue(hmax)
  hausd = Internal.getNodeFromName1(CAD, 'hausd')
  if hausd is not None: hausd = Internal.getValue(hausd)
  return [hmax, hausd]

# build or update edges:FaceList
def _updateEdgesFaceList__(t):

  # build edgeOfFaces
  edgeOfFaces = {}
  b = Internal.getNodeFromName1(t, 'EDGES')
  for e in Internal.getZones(b):
    edgeno = getNo(e)
    edgeOfFaces[edgeno] = []
    
  # update edgeOfFaces
  b = Internal.getNodeFromName1(t, 'FACES')
  for f in Internal.getZones(b):
    faceno = getNo(f)
    cad = Internal.getNodeFromName1(f, 'CAD')
    edgeList = Internal.getNodeFromName1(cad, 'edgeList')
    edgeList = edgeList[1]
    for i in edgeList: edgeOfFaces[i].append(faceno)

  # update EGDES:faceList
  b = Internal.getNodeFromName1(t, 'EDGES')
  for e in Internal.getZones(b):
    edgeno = getNo(e)
    cad = Internal.getNodeFromName1(e, 'CAD')
    faces = edgeOfFaces[edgeno]
    n = numpy.array(faces, dtype=Internal.E_NpyInt)
    Internal._createUniqueChild(cad, 'faceList', 'DataArray_t', value=n)
  return None

# mesh all faces or a subset from edges U
def _meshAllFacesTri(hook, t, metric=True, faceList=None, hList=[], hmax=-1, hausd=-1):

  b = Internal.getNodeFromName1(t, 'EDGES')
  dedges = []
  for z in Internal.getZones(b):
    pf = Internal.getNodeFromName2(z, 'u')
    if pf is None: print("Error: meshAllFacesTri: u field missing in edges.")
    e = C.getFields([Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__], z, api=2)[0]
    dedges.append(e)

  nbFaces = occ.getNbFaces(hook)
  if faceList is None:
    N = nbFaces // Cmpi.size
    nstart = Cmpi.rank*N
    nend = nstart+N
    if Cmpi.rank == Cmpi.size-1: nend = nbFaces
    faceList = range(nstart+1, nend+1)

  if hList is None or hList == []:
    if hausd < 0 and hmax > 0: hList = [(hmax,hmax,hausd)]*len(faceList)
    elif hausd > 0 and hmax < 0: hList = [(1.e-5,10000.,hausd)]*len(faceList)
    else: hList = [(hmax*0.8,hmax*1.2,hausd)]*len(faceList)

  faces = OCC.meshAllFacesTri(hook, dedges, metric, faceList, hList)

  b = Internal.getNodeFromName1(t, 'FACES')
  if b is None: b = Internal.newCGNSBase('FACES', parent=t)
  for c, f in enumerate(faces):
    if f is None: continue # Failed face
    noface = faceList[c]
    z = Internal.createZoneNode('face%03d'%(noface), f, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    edgeNo = OCC.occ.getEdgeNoByFace(hook, noface)
    # conserve hook, name, type
    r = Internal.createChild(z, "CAD", "UserDefinedData_t")
    Internal._createChild(r, "name", "DataArray_t", value="face%03d"%(noface))
    Internal._createChild(r, "type", "DataArray_t", value="face")
    Internal._createChild(r, "no", "DataArray_t", value=noface)
    Internal._createChild(r, "edgeList", "DataArray_t", value=edgeNo)
    Internal._createChild(r, "hsize", "DataArray_t", value=hList[c])
    #Internal._createChild(r, "hook", "UserDefinedData_t", value=hook)
    b[2].append(z)

  _updateEdgesFaceList__(t)
  _setLonelyEdgesColor(t)

  return None

# mesh all faces or a subset from edges U
def _meshAllFacesStruct(hook, t, faceList=None):

  b = Internal.getNodeFromName1(t, 'EDGES')
  dedges = []
  for z in Internal.getZones(b):
    pf = Internal.getNodeFromName2(z, 'u')
    if pf is None: print("Error: meshAllFaces: u field missing in edges.")
    e = C.getFields([Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__], z, api=2)[0]
    dedges.append(e)

  nbFaces = occ.getNbFaces(hook)
  if faceList is None:
    N = nbFaces // Cmpi.size
    nstart = Cmpi.rank*N
    nend = nstart+N
    if Cmpi.rank == Cmpi.size-1: nend = nbFaces
    faceList = range(nstart+1, nend+1)

  faces, noloct, nofacet = OCC.meshAllFacesStruct(hook, dedges, faceList)

  b = Internal.getNodeFromName1(t, 'FACES')
  if b is None: b = Internal.newCGNSBase('FACES', parent=t)
  for c, f in enumerate(faces):
    if f is None: continue # Failed face
    noface = nofacet[c]
    noloc = noloct[c]
    if noloc == 1: zname = 'face%03d'%(noface)
    else: zname = 'face%03d_%03d'%(noface, noloc)
    z = Internal.createZoneNode(zname, f, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    edgeNo = OCC.occ.getEdgeNoByFace(hook, noface)
    # conserve hook, name, type
    r = Internal.createChild(z, "CAD", "UserDefinedData_t")
    Internal._createChild(r, "name", "DataArray_t", value="face%03d"%(noface))
    Internal._createChild(r, "type", "DataArray_t", value="face")
    Internal._createChild(r, "no", "DataArray_t", value=noface)
    Internal._createChild(r, "edgeList", "DataArray_t", value=edgeNo)
    Internal._createChild(r, "hsize", "DataArray_t", value=-1)
    #Internal._createChild(r, "hook", "UserDefinedData_t", value=hook)
    b[2].append(z)

  _updateEdgesFaceList__(t)
  _setLonelyEdgesColor(t)

  return None

# set color red to lonelyEdges
def _setLonelyEdgesColor(t):
  import CPlot.PyTree as CPlot
  b = Internal.getNodeFromName1(t, 'EDGES')
  if b is None: return None
  zones = Internal.getZones(b)
  for ze in zones:
    CAD = Internal.getNodeFromName1(ze, 'CAD')
    faceList = Internal.getNodeFromName1(CAD, 'faceList')
    if faceList is not None:
      size = faceList[1].size
      if size == 2: # ok: blue
        CPlot._addRender2Zone(ze, color='Green')
      elif size == 1: # lonely: red
        CPlot._addRender2Zone(ze, color='Red')
      else: # strange!!
        CPlot._addRender2Zone(ze, color='Red')
  return None

# Return the number of lonely edges in t
def getNbLonelyEdges(t):
  b = Internal.getNodeFromName1(t, 'EDGES')
  if b is None: return 0
  nb = 0
  zones = Internal.getZones(b)
  for ze in zones:
    CAD = Internal.getNodeFromName1(ze, 'CAD')
    faceList = Internal.getNodeFromName1(CAD, 'faceList')
    if faceList is not None:
      size = faceList[1].size
      if size == 1: nb += 1
  return nb

#===========================================================================================
# Build interpData with ghostcells from a CAD PyTree t and its
# corresponding connectivity tree tc (in place)
def _setInterpData(t, tc):
  pos = getAllPos(t)
  nCADFaces = len(pos[1])
  print("pos[0]: {} ...".format(pos[0]))
  print("pos[1]: {} ...".format(pos[1]))
  print("nCADFaces: {} ...".format(nCADFaces))

  for faceNo in range(1, nCADFaces+1):
    print("faceNo: {} ...".format(faceNo))
    face = getFace(t, pos, faceNo)
    faceName = Internal.getName(face)
    # Get list of edges for that face, their range, and the total no of points
    edgeList = getEdgeListOfFace(t, pos, faceNo)
    r = getEdgeRangeOfFace(t, pos, faceNo, edgeList)
    print("  - edgeList: {}".format(edgeList))
    print("  - r: {}".format(r))
    print("  - nptsOnEdgesOfFace: {}".format(getNPtsOnEdgesOfFace(t, pos, faceNo, edgeList)))

    # Loop over edge indices of that face
    for edgeNo in edgeList:
      nptsFace = Internal.getValue(face)[0][0]
      # faceOpp via edge
      faceOppNo = getFaceNoOppOfEdge(t, pos, edgeNo, faceNo)
      if faceOppNo == -1: continue
      print("    + nptsFace: {}".format(nptsFace))
      print("    + edgeNo: {}".format(edgeNo))
      print("      * faceOppNo: {}".format(faceOppNo))
      faceOpp = getFace(t, pos, faceOppNo)
      faceOppName = Internal.getName(faceOpp)
      # Extract number of elements of faceOpp prior to adding any ghost cells
      rindOpp = Internal.getNodeFromType1(faceOpp, 'Rind_t')
      if rindOpp is None: neltsOpp = Internal.getValue(faceOpp)[0][1]
      else: neltsOpp = Internal.getValue(rindOpp)[0] - 1
      print("      * rindOpp: {}".format(rindOpp))
      print("      * neltsOpp: {}".format(neltsOpp))
      print("      * r: {}".format(getEdgeRangeOfFace(t, pos, faceNo, edgeList)[edgeNo]))
      print("      * rOpp: {}".format(getEdgeRangeOfFace(t, pos, faceOppNo)[edgeNo]))
      faceOpp = C.getAllFields(faceOpp, 'nodes')[0]

      # Extract vertex indices of edge
      vIdx = getEdgeVerticesOfFace(t, pos, faceNo).get(edgeNo)
      vIdxOpp = getEdgeVerticesOfFace(t, pos, faceOppNo).get(edgeNo)
      print("      * vIdx: {}".format(vIdx))
      print("      * vIdxOpp: {}".format(vIdxOpp))
      oppData, ptList, ptListDonor = occ.getOppData(faceOpp, vIdx, vIdxOpp, nptsFace, neltsOpp)
      oppData = C.convertArrays2ZoneNode('Zone', [oppData])
      # Append vertices and connectivity of opp. face to current face
      _addOppFaceData2Face(face, oppData, ptList, r[edgeNo])
      
      # Update connectivity tree using newly inserted indices in current face
      # (receiver) and their corresponding indices in opposite face (donor)
      ptList = ptList[-len(ptListDonor):]
      print("      * ptList: {}".format(ptList))
      print("      * ptListDonor: {}".format(ptListDonor))
      _updateConnectivityTree(tc, name=faceName, nameDonor=faceOppName,
                              ptList=ptList-1, ptListDonor=ptListDonor-1)
  return None

# Retourne l'edge a partir de edgeNo (numero global CAD)
def getEdge(t, pos, edgeNo):
  be = Internal.getNodeFromName1(t, 'EDGES')
  ze = be[2][pos[0][edgeNo]]
  return ze

# Retourne les edges a partir d'une liste de edgeNos (numero global CAD)
def getEdges(t, pos, edgeNos):
  be = Internal.getNodeFromName1(t, 'EDGES')
  locEdgeNos = [pos[0].get(k) for k in edgeNos]
  ze = [be[2][i] for i in locEdgeNos]
  return ze

# Retourne la face de faceNo (numero global CAD)
def getFace(t, pos, faceNo):
  bf = Internal.getNodeFromName1(t, 'FACES')
  zf = bf[2][pos[1][faceNo]]
  return zf

# Return the position of edgeNo in faceNo
def getEdgePosInFace(t, pos, faceNo, edgeNo):
  edgeList = getEdgeListOfFace(t, pos, faceNo)
  for c, e in enumerate(edgeList):
    if e == edgeNo: return c
  return -1 # not found

# Get edge list from faceNo
def getEdgeListOfFace(t, pos, faceNo):
  zf = getFace(t, pos, faceNo)
  CAD = Internal.getNodeFromName1(zf, 'CAD')
  edgeList = Internal.getNodeFromName1(CAD, 'edgeList')
  return edgeList[1]

# Get face list from edgeNo
def getFaceListOfEdge(t, pos, edgeNo):
  ze = getEdge(t, pos, edgeNo)
  CAD = Internal.getNodeFromName1(ze, 'CAD')
  faceList = Internal.getNodeFromName1(CAD, 'faceList')
  return faceList[1]

# Return the ranges of edges of faceNo
def getEdgeRangeOfFace(t, pos, faceNo, edgeList=None):
  if edgeList is None: edgeList = getEdgeListOfFace(t, pos, faceNo)
  ranges = {}
  c = 0
  for e in edgeList:
    ze = getEdge(t, pos, e)
    npts = C.getNPts(ze)
    ranges[e] = [c, c+npts]
    c += npts - 1 # last point is the 1st point of next edge
  return ranges

# Return the number of points on the edges of faceNo
def getNPtsOnEdgesOfFace(t, pos, faceNo, edgeList=None):
  if edgeList is None: edgeList = getEdgeListOfFace(t, pos, faceNo)
  npts = 0
  for e in edgeList:
    ze = getEdge(t, pos, e)
    npts += C.getNPts(ze) - 1 # last point is the 1st point of next edge
  return npts

# Return the vertex indices for all edges of faceNo
def getEdgeVerticesOfFace(t, pos, faceNo, edgeList=None):
  if edgeList is None: edgeList = getEdgeListOfFace(t, pos, faceNo)
  ranges = getEdgeRangeOfFace(t, pos, faceNo)
  values = [numpy.arange(*ranges.get(k), dtype=Internal.E_NpyInt) for k in edgeList]
  values[-1][-1] = values[0][0] # last index is a repetition of the first
  return dict((k, v) for k, v in zip(edgeList, values))

# Get the face opp of edgeNo belonging to faceNo
def getFaceNoOppOfEdge(t, pos, edgeNo, faceNo):
  faceList = getFaceListOfEdge(t, pos, edgeNo)
  for f in faceList:
    if f != faceNo: return f
  return -1

# Add opposite face vertices and connectivity to a face (in place)
def _addOppFaceData2Face(z, zOpp, ptList, rgEdge):
  # Get dimensions of current face
  nptsEdge = rgEdge[1] - rgEdge[0]
  dims = Internal.getValue(z)[0]
  
  # Add missing opposite vertices to current
  # NB: the first `nptsEdge` points in 'Opp' should be skipped
  coords = C.getFields('coords', z)
  coordsOpp = C.getFields('coords', zOpp)
  coordsv = Internal.getValue(coords[0])
  coordsOppv = Internal.getValue(coordsOpp[0])
  coords[0][1] = numpy.hstack((coordsv, coordsOppv[:,nptsEdge:]))
  z = C.setFields(coords, z, 'nodes')
  
  # Edit opposite connectivity and append to current
  n = Internal.getNodesFromType1(z, 'Elements_t')
  c = Internal.getNodeFromType1(n[0], 'DataArray_t')
  cn = Internal.getValue(c)
  nOpp = Internal.getNodesFromType1(zOpp, 'Elements_t')
  cOpp = Internal.getNodeFromType1(nOpp[0], 'DataArray_t')
  cnOpp = Internal.getValue(cOpp)
  c[1] = numpy.concatenate((cn, ptList[cnOpp-1]))
  print("      * cnOpp: {}".format(ptList[cnOpp-1]))
  
  # Edit current Rind using offsetted opposite element range
  elRgOpp = Internal.getNodeFromType1(nOpp[0], 'IndexRange_t')
  # offset elt count in opp using current
  elRgOppv = Internal.getValue(elRgOpp) + dims[1]
  ntotElts = elRgOppv[1] # shorthand
  rind = Internal.getNodeFromType1(z, 'Rind_t')
  if rind is None:
    Internal.newRind(value=elRgOppv, parent=z)
  else:
    rindv = Internal.getValue(rind)
    rindv[1] = ntotElts
    
  # Edit element range of merged connectivity
  elRg = Internal.getNodeFromType1(n[0], 'IndexRange_t')
  elRgv = Internal.getValue(elRg)
  elRgv[1] = ntotElts
  
  # Edit number of elements of current face
  dims = Internal.getValue(z)[0]
  dims[1] = ntotElts
  return None

# Update connectivity tree tc of donor face `nameDonor` using data from
# the opposite/receiver face `name` (in place)
def _updateConnectivityTree(tc, name, nameDonor, ptList, ptListDonor):
  bf = Internal.getNodeFromName1(tc, 'Base')
  z = Internal.getNodeFromName1(bf, nameDonor)
  if z is None:
    z = Internal.createNode(nameDonor, 'Zone_t', value=[0,0,0], parent=bf)
  Internal.createNode('ZoneType', 'ZoneType_t', value='Unstructured', parent=z)

  zsr = Internal.createNode('ID_'+name, 'ZoneSubRegion_t', value=name, parent=z)
  Internal.createNode('ZoneRole', 'DataArray_t', value='Donor', parent=zsr)
  Internal.createNode('GridLocation', 'GridLocation_t', value='Vertex', parent=zsr)
  
  npts = len(ptList)
  # indices des points receveur
  Internal.createNode('PointList', 'IndexArray_t', value=ptListDonor, parent=zsr)
  # indices des points donneur
  Internal.createNode('PointListDonor', 'IndexArray_t', value=ptList, parent=zsr)
  # coefficient a 1 (injection)
  data = numpy.ones((npts), dtype=numpy.float64)
  Internal.createNode('InterpolantsDonor', 'DataArray_t', value=data, parent=zsr)
  # type d'interpolation a 1 (injection)
  data = numpy.ones((npts), dtype=Internal.E_NpyInt)
  Internal.createNode('InterpolantsType', 'DataArray_t', value=data, parent=zsr)
  return None

#=============================================================================
# CAD fixing
#=============================================================================

# return ordered edgeList with possible negative number (meaning to be reversed)
# edges: list of arrays
# return edgeList: list of edge numbers in CAD
def orderEdgeList(edges, tol=1.e-10):
  """Order edges in a loop."""
  import Transform.PyTree as T
  import KCore.Vector as Vector
  out = [] # ordered list of edges
  pool = edges[:] # list copy
  cur = pool[0]; pool.pop(0)
  out.append(cur)
  P1p = C.getValue(cur, 'GridCoordinates', -1)
  while len(pool) > 0:
    found = False
    for c, p in enumerate(pool):
      P0 = C.getValue(p, 'GridCoordinates', 0)
      P1 = C.getValue(p, 'GridCoordinates', -1)
      if Vector.squareDist(P1p,P0) < tol*tol:
        cur = p; out.append(cur); P1p = P1; pool.pop(c); found=True; break
      if Vector.squareDist(P1p,P1) < tol*tol:
        cur = T.reorder(p,(-1,2,3))
        CAD = Internal.getNodeFromName1(p, 'CAD')
        no = Internal.getNodeFromName2(CAD, 'no')
        nov = Internal.getValue(no)
        Internal.setValue(no, -nov)
        cur.append(CAD)
        out.append(cur); P1p = P0; pool.pop(c); found=True; break
    if not found: break
  outno = []
  for e in out:
    CAD = Internal.getNodeFromName1(e, 'CAD')
    no = Internal.getNodeFromName1(CAD, 'no')
    no = Internal.getValue(no)
    outno.append(no)
  return outno

# sew a set of faces
# faces: face list numbers
def _sewing(hook, faces, tol=1.e-6):
  OCC.occ.sewing(hook, faces, tol)
  return None

# add fillet from edges with given radius
def _addFillet(hook, edges, radius, new2OldEdgeMap=[], new2OldFaceMap=[]):
  OCC.occ.addFillet(hook, edges, radius, new2OldEdgeMap, new2OldFaceMap)
  print(new2OldEdgeMap, new2OldFaceMap)
  return None

# edgeMap and faceMap are new2old maps
def _removeFaces(hook, faces, new2OldEdgeMap=[], new2OldFaceMap=[]):
  OCC.occ.removeFaces(hook, faces, new2OldEdgeMap, new2OldFaceMap)
  return None

# fill hole from edges
# edges: edge list numbers (must be ordered)
def _fillHole(hook, edges):
  OCC.occ.fillHole(hook, edges)
  return None

# Return the number of edges in CAD hook
def getNbEdges(hook):
  """Return the number of edges in CAD hook."""
  return OCC.occ.getNbEdges(hook)

# Return the number of faces in CAD hook
def getNbFaces(hook):
  """Return the number of faces in CAD hook."""
  return OCC.occ.getNbFaces(hook)

# Return the file and format used to load CAD in hook
def getFileAndFormat(hook):
  return OCC.occ.getFileAndFormat(hook)
  
# IN: new2old: new2old map
# IN: Nold: size of old entities 
# OUT: odl2new array
def getOld2NewMap(Nold, new2old):
  old2new = numpy.zeros( (Nold), dtype=Internal.E_NpyInt)
  old2new[:] = -1
  Nnew = len(new2old)
  for n in range(Nnew): 
    v = new2old[n]
    if v > 0: old2new[v-1] = n+1
  return old2new
 
# update numbering in persistent zones (common part)
def _updateCADNumbering__(b, old2new, name):
  zones = Internal.getZones(b)
  for z in zones:
    CAD = Internal.getNodeFromName1(z, 'CAD')
    if CAD is not None:
      n = Internal.getNodeFromName1(CAD, 'name')
      no = Internal.getValue(n)
      no = no.replace(name, '')
      no = int(no)
      no = old2new[no-1]
      Internal._setValue(n, name+'%03d'%no)
      z[0] = name+'%03d'%no
      n = Internal.getNodeFromName1(CAD, 'no')
      no = Internal.getValue(n)
      no = old2new[no-1]
      Internal._setValue(n, no)
  return None

# update numbering in persistent zones (edge part)
def _updateEdgesNumbering__(t, old2NewEdgeMap, old2NewFaceMap):
  b = Internal.getNodeFromName1(t, 'EDGES')
  _updateCADNumbering__(b, old2NewEdgeMap, 'edge')
  for z in Internal.getZones(b):
    CAD = Internal.getNodeFromName1(z, 'CAD')
    if CAD is not None:
      n = Internal.getNodeFromName1(CAD, 'faceList')
      if n is not None:
        pn = n[1]; index = []
        for i in range(pn.size): 
          pn[i] = old2NewFaceMap[pn[i]-1]
          if pn[i] == -1: index.append(i)
        n[1] = numpy.delete(pn, index)
  return None

# update numbering in persistent zones (faces part)
def _updateFacesNumbering__(t, old2NewEdgeMap, old2NewFaceMap):
  b = Internal.getNodeFromName1(t, 'FACES')
  _updateCADNumbering__(b, old2NewFaceMap, 'face')
  for z in Internal.getZones(b):
    CAD = Internal.getNodeFromName1(z, 'CAD')
    if CAD is not None:
      n = Internal.getNodeFromName1(CAD, 'edgeList')
      if n is not None:
        pn = n[1]; index = []
        for i in range(pn.size):
          pn[i] = old2NewEdgeMap[pn[i]-1]
          if pn[i] == -1: index.append(i)
        n[1] = numpy.delete(pn, index)
  return None

# update numbering (edges and faces)
def _updateNumbering(t, old2NewEdgeMap, old2NewFaceMap):
  _updateEdgesNumbering__(t, old2NewEdgeMap, old2NewFaceMap)
  _updateFacesNumbering__(t, old2NewEdgeMap, old2NewFaceMap)
  return None

# suppress edge or face depending on old2new
def _suppressZones__(b, old2new):
  rme = []
  for i in range(old2new.size):
    if old2new[i] == -1: rme.append(i+1)
  pos, posi = getPos(b)
  rme = list(reversed(rme))
  for c, f in enumerate(rme):
    cd = pos[f]
    del b[2][cd]
  return None

# add zones depending on new2old
def _addZones__(b, new2old):
  addme = []
  for i in range(new2old.size):
    if new2old[i] == -1: addme.append(i+1)
  # mesh and add zone...
  return None

# update tree when CAD hook is modified
def _updateTree(t, oldNbEdges, oldNbFaces, new2OldEdgeMap, new2OldFaceMap):
  """Update tree when CAD is mmodified."""
  old2NewEdgeMap = getOld2NewMap(oldNbEdges, new2OldEdgeMap)
  old2NewFaceMap = getOld2NewMap(oldNbFaces, new2OldFaceMap)
  
  # remove zones that are no longer in CAD
  b = Internal.getNodeFromName1(t, 'EDGES')
  _suppressZones__(b, old2NewEdgeMap)
  b = Internal.getNodeFromName1(t, 'FACES')
  _suppressZones__(b, old2NewFaceMap)

  # update numbering of persistent zones
  _updateNumbering(t, old2NewEdgeMap, old2NewFaceMap)

  # sort zone by names
  b = Internal.getNodeFromName1(t, 'EDGES')
  Internal._sortByName(b, recursive=False)
  b = Internal.getNodeFromName1(t, 'FACES')
  Internal._sortByName(b, recursive=False)

  # New edges and faces must be added manually after this function call

  # recolor
  _setLonelyEdgesColor(t)
  return None