"""OpenCascade definition module (pyTree).
"""
try:
    import OCC
    import Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
except ImportError: 
  raise ImportError("OCC.PyTree: requires Converter module.")

__version__ = OCC.__version__

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
def readCAD(fileName, format='fmt_iges'):
  """Read CAD and return a CAD hook."""
  return OCC.occ.readCAD(fileName, format)

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

# a mettre en convertCAD2PyTree?
def getTree(hook, N=11, hmax=-1, hausd=-1.):
  """Get a first TRI meshed tree linked to CAD."""
  
  t = C.newPyTree(['EDGES', 'SURFACES'])

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
    # Conserve hook, name, type et no de l'edge dans la CAD
    r = Internal.createChild(z, "CAD", "UserDefined_t")
    Internal._createChild(r, "name", "UserDefined_t", value="edge%03d"%(c+1))
    Internal._createChild(r, "type", "UserDefined_t", value="edge")
    Internal._createChild(r, "no", "UserDefined_t", value=(c+1))
    Internal._createChild(r, "hook", "UserDefined_t", value=hook)
    b[2].append(z)

  # Faces
  b = Internal.getNodeFromName1(t, 'SURFACES')
  faceNo = []
  m = OCC.meshTRI__(hook, N=N, hmax=hmax, hausd=hausd, faceNo=faceNo)
  #m = OCC.meshSTRUCT__(hook, N=N, faceNo=faceNo)

  for c, f in enumerate(m):
    noface = faceNo[c]
    z = Internal.createZoneNode(C.getZoneName('face%03d'%noface), f, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    edgeNo = OCC.occ.getEdgeNoByFace(hook, noface)
    # conserve hook, name, type
    r = Internal.createChild(z, "CAD", "UserDefined_t")
    Internal._createChild(r, "name", "UserDefined_t", value="face%03d"%noface)
    Internal._createChild(r, "type", "UserDefined_t", value="face")
    Internal._createChild(r, "no", "UserDefined_t", value=noface)
    Internal._createChild(r, "edgeList", "UserDefined_t", value=edgeNo)
    Internal._createChild(r, "hook", "UserDefined_t", value=hook)
    b[2].append(z)

  return t

# remesh tree from edges with new ue from EDGES
def remeshTreeFromEdges(hook, tp):

  t = C.newPyTree(['EDGES', 'SURFACES'])

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
    # Conserve hook, name, type et no de l'edge dans la CAD
    r = Internal.createChild(z, "CAD", "UserDefined_t")
    Internal._createChild(r, "name", "UserDefined_t", value="edge%03d"%(c+1))
    Internal._createChild(r, "type", "UserDefined_t", value="edge")
    Internal._createChild(r, "no", "UserDefined_t", value=(c+1))
    Internal._createChild(r, "hook", "UserDefined_t", value=hook)
    b[2].append(z)

  # Faces
  b = Internal.getNodeFromName1(t, 'SURFACES')
  faceNo = []
  m = OCC.meshTRIU__(hook, arrays, faceNo=faceNo)
  for c, f in enumerate(m):
    noface = faceNo[c]
    z = Internal.createZoneNode(C.getZoneName('face%03d'%noface), f, [],
                                Internal.__GridCoordinates__,
                                Internal.__FlowSolutionNodes__,
                                Internal.__FlowSolutionCenters__)
    edgeNo = OCC.occ.getEdgeNoByFace(hook, noface)
    # conserve hook, name, type
    r = Internal.createChild(z, "CAD", "UserDefined_t")
    Internal._createChild(r, "name", "UserDefined_t", value="face%03d"%noface)
    Internal._createChild(r, "type", "UserDefined_t", value="face")
    Internal._createChild(r, "no", "UserDefined_t", value=noface)
    Internal._createChild(r, "edgeList", "UserDefined_t", value=edgeNo)
    Internal._createChild(r, "hook", "UserDefined_t", value=hook)
    b[2].append(z)

  return t
