# -- Internal --
# -- PyTree node manipulations --
from sys import version_info
try: range = xrange
except: pass

import numpy
import fnmatch # unix wildcards
import KCore
from . import converter

# Containeurs
__GridCoordinates__ = 'GridCoordinates'
__FlowSolutionNodes__ = 'FlowSolution'
__FlowSolutionCenters__ = 'FlowSolution#Centers'

# Pour les BCs non structurees (nommage de transition courant/derniere norme)
__ELEMENTRANGE__ = 'ElementRange' # ElementRange/PointRange
__ELEMENTLIST__ = 'PointList' # ElementList/PointList
__FACELIST__ = 'PointList' # FaceList/PointList
__FACERANGE__ = 'PointRange' # FaceRange/PointRange
__POINTLIST__ = 'PointList' # Vertex list

# Separators
# Separateur intra-nom (pour le reperage des BCs qui deviennent des zones)
SEP1 = '/'
# Separateur intra-BCs (pour les exports getBCFaces et addBCFaces)
SEP2 = '@'

# other var names to CGNS var names
name2CGNS = { \
'x'              :'CoordinateX'                   , \
'y'              :'CoordinateY'                   , \
'z'              :'CoordinateZ'                   , \
'X'              :'CoordinateX'                   , \
'Y'              :'CoordinateY'                   , \
'Z'              :'CoordinateZ'                   , \
'ro'             :'Density'                       , \
'rou'            :'MomentumX'                     , \
'rov'            :'MomentumY'                     , \
'row'            :'MomentumZ'                     , \
'rovx'           :'MomentumX'                     , \
'rovy'           :'MomentumY'                     , \
'rovz'           :'MomentumZ'                     , \
'roe'            :'EnergyStagnationDensity'       , \
'roE'            :'EnergyStagnationDensity'       , \
'rok'            :'TurbulentEnergyKineticDensity' , \
'ronutilde'      :'TurbulentSANuTildeDensity'    , \
'roeps'          :'TurbulentDissipationDensity'   , \
'roomega'        :'TurbulentDissipationRateDensity', \
'mach'           :'Mach'                          , \
'psta'           :'Pressure'                      , \
'tsta'           :'Temperature'                   , \
'viscrapp'       :'Viscosity_EddyMolecularRatio'  , \
'walldistance'   :'TurbulentDistance'             , \
'wallglobalindex':'TurbulentDistanceIndex'        , \
}

# Known BCs
KNOWNBCS = ['BCWall', 'BCWallInviscid','BCWallViscous', 'BCWallViscousIsothermal',
            'BCFarfield', 'BCExtrapolate', 'BCInflow', 'BCInflowSubsonic', 'BCOutflow',
            'BCMatch', 'BCNearMatch', 'BCOverlap', 'BCSymmetryPlane', 
            'BCDegenerateLine', 'BCDegeneratePoint', 'BCStage',
            'UserDefined']

#import math
__DEG2RAD__ = 0.017453292519943295 #math.pi/180.
__RAD2DEG__ = 57.29577951308232    #180./math.pi

#==============================================================================
# -- is? --
#==============================================================================

# -- is node a top tree or something else?
def isTopTree(node):
    """Return True if node corresponds to a top tree node (CGSNTree_t)."""
    if len(node) != 4: return False
    if node[3] == 'CGNSTree_t': return True
    return False

# -- is node a standard node?
# Retourne -2 si node n'est pas un noeud standard
# Retourne -1 si node est un noeud standard de l'arbre
# Retourne 0 si node est une liste de noeuds standards (meme vide)
def isStdNode(node):
    """Return 0 if node is a list of standard pyTree nodes, -1 if node is a standard pyTree node, -2 otherwise."""
    if not isinstance(node, list): return -2
    if len(node) == 0: return 0
    node0 = node[0]
    if (isinstance(node0, str) and len(node) == 4 and
        isinstance(node[2], list)): return -1
    if not isinstance(node0, list): return -2
    if (isinstance(node0[0], str) and len(node0) == 4
        and isinstance(node0[2], list)): return 0
    return -2

# -- typeOfNode
# Retourne 1 si node est une zone
# Retourne 2 si node est une liste de zones
# Retourne 3 si node est un arbre
# Retourne 4 si node est une base
# Retourne 5 si node est une liste de base
# Retourne -1 sinon
def typeOfNode(node):
    """Return the type of node as an integer."""
    l = len(node)
    if l == 4:
        ntype = node[3]
        if ntype == 'Zone_t': return 1
        elif ntype == 'CGNSTree_t': return 3
        elif ntype == 'CGNSBase_t': return 4
        else:
            node1 = node[0]
            if len(node1) == 4:
                if node1[3] == 'Zone_t': return 2
                elif node1[3] == 'CGNSBase_t': return 5
            else: return -1
    elif l > 0:
        node1 = node[0]
        if len(node1) == 4:
            if node1[3] == 'Zone_t': return 2
            elif node1[3] == 'CGNSBase_t': return 5
        else: return -1
    else: return -1

# -- isType
# Compare given type and node type. Accepts widlcards.
# IN: ntype: string (Zone_t...)
# IN: node: pyTree node
def isType(node, ntype):
    """Return True if node is of given type."""
    tnode = node[3]
    if ('*' in ntype)|('?' in ntype)|('[' in ntype): return fnmatch.fnmatch(tnode, ntype)
    else: return tnode == ntype

# -- isName
# Compare given name and node name. Accepts widlcards.
# IN: name: string ou numpy char array
# IN: node: pyTree node
def isName(node, name):
    """Return True if node has given name."""
    if isinstance(name, numpy.ndarray): sname = name.tostring()
    else: sname = str(name)
    snode = node[0]
    if ('*' in sname)|('?' in sname)|('[' in sname): return fnmatch.fnmatch(snode, sname)
    else: return snode == sname

# -- isValue
# Compare given value and node value. Accepts wildcards.
def isValue(node, value):
    """Return True if node has given value."""
    # None case
    if node[1] is None: return (node[1] == value)
    # bool for string value
    isStrValue = isinstance(value, str)
    isNStrValue= (isinstance(value, numpy.ndarray) and value.dtype == numpy.array('a').dtype)
    # Not string case
    if not(isStrValue or isNStrValue):
        # Value to ndarray
        if not isinstance(value, numpy.ndarray): value = numpy.array([value]) # float or int
        # Node Value to ndarray
        if not isinstance(node[1], numpy.ndarray): nodeValue = numpy.array([node[1]]) # float or int or str not CGNS/Python compliant
        else: nodeValue = node[1]
        # Comparison on ndarray
        if value.dtype != nodeValue.dtype: return False
        elif value.shape != nodeValue.shape: return False 
        res = (value == nodeValue)
        return res.all()
    # String Case
    else:
        # bool for string node value
        isStrNode = isinstance(node[1], str)
        isNStrNode=(isinstance(node[1], numpy.ndarray) and node[1].dtype == numpy.array('a').dtype)
        # node value not string
        if not(isStrNode or isNStrNode): return False
        # Value  to string if ndarray string
        if isNStrValue: value = value.tostring()
        # Node value to string if ndarray string
        if isNStrNode: nodeValue = node[1].tostring()
        else: nodeValue = node[1] # not CGNS/Compliant
        # Comparison
        if ('*' in value)|('?' in value)|('[' in value): res = fnmatch.fnmatch(nodeValue, value)
        else: res = (nodeValue == value)
        return res

# -- isNameAndType
# Compare type and name of node with given (name, ntype)
def isNameAndType(node, name, ntype):
    return isName(node, name) and isType(node, ntype)

# -- is child
# Return true if idNode is a child of start (even at deep levels)
def isChild(start, node):
    """Return True if node is a child of node start."""
    if id(start) == id(node): return True
    else:
       for n in start[2]: 
        ret = isChild(n, node)
        if ret: return True
    return False

# Only first level
def isChild1(start, node):
    if id(start) == id(node): return True
    for n in start[2]:
        if id(n) == id(node): return True
    return False

# Only second level
def isChild2(start, node):
    if id(start) == id(node): return True
    for n in start[2]:
        if id(n) == id(node): return True
        for p in n[2]:
            if id(p) == id(node): return True
    return False
    
#==============================================================================
# -- Set/create nodes --
#==============================================================================

# -- setValue
# set value in node (modify node)
def setValue(node, value=None):
    """Set given value in node."""
    if node[3] == 'CGNSLibraryVersion_t':
        if isinstance(value, int) or isinstance(value, float): node[1] = numpy.array([value], 'f')
        elif isinstance(value, numpy.ndarray): node[1] = numpy.array(value, 'f')
        else: raise TypeError("setValue: CGNSLibraryVersion node value should be a float.")
    else:
        if value is None: node[1] = None
        elif isinstance(value, numpy.ndarray): node[1] = value
        elif isinstance(value, int) or isinstance(value, numpy.int32) or isinstance(value,numpy.int64): node[1] = numpy.array([value],'i')
        elif isinstance(value, float) or isinstance(value, numpy.float32) or isinstance(value, numpy.float64): node[1] = numpy.array([value],'d')
        elif isinstance(value, str): node[1] = numpy.array([c for c in value],'c')
        elif isinstance(value, list):
            testValue = value
            while isinstance(testValue, list): testValue = testValue[0]
            if isinstance(testValue, float) or isinstance(testValue, numpy.float32) or isinstance(testValue, numpy.float64):
                node[1] = numpy.array(value, dtype=numpy.float64, order='F')
            elif isinstance(testValue, int) or isinstance(testValue, numpy.int32) or isinstance(testValue, numpy.int64):
                node[1] = numpy.array(value, dtype=numpy.int32, order='F')
            elif isinstance(testValue, str):
                if isinstance(value[0], str):
                    v = numpy.empty( (32,len(value) ), dtype='c', order='F')
                    for c, i in enumerate(value): 
                        s = min(len(i),32)
                        v[:,c] = ' '
                        v[0:s,c] = i[0:s]
                    node[1] = v
                else:
                    v = numpy.empty( (32,len(value[0]),len(value) ), dtype='c', order='F')
                    for c in range(len(value)):
                        for d in range(len(value[c])):
                            s = min(len(value[c][d]),32)
                            v[:,d,c] = ' '
                            v[0:s,d,c] = value[c][d][0:s]
                    node[1] = v
        elif isinstance(value, tuple):
            testValue = value
            while isinstance(testValue, tuple): testValue = testValue[0]
            if isinstance(testValue, float) or isinstance(testValue, numpy.float32) or isinstance(testValue, numpy.float64):
                node[1] = numpy.array(value, dtype=numpy.float64, order='F')
            elif isinstance(testValue, int) or isinstance(testValue, numpy.int32) or isinstance(testValue, numpy.int64):
                node[1] = numpy.array(value, dtype=numpy.int32, order='F')
            elif isinstance(testValue, str):
                if isinstance(value[0], str):
                    v = numpy.empty( (32,len(value) ), dtype='c', order='F')
                    for c, i in enumerate(value): v[:,c] = i[:]
                    node[1] = v
                else:
                    v = numpy.empty( (32,len(value[0]),len(value) ), dtype='c', order='F')
                    for c in range(len(value)):
                        for d in range(len(value[c])): v[:,d, c] = value[c][d][:]
                    node[1] = v
        else: node[1] = numpy.array([value])
    return None

def _setValue(node, value=None):
    setValue(node, value)
    return None

# -- setName
# set node name (modify node)
def setName(node, name):
    """Set name in node."""
    if isinstance(name, str): node[0] = name
    elif isinstance(name, numpy.ndarray):
        if name.dtype.char == 'S': node[0] = name.tostring()
        elif name.dtype.char == 'c': node[0] = name.tostring()
        else:
            raise TypeError("setName: name of node must be a string(%s)"%(name.__repr__()[:min(len(name.__repr__()),60)]))
    else: raise TypeError("setName: name of node must be a string(%s)"%(name.__repr__()[:min(len(name.__repr__()),60)]))
    return None

def _setName(node, name):
    setName(node, name)
    return None

# -- setType
# Set node type (modify node)
def setType(node, ntype):
    """Set type in node."""
    if isinstance(ntype, str): node[3] = ntype
    else: raise TypeError("setType: type of a node must be a string (%s)"%(ntype.__repr__()[:min(len(ntype.__repr__()),60)]))
    return None

def _setType(node, ntype):
    setType(node, ntype)
    return None

# -- createNode
# create a node
def createNode(name, ntype, value=None, children=None, parent=None):
    """Create a pyTree node."""
    if children is None: node = [None, None, [], None]
    else: node = [None, None, children, None]
    setName(node, name)
    setType(node, ntype)
    setValue(node, value)
    if parent is not None: parent[2].append(node)
    return node

# -- addChild (modify node)
def addChild(node, child, pos=-1):
    """Add a child node to node's children."""
    # On devrait ajouter le n=copyRef(node)
    _addChild(node, child, pos=pos)
    return child

def _addChild(node, child, pos=-1):
    isStd = isStdNode(child)
    node2 = node[2]
    if isStd == -1:
      if pos == -1: node2.append(child)
      elif pos == 0: node2.insert(0, child)
      else: node2.insert(pos, child)
    elif isStd == 0:
      if pos == -1: node2.extend(child)
      elif pos == 0: node2[:] = child+node2
      else: node2[:] = node2[:pos]+child+node2[pos:]
    return None

# -- createChild (modify node)
def createChild(node, name, ntype, value=None, children=None, pos=-1):
    """Create a node and attach it to node's children."""
    return addChild(node, createNode(name, ntype, value, children), pos)

def _createChild(node, name, ntype, value=None, children=None, pos=-1):
    _addChild(node, createNode(name, ntype, value, children), pos)
    return None

# -- createUniqueChild (modify node)
def createUniqueChild(node, name, ntype, value=None, children=None, pos=-1):
    """Create a unique child in node's children."""
    # name exists?
    e = -1; i = 0
    for n in node[2]:
        if name == n[0]: e = i; break
        i += 1
    if e == -1: child = createChild(node, name, ntype, value=value,
                                      children=children, pos=pos)
    else: # replace value, ntype, children
        child = node[2][e]
        setValue(child, value)
        setType(child, ntype)
        if children is not None: child[2] = children
    return child

def _createUniqueChild(node, name, ntype, value=None, children=None, pos=-1):
    createUniqueChild(node, name, ntype, value, children, pos)
    return None

# -- append
def append(t, node, path):
    """Append a node to t specifying its path in t."""
    tp = copyRef(t)
    _append(tp, node, path)
    return tp

def _append(t, node, path):
    root = getNodeFromPath(t, path)
    # already exists?
    child = getNodeFromPath(root, '%s'%node[0])
    if child is None: _addChild(root, node)
    elif node[1] is not None:
        child[1] = node[1]; child[3] = node[3]
        for child in node[2]: # append les enfants de node uniquement
            _append(root, child, '%s'%node[0])
    else:
        for child in node[2]: # append les enfants de node uniquement
            _append(root, child, '%s'%node[0])
    return None

# -- create empty nodes from a path
def createNodesFromPath(t, path):
    tp = copyRef(t)
    _createNodesFromPath(t, path)
    return tp

def _createNodesFromPath(t, path):
    paths = path.split('/')
    p = t
    for i in paths:
        if i != '':
            r = getNodeFromName1(p, i)
            if r is None:
                p = createChild(p, i, 'DataArray_t', None)
            else: p = r
    return None

#==============================================================================
# -- specific CGNS nodes creation --
#==============================================================================

# -- Create a CGNS version node
def createCGNSVersionNode():
    version = numpy.array([3.1], dtype=numpy.float64)
    return ['CGNSLibraryVersion', version, [], 'CGNSLibraryVersion_t']

# -- Cree le noeud root de l'arbre
def createRootNode(name='CGNSTree', children=None):
    if children is None: children = []
    return [name, None, children, 'CGNSTree_t']

# -- Create base node named name with dim
# cellDim=2 (cells surfaciques), cellDim=3 (cells volumiques)
def createBaseNode(name, cellDim):
    a = numpy.empty((2), numpy.int32)
    a[0] = cellDim; a[1] = 3
    return [name, a, [], 'CGNSBase_t']

# -- newCGNSTree
# create a pyTree
def newCGNSTree():
    """Create a new pyTree."""
    t = createRootNode()
    t[2].append(createCGNSVersionNode())
    return t

# -- newCGNSBase
def newCGNSBase(name='Base', cellDim=3, physDim=3, parent=None):
    """Create a new Base node."""
    if parent is None:
        node = createNode(name, 'CGNSBase_t', value=[cellDim,physDim])
    else:
        node = createUniqueChild(parent, name, 'CGNSBase_t',
                                 value=[cellDim,physDim])
    return node

# -- newZone
def newZone(name='Zone', zsize=None, ztype='Structured',
            family=None, parent=None):
    """Create a new Zone node."""
    ZoneType_l = ['Null', 'UserDefined', 'Structured', 'Unstructured']
    if ztype not in ZoneType_l:
        raise ValueError('newZone: ztype must be in %s.'%str(ZoneType_l))
    if parent is None:
        node = createNode(name, 'Zone_t', value=zsize)
    else:
        node = createUniqueChild(parent, name, 'Zone_t', value=zsize)
    createUniqueChild(node, 'ZoneType', 'ZoneType_t', value=ztype)
    if family is not None:
        createUniqueChild(node, 'FamilyName', 'FamilyName_t', value=family)
    return node

# -- newGridCoordinates
def newGridCoordinates(name=__GridCoordinates__, parent=None):
    """Create a GridCoordinates node."""
    if parent is None:
        node = createNode(name, 'GridCoordinates_t')
    else: node = createUniqueChild(parent, name, 'GridCoordinates_t')
    return node

# -- newDataArray
def newDataArray(name='Data', value=None, parent=None):
    """Create a new DataArray node."""
    if parent is None:
        node = createNode(name, 'DataArray_t', value=value)
    else: node = createUniqueChild(parent, name, 'DataArray_t', value=value)
    return node

# -- newDataClass
def newDataClass(value='UserDefined', parent=None):
    """Create a new DataClass node."""
    if parent is None:
        node = createNode('DataClass', 'DataClass_t', value=value)
    else: node = createUniqueChild(parent, 'DataClass',
                                   'DataClass_t', value=value)
    return node

def newDimensionalUnits(massUnit='Kilogram', lengthUnit='Meter', timeUnit='Second', temperatureUnit='Kelvin',
                        angleUnit='Radian', parent=None):
    """Create a new DimensionalUnits node."""
    #unit in ['Null', 'UserDefined', 'Kilogram', 'Gram', 'Slug', 'PoundMass',
    #'Meter', 'Centimeter', 'Millimeter', 'Foot', 'Inch', 'Second',
    #'Kelvin', 'Celsius', 'Rankine', 'Radian', 'Degree']
    if parent is None:
        node = createNode('DimensionalUnits', 'DimensionalUnits_t', value=[massUnit, lengthUnit, 
                          timeUnit, temperatureUnit, angleUnit])
    else:
        node = createUniqueChild(parent, 'DimensionalUnits', 'DimensionalUnits_t', value=[massUnit, lengthUnit, 
                          timeUnit, temperatureUnit, angleUnit])
    return node

# -- newDimensionalExponents
def newDimensionalExponents(massExponent=0., lengthExponent=0.,
                            timeExponent=0., temperatureExponent=0.,
                            angleExponent=0., parent=None):
    """Create a new DimensionalExponents node."""
    value = [massExponent, lengthExponent, timeExponent,
             temperatureExponent, angleExponent]
    if parent is None:
        node = createNode('DimensionalExponents', 'DimensionalExponents_t',
                          value=value)
    else:
        node = createUniqueChild(parent, 'DimensionalExponents',
                                 'DimensionalExponents_t', value=value)
    return node

# -- newDataConversion
def newDataConversion(conversionScale=1., conversionOffset=0., parent=None):
    """Create a new DataConversion node."""
    if parent is None:
        node = createNode('DataConversion', 'DataConversion_t',
                          value=[conversionScale, conversionOffset])
    else: node = createUniqueChild(parent, 'DataConversion', 'DataConversion_t',
                                   value=[conversionScale, conversionOffset])
    return node

# -- newDescriptor
def newDescriptor(name='Descriptor', value='', parent=None):
    """Create a new Descriptor node."""
    if parent is None:
        node = createNode(name, 'Descriptor_t', value=value)
    else: node = createUniqueChild(parent, name, 'Descriptor_t', value=value)
    return node

# -- newGridLocation
def newGridLocation(value='CellCenter', parent=None):
    """Create a new GridLocation node."""
    GridLocation_l = ['Null', 'UserDefined', 'Vertex', 'CellCenter',
                      'FaceCenter', 'IFaceCenter', 'JFaceCenter',
                      'KFaceCenter', 'EdgeCenter']
    if value not in GridLocation_l:
        raise ValueError('newGridLocation: value must be in %s.'%str(GridLocation_l))
    if parent is None:
        node = createNode('GridLocation', 'GridLocation_t', value=value)
    else: node = createUniqueChild(parent, 'GridLocation', 'GridLocation_t',
                                   value=value)
    return node

# -- newIndexArray
def newIndexArray(name='Index', value=None, parent=None):
    """Create a new IndexArray node."""
    if parent is None:
        node = createNode(name, 'IndexArray_t', value=value)
    else: node = createUniqueChild(parent, name, 'IndexArray_t',
                                   value=value)
    return node

# -- newPointList
def newPointList(name='PointList', value=None, parent=None):
    """Create a new PointList node."""
    if parent is None:
        node = createNode(name, 'IndexArray_t', value=value)
    else: node = createUniqueChild(parent, name, 'IndexArray_t',
                                   value=value)
    return node

# -- newPointRange
def newPointRange(name='PointRange', value=None, parent=None):
    """Create a new PointRange node."""
    if parent is None:
        node = createNode(name, 'IndexRange_t', value=value)
    else: node = createUniqueChild(parent, name, 'IndexRange_t',
                                   value=value)
    return node

# -- newRind
def newRind(value=None, parent=None):
    """Create a new Rind node."""
    if parent is None:
        node = createNode('Rind', 'Rind_t', value=value)
    else: node = createUniqueChild(parent, 'Rind', 'Rind_t', value=value)
    return node

# -- newSimulation
def newSimulationType(value='TimeAccurate', parent=None):
    """Create a new SimulationType node."""
    SimulationType_l = ['Null', 'UserDefined', 'TimeAccurate', 'NonTimeAccurate']
    if value not in SimulationType_l:
        raise ValueError('newSimulationType: value must be in %s.'%str(SimulationType_l))
    if parent is None:
        node = createNode('SimulationType', 'SimulationType_t', value=value)
    else: node = createUniqueChild(parent, 'SimulationType',
                                   'SimulationType_t', value=value)
    return node

# -- newOrdinal
def newOrdinal(value=0, parent=None):
    """Create a new Ordinal node."""
    if parent is None:
        node = createNode('Ordinal', 'Ordinal_t', value=value)
    else: node = createUniqueChild(parent, 'Ordinal',
                                   'Ordinal_t', value=value)
    return node

# -- newDiscreteData
def newDiscreteData(name='DiscreteData', parent=None):
    """Create a new DiscreteData node."""
    if parent is None:
        node = createNode(name, 'DiscreteData_t')
    else: node = createUniqueChild(parent, name, 'DiscreteData_t')
    return node

# -- newIntegralData
def newIntegralData(name='Integral', parent=None):
    """Create a new IntegralData node."""
    if parent is None:
        node = createNode(name, 'IntegralData_t')
    else: node = createUniqueChild(parent, name, 'IntegralData_t')
    return node

# -- newElements
def newElements(name='Elements', etype='UserDefined',
                econnectivity=None,
                erange=None, eboundary=0, parent=None):
    """Create a new Elements node."""
    if isinstance(etype, int): etp = etype
    else: etp, nnodes = eltName2EltNo(etype)
    if parent is None:
        node = createNode(name, 'Elements_t', value=[etp,eboundary])
    else: node = createUniqueChild(parent, name, 'Elements_t',
                                   value=[etp,eboundary])
    newDataArray('ElementConnectivity', econnectivity, parent=node)
    #if erange is None: erange = numpy.ones(2, dtype=numpy.int32)
    newPointRange('ElementRange', erange, parent=node)
    return node

# -- newZoneBC
def newZoneBC(parent=None):
    """Create a new ZoneBC node."""
    if parent is None:
        node = createNode('ZoneBC', 'ZoneBC_t')
    else: node = createUniqueChild(parent, 'ZoneBC', 'ZoneBC_t')
    return node

# -- newBoundary
def newBC(name='BC', pointRange=None, pointList=None,
          btype='Null', family=None, parent=None):
    """Create a new BC node."""
    BCType_l = ['Null', 'UserDefined', 'BCGeneral', 'BCDirichlet',
                'BCNeumann', 'BCExtrapolate', 'BCWallInviscid',
                'BCWallViscousHeatFlux', 'BCWallViscousIsothermal',
                'BCWallViscous', 'BCWall',
                'BCInflowSubsonic', 'BCInflowSupersonic',
                'BCOutflowSubsonic', 'BCOutflowSupersonic',
                'BCTunnelInflow', 'BCTunnelOutflow',
                'BCDegenerateLine', 'BCDegeneratePoint',
                'BCSymmetryPlane', 'BCSymmetryPolar',
                'BCAxisymmetricWedge', 'FamilySpecified',
                'BCFarfield', 'BCInflow', 'BCOutflow']
    if btype not in BCType_l:
        raise ValueError('newBC: btype must be in %s.'%str(BCType_l))
    if parent is None:
        node = createNode(name, 'BC_t', value=btype)
    else: node = createUniqueChild(parent, name, 'BC_t', value=btype)
    if pointRange is not None: newPointRange(value=pointRange, parent=node)
    if pointList is not None: newPointList(value=pointList, parent=node)
    if family is not None:
        createChild(node, 'FamilyName', 'FamilyName_t', value=family)
    return node

# -- newBCDataSet
def newBCDataSet(name='BCDataSet', value='Null', gridLocation=None, parent=None):
    """Create BCDataSet node."""
    BCType_l = ['Null', 'UserDefined', 'BCGeneral', 'BCDirichlet',
                'BCNeumann', 'BCExtrapolate', 'BCWallInviscid',
                'BCWallViscousHeatFlux', 'BCWallViscousIsothermal',
                'BCWallViscous', 'BCWall',
                'BCInflowSubsonic', 'BCInflowSupersonic',
                'BCOutflowSubsonic', 'BCOutflowSupersonic',
                'BCTunnelInflow', 'BCTunnelOutflow',
                'BCDegenerateLine', 'BCDegeneratePoint',
                'BCSymmetryPlane', 'BCSymmetryPolar',
                'BCAxisymmetricWedge', 'FamilySpecified',
                'BCFarfield', 'BCInflow', 'BCOutflow']
    if value not in BCType_l:
        raise ValueError('newBCDataSet: value must be in %s.'%str(BCType_l))
    if parent is None:
        node = createNode(name, 'BCDataSet_t', value=value)
    else: node = createUniqueChild(parent, name, 'BCDataSet_t', value=value)
    if gridLocation is not None: newGridLocation(gridLocation, parent=node)
    return node

# -- newBCData
def newBCData(name='BCData', parent=None):
    """Create a new BCData node."""
    if parent is None:
        node = createNode(name, 'BCData_t')
    else: node = createUniqueChild(parent, name, 'BCData_t')
    return node

# -- newBCProperty
def newBCProperty(wallFunction='Null', area='Null', parent=None):
    """Create BCProperty node."""
    if parent is None:
        node = createNode('BCProperty', 'BCProperty_t')
    else: node = createUniqueChild(parent, 'BCProperty', 'BCProperty_t')
    createUniqueChild(node, 'WallFunctionType', 'WallFunctionType_t',
                      value=wallFunction)
    createUniqueChild(node, 'Area', 'Area_t', value=area)
    return node

# -- newAxiSymmetry
def newAxiSymmetry(referencePoint=[0.,0.,0.],
                   axisVector=[0.,0.,0.], parent=None):
    """Create AxiSymetry node."""
    if parent is None:
        node = createNode('AxiSymmetry', 'AxiSymmetry_t')
    else: node = createUniqueChild(parent, 'AxiSymmetry', 'AxiSymmetry_t')
    newDataArray('AxiSymmetryReferencePoint', value=referencePoint, parent=node)
    newDataArray('AxiSymmetryAxisVector', value=axisVector, parent=node)
    return node

# -- newRotatingCoordinates
def newRotatingCoordinates(rotationCenter=[0.,0.,0.],
                           rotationRateVector=[0.,0.,0.], parent=None):
    """Create a new RotatingCoordinates node."""
    if parent is None:
        node = createNode('RotatingCoordinates', 'RotatingCoordinates_t')
    else: node = createUniqueChild(parent, 'RotatingCoordinates',
                                   'RotatingCoordinates_t')
    newDataArray('RotationCenter', value=rotationCenter, parent=node)
    newDataArray('RotationRateVector', value=rotationRateVector, parent=node)
    return node

# -- newFlowSolution
def newFlowSolution(name=__FlowSolutionNodes__,
                    gridLocation='Vertex', parent=None):
    """Create a new FlowSolution node."""
    if parent is None:
        node = createNode(name, 'FlowSolution_t')
    else: node = createUniqueChild(parent, name, 'FlowSolution_t')
    if gridLocation is not None: newGridLocation(gridLocation, node)
    return node

# -- newZoneGridConnectivity
def newZoneGridConnectivity(name='ZoneGridConnectivity', parent=None):
    """Create a new ZoneGridConnectivity node."""
    if parent is None:
        node = createNode(name, 'ZoneGridConnectivity_t')
    else: node = createUniqueChild(parent, name, 'ZoneGridConnectivity_t')
    return node

# -- newGridConnectivity1to1
def newGridConnectivity1to1(name='Match', donorName=None,
                            pointRange=None, pointList=None,
                            pointRangeDonor=None, pointListDonor=None,
                            transform=None, parent=None):
    """Create a new GridConnectivity1to1 node."""
    if parent is None:
        node = createNode(name, 'GridConnectivity1to1_t', value=donorName)
    else: node = createUniqueChild(parent, name, 'GridConnectivity1to1_t',
                                   value=donorName)
    if pointRange is not None: newPointRange(value=pointRange, parent=node)
    if pointList is not None: newPointList(value=pointList, parent=node)
    if pointRangeDonor is not None:
        newPointRange(name='PointRangeDonor', value=pointRangeDonor, parent=node)
    if pointListDonor is not None:
        newPointList(name='PointListDonor', value=pointListDonor, parent=node)
    return node

# -- newGridConnectivity
def newGridConnectivity(name='Overlap', donorName=None,
                        ctype='Overset', parent=None):
    """Create a new GridConnectivity node."""
    if parent is None:
        node = createNode(name, 'GridConnectivity_t', value=donorName)
    else: node = createUniqueChild(parent, name,
                                   'GridConnectivity_t', value=donorName)
    newGridConnectivityType(ctype, node)
    return node

# -- newGridConnectivityType
def newGridConnectivityType(ctype='Overset', parent=None):
    """Create a new GridConnectivityType node."""
    if parent is None:
        node = createNode('GridConnectivityType', 'GridConnectivityType_t', value=ctype)
    else: node = createUniqueChild(parent, 'GridConnectivityType',
                                   'GridConnectivityType_t', value=ctype)
    return node

# -- newGridConnectivityProperty
def newGridConnectivityProperty(parent=None):
    """Create a new GridConnectivityType node."""
    if parent is None:
        node = createNode('GridConnectivityProperty',
                          'GridConnectivityProperty_t')
    else: node = createUniqueChild(parent, 'GridConnectivityProperty',
                                   'GridConnectivityProperty_t')
    return node

# -- newPeriodic
def newPeriodic(rotationCenter=[0.,0.,0.],
                rotationAngle=[0.,0.,0.], translation=[0.,0.,0.], parent=None):
    """Create a new Periodic node."""
    if parent is None:
        node = createNode('Periodic', 'Periodic_t')
    else: node = createUniqueChild(parent, 'Periodic', 'Periodic_t')
    newDataArray('RotationCenter', value=rotationCenter, parent=node)
    newDataArray('RotationAngle', value=rotationAngle, parent=node)
    newDataArray('Translation', value=translation, parent=node)
    return node

# -- newAverageInterface
def newAverageInterface(value='Null', parent=None):
    """Create a new AverageInterface node."""
    #if value not in AverageInterface_l:
    #    raise ValueError('newAverageInterface: value must be in %s.'%str(AverageInterface_l))
    if parent is None:
        node = createNode('AverageInterface', 'AverageInterface_t', value=value)
    else: node = createUniqueChild(parent, 'AverageInterface',
                                   'AverageInterface_t', value=value)
    return node

# -- newOversetHoles
def newOversetHoles(name='OversetHoles', pointRange=None,
                    pointList=None, parent=None):
    """Create a new OversetHoles node."""
    if parent is None:
        node = createNode(name, 'OversetHoles_t')
    else: node = createUniqueChild(parent, name, 'OversetHoles_t')
    if pointRange is not None: newPointRange(value=pointRange, parent=node)
    if pointList is not None: newPointList(value=pointList, parent=node)
    return node

# -- newFlowEquationSet
def newFlowEquationSet(parent=None):
    """Create a new FlowEquationSet node."""
    if parent is None:
        node = createNode('FlowEquationSet', 'FlowEquationSet_t')
    else: node = createUniqueChild(parent, 'FlowEquationSet',
                                   'FlowEquationSet_t')
    return node

# -- newGoverningEquations
def newGoverningEquations(value='Euler', parent=None):
    """Create a new GoverningEquation node."""
    GoverningEquations_l = ['Null', 'UserDefined', 'FullPotential',
                            'Euler', 'NSLaminar','NSTurbulent',
                            'NSLaminarIncompressible',
                            'NSTurbulentIncompressible']
    if value not in GoverningEquations_l:
        raise ValueError('newGoverningEquation: value must be in %s.'%str(GoverningEquations_l))
    if parent is None:
        node = createNode('GoverningEquations', 'GoverningEquations_t',
                          value=value)
    else: node = createUniqueChild(parent, 'GoverningEquations',
                                   'GoverningEquations_t', value=value)
    return node

# -- newGasModel
def newGasModel(value='Ideal', parent=None):
    """Create a new GasModel node."""
    GasModel_l = ['Null', 'UserDefined', 'Ideal', 'VanderWaals',
                  'CaloricallyPerfect', 'ThermallyPerfect',
                  'ConstantDensity', 'RedlichKwong']
    if value not in GasModel_l:
        raise ValueError('newGasModel: value must be in %s.'%str(GasModel_l))
    if parent is None:
        node = createNode('GasModel', 'GasModel_t', value=value)
    else: node = createUniqueChild(parent, 'GasModel',
                                   'GasModel_t', value=value)
    return node

# -- newThermalConductivityModel
def newThermalConductivityModel(value='SutherlandLaw', parent=None):
    """Create a new ThermalConductivity node."""
    ThermalConductivityModel_l = ['Null', 'UserDefined', 'ConstantPrandtl',
                                  'PowerLaw', 'SutherlandLaw']
    if value not in ThermalConductivityModel_l:
        raise ValueError('newThermalConductivityModel: value must be in %s'%str(ThermalConductivityModel_l))
    if parent is None:
        node = createNode('ThermalConductivityModel',
                          'ThermalConductivityModel_t', value=value)
    else: node = createUniqueChild(parent, 'ThermalConductivityModel',
                                   'ThermalConductivityModel_t', value=value)
    return node

# -- newViscosityModel
def newViscosityModel(value='SutherlandLaw', parent=None):
    """Create a new ViscosityModel node."""
    ViscosityModel_l = ['Null', 'UserDefined', 'Constant', 'PowerLaw',
                        'SutherlandLaw']
    if value not in ViscosityModel_l:
        raise ValueError('newViscosityModel: value must be in %s.'%str(ViscosityModel_l))
    if parent is None:
        node = createNode('ViscosityModel', 'ViscosityModel_t',
                          value=value)
    else: node = createUniqueChild(parent, 'ViscosityModel',
                                   'ViscosityModel_t', value=value)
    return node

# -- newTurbulenceClosure
def newTurbulenceClosure(value='Null', parent=None):
    """Create a new TurbulenceClosure node."""
    TurbulenceClosure_l = ['Null', 'UserDefined', 'EddyViscosity',
                           'ReynoldsStress', 'ReynoldsStressAlgebraic']
    if value not in TurbulenceClosure_l:
        raise ValueError('newTurbulenceColsure: value must be in %s.'%str(TurbulenceClosure_l))
    if parent is None:
        node = createNode('TurbulenceClosure', 'TurbulenceClosure_t',
                          value=value)
    else: node = createUniqueChild(parent, 'TurbulenceClosure',
                                   'TurbulenceClosure_t', value=value)
    return node

# -- newTurbulenceModel
def newTurbulenceModel(value='OneEquation_SpalartAllmaras', parent=None):
    """Create a new TurbulenceModel node."""
    TurbulenceModel_l = ['Null', 'UserDefined', 'Algebraic_BaldwinLomax',
                         'Algebraic_CebeciSmith', 'HalfEquation_JohnsonKing',
                         'OneEquation_BaldwinBarth',
                         'OneEquation_SpalartAllmaras',
                         'TwoEquation_JonesLaunder',
                         'TwoEquation_MenterSST', 'TwoEquation_Wilcox']
    if value not in TurbulenceModel_l:
        raise ValueError('newTurbulenceModle: value must be in %s.'%str(TurbulenceModel_l))
    if parent is None:
        node = createNode('TurbulenceModel', 'TurbulenceModel_t',
                          value=value)
    else: node = createUniqueChild(parent, 'TurbulenceModel',
                                   'TurbulenceModel_t', value=value)
    return node

# -- newThermalRelaxationModel
def newThermalRelaxationModel(value='Null', parent=None):
    """Create a new ThermalRelaxationModel node."""
    ThermalRelaxationModel_l = ['Null', 'UserDefined', 'Frozen',
                                'ThermalEquilib', 'ThermalNonequilib']
    if value not in ThermalRelaxationModel_l:
        raise ValueError('newThermalRelaxationModel: value must be in %s.'%str(ThermalRelaxationModel_l))
    if parent is None:
        node = createNode('ThermalRelaxationModel', 'ThermalRelaxationModel_t',
                          value=value)
    else: node = createUniqueChild(parent, 'ThermalRelaxationModel',
                                   'ThermalRelaxationModel_t', value=value)
    return node

# -- newChemicalKineticsModel
def newChemicalKineticsModel(value='Null', parent=None):
    """Create a new ChemicalKineticsModel node."""
    ChemicalKineticsModel_l = ['Null', 'UserDefined', 'Frozen',
                                'ChemicalEquilibCurveFit',
                                'ChemicalEquilibMinimization',
                                'ChemicalNonequilib']
    if value not in ChemicalKineticsModel_l:
        raise ValueError('newChemicalFineticsModel: value must be in %s.'%str(ChemicalKineticsModel_l))
    if parent is None:
        node = createNode('ChemicalKineticsModel', 'ChemicalKineticsModel_t',
                          value=value)
    else: node = createUniqueChild(parent, 'ChemicalKineticsModel',
                                   'ChemicalKineticsModel_t', value=value)
    return node

# -- newEMElectricFieldModel
def newEMElectricFieldModel(value='Null', parent=None):
    """Create a new EMElectricFieldModel node."""
    EMElectricFieldModel_l = ['Null', 'UserDefined', 'Constant', 'Frozen',
                              'Interpolated', 'Voltage']
    if value not in EMElectricFieldModel_l:
        raise ValueError('newEMElectricFieldModel: value must be in %s.'%str(EMElectricFieldModel_l))
    if parent is None:
        node = createNode('EMElectricFieldModel', 'EMElectricFieldModel_t',
                          value=value)
    else: node = createUniqueChild(parent, 'EMElectricFieldModel',
                                   'EMElectricFieldModel_t', value=value)
    return node

# -- newEMMagneticFieldModel
def newEMMagneticFieldModel(value='Null', parent=None):
    """Create a new EMMagneticFieldModel node."""
    EMMagneticFieldModel_l = ['Null', 'UserDefined',
                              'Constant', 'Frozen', 'Interpolated']
    if value not in EMMagneticFieldModel_l:
        raise ValueError('newEMMagneticFieldModel: value must be in %s.'%str(EMMagneticFieldModel_l))
    if parent is None:
        node = createNode('EMMagneticFieldModel', 'EMMagneticFieldModel_t',
                          value=value)
    else: node = createUniqueChild(parent, 'EMMagneticFieldModel',
                                   'EMMagneticFieldModel_t', value=value)
    return node

# -- newEMConductivityModel
def newEMConductivityModel(value='Null', parent=None):
    """Create a new EMConductivityModel node."""
    EMConductivityModel_l = ['Null', 'UserDefined', 'Constant', 'Frozen',
                             'Equilibrium_LinRessler', 'Chemistry_LinRessler']
    if value not in EMConductivityModel_l:
        raise ValueError('newEMConductivityModel: value must be in %s.'%str(EMConductivityModel_l))
    if parent is None:
        node = createNode('EMConductivityModel', 'EMConductivityModel_t',
                          value=value)
    else: node = createUniqueChild(parent, 'EMConductivityModel',
                                   'EMConductivityModel_t', value=value)
    return node

# -- newBaseIterativeData
def newBaseIterativeData(name='BaseIterativeData', nsteps=0,
                         itype='IterationValues', parent=None):
    """Create a new BaseIterativeData node."""
    if parent is None:
        node = createNode(name, 'BaseIterativeData_t', value=nsteps)
    else: node = createUniqueChild(parent, name, 'BaseIterativeData_t',
                                   value=nsteps)
    newDataArray(itype, value=numpy.arange(1,nsteps+1,dtype='int32'), parent=node)
    return node

# -- newZoneIterativeData
def newZoneIterativeData(name='ZoneIterativeData', parent=None):
    """Create a new ZoneIterativeData node."""
    if parent is None:
        node = createNode(name, 'ZoneIterativeData_t')
    else: node = createUniqueChild(parent, name, 'ZoneIterativeData_t')
    return node

# -- newRigidGridMotion
def newRigidGridMotion(name='Motion',
                       origin=[0.,0.,0.], mtype='Null', parent=None):
    """Create a new RigidMotion node."""
    if parent is None:
        node = createNode(name, 'RigidGridMotion_t')
    else: node = createUniqueChild(parent, name, 'RigidGridMotion_t')
    newDataArray('OriginLocation', value=origin, parent=node)
    newRigidGridMotionType(value=mtype, parent=node)
    return node

# -- newRigidGridMotionType
def newRigidGridMotionType(value='Null', parent=None):
    """Create a new RigidGridMotionType node."""
    RigidGridMotionType_l = ['Null', 'UserDefined', 'ConstantRate', 'VariableRate']
    if value not in RigidGridMotionType_l:
        raise ValueError('newRigidGridMotionType: value must be in %s.'%str(RigidGridMotionType_l))
    if parent is None:
        node = createNode('RigidGridMotionType', 'RigidGridMotionType_t', value=value)
    else: node = createUniqueChild(parent, 'RigidGridMotionType', 'RigidGridMotionType_t', value=value)
    return node

# -- newReferenceState
def newReferenceState(name='ReferenceState', parent=None):
    """Create a new Reference State node."""
    if parent is None:
        node = createNode(name, 'ReferenceState_t')
    else: node = createUniqueChild(parent, name, 'ReferenceState_t')
    return node

# -- newConvergenceHistory
def newConvergenceHistory(name='GlobalConvergenceHistory',
                          value=0, parent=None):
    """Create a new ConvergenceHistory node."""
    ConvergenceHistory_l = ['ZoneConvergenceHistory',
                            'GlobalConvergenceHistory']
    if name not in ConvergenceHistory_l:
        raise ValueError('newConvergenceHistory: name must be in %s.'%str(ConvergenceHistory_l))
    if parent is None:
        node = createNode(name, 'ConvergenceHistory_t', value=value)
    else: node = createUniqueChild(parent, name, 'ConvergenceHistory_t',
                                   value=value)
    return node

# -- newFamily (zones)
def newFamily(name='Family', parent=None):
    """Create a new Family node."""
    if parent is None:
        node = createNode(name, 'Family_t')
    else: node = createUniqueChild(parent, name, 'Family_t')
    return node

# -- newFamilyName (zones)
def newFamilyName(name='Family', value='myFamily', parent=None):
    """Create a new FamilyName node."""
    if parent is None:
        node = createNode(name, 'FamilyName_t', value=value)
    else: node = createUniqueChild(parent, name, 'FamilyName_t', value=value)
    return node

# -- newFamilyBC
def newFamilyBC(value='UserDefined', parent=None):
    """Create a new FamilyBC node."""
    BCType_l = ['Null', 'UserDefined', 'BCGeneral', 'BCDirichlet',
                'BCNeumann', 'BCExtrapolate', 'BCWallInviscid',
                'BCWallViscousHeatFlux', 'BCWallViscousIsothermal',
                'BCWallViscous', 'BCWall',
                'BCInflowSubsonic', 'BCInflowSupersonic',
                'BCOutflowSubsonic', 'BCOutflowSupersonic',
                'BCTunnelInflow', 'BCTunnelOutflow',
                'BCDegenerateLine', 'BCDegeneratePoint',
                'BCSymmetryPlane', 'BCSymmetryPolar',
                'BCAxisymmetricWedge', 'FamilySpecified',
                'BCFarfield', 'BCInflow', 'BCOutflow']
    if value not in BCType_l:
        raise ValueError('newFamilyBC: value must be in %s.'%str(BCType_l))
    if parent is None:
        node = createNode('FamilyBC', 'FamilyBC_t', value=value)
    else: node = createUniqueChild(parent, 'FamilyBC',
                                   'FamilyBC_t', value=value)
    return node

# -- newGeometryReference
def newGeometryReference(value='UserDefined', file='MyCAD.iges', parent=None):
    """Create a new GeometryReference node."""
    GeometryFormat_l = ['Null', 'NASA-IGES', 'SDRC', 'Unigraphics',
                        'ProEngineer', 'ICEM-CFD', 'UserDefined']
    if value not in GeometryFormat_l:
        raise ValueError('newGeometryReference: value must be in %s.'%str(GeometryFormat_l))
    if parent is None:
        node = createNode('GeometryReference', 'GeometryReference_t')
    else: node = createUniqueChild(parent, 'GeometryReference',
                                   'GeometryReference_t')
    createUniqueChild(node, 'GeometryFormat', 'GeometryFormat_t',
                      value=value)
    createUniqueChild(node, 'GeometryFile', 'GeometryFile_t',
                      value=file)
    return node

# -- newArbitraryGridMotion
def newArbitraryGridMotion(name='Motion', value='Null', parent=None):
    """Create anew ArbitrayGridMotion node."""
    ArbitraryGridMotion_l = ['Null', 'UserDefined',
                             'NonDeformingGrid', 'DeformingGrid']
    if value not in ArbitraryGridMotion_l:
        raise ValueError('newArbitraryGridMotion: value must be in %s.'%str(ArbitraryGridMotion_l))
    if parent is None:
        node = createNode(name, 'ArbitraryGridMotion_t')
    else: node = createUniqueChild(parent, name, 'ArbitraryGridMotion_t')
    createUniqueChild(node, 'ArbitraryGridMotion',
		      'ArbitraryGridMotionType_t', value=value)
    return node

# -- newUserDefinedData
def newUserDefinedData(name='UserDefined', value=None, parent=None):
    """Create a new UserDefinedData node."""
    if parent is None:
        node = createNode(name, 'UserDefinedData_t', value=value)
    else: node = createUniqueChild(parent, name, 'UserDefinedData_t',
                                   value=value)
    return node

# -- newGravity
def newGravity(value=[0.,0.,9.81], parent=None):
    """Create a new Gravity node."""
    if parent is None:
        node = createNode('Gravity', 'Gravity_t')
    else: node = createUniqueChild(parent, 'Gravity', 'Gravity_t')
    newDataArray('GravityVector', value=value, parent=node)
    return node

# -- newParentElements
def newParentElements(value=None, parent=None):
    """Create a new ParentElements node."""
    node = newDataArray('ParentElements', value=value, parent=parent)
    return node

# -- newParentElementsPosition
def newParentElementsPosition(value=None, parent=None):
    """Create a new ParentElementsPosition node."""
    node = newDataArray('ParentElementsPosition', value=value, parent=parent)
    return node

#==============================================================================
# -- Node access --
#==============================================================================

# -- Retourne le path d'un noeud
# IN: t: noeud de depart (generalement pyTree)
# IN: node: le noeud dont on cherche le chemin par rapport a t
# OUT: le chemin du noeud ou None si not found
def getPath(t, node, pyCGNSLike=False):
    """Return the path of node."""
    if t is node: return ''
    path = ''; found = []
    getPath__(t, node, path, found)
    if len(found) == 0: return None
    spath = found[0][0:-1]
    if pyCGNSLike: spath = spath.replace('CGNSTree', '')
    return spath

def getPath__(n, node, path, found):
    path += n[0]+'/'
    if n is node: 
        found.append(path); return
    for c in n[2]:
        getPath__(c, node, path, found)

# -- Retourne un noeud d'apres son path
# IN: node: noeud de depart pour la recherche
# IN: path: chemin relatif par rapport a ce noeud
def getNodeFromPath(t, path):
    """Return a node from a path."""
    if path == '' or path == '/': return t
    if path[0] == '/': path = path[1:]
    if path[-1] == '/': path = path[:-1]
    if t[0] == path: return t
    if t[3] == 'CGNSTree_t': p = path.replace(t[0]+'/','')
    else: p = path
    p = p.split('/')
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t[isStd:]:
            if c[0] == p[0]:
                if len(p) == 1: return c
                return getNodeFromPath__(c, p[1:])
    else: return getNodeFromPath__(t, p)

getNodeByPath = getNodeFromPath # alias

def getNodeFromPath__(node, path):
    for c in node[2]:
        if c[0] == path[0]:
            if len(path) == 1: return c
            return getNodeFromPath__(c, path[1:])
    return None

# -- Retourne une liste des noeuds qui ont le type 'ntype'.
# On demarre le parcours a partir de node. Parcours complet de l'arbre.
def getNodesFromType(t, ntype):
    """Return a list of nodes matching given type."""
    result = []
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t[isStd:]: getNodesFromType__(c, ntype, result)
    else: getNodesFromType__(t, ntype, result)
    return result

getNodesByType = getNodesFromType # alias

def getNodesFromType__(node, ntype, result):
    if node[3] == ntype: result.append(node)
    for c in node[2]: getNodesFromType__(c, ntype, result)

# Parcours 1 niveau de recursivite seulement
def getNodesFromType1(node, ntype):
    result = []
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getNodesFromType1__(c, ntype, result)
    else: getNodesFromType1__(node, ntype, result)
    return result

getNodesByType1 = getNodesFromType1 # alias

def getNodesFromType1__(node, ntype, result):
    if node[3] == ntype: result.append(node)
    for c in node[2]:
        if c[3] == ntype: result.append(c)

# Parcours 2 niveaux de recursivite seulement
def getNodesFromType2(node, ntype):
    result = []
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getNodesFromType2__(c, ntype, result)
    else: getNodesFromType2__(node, ntype, result)
    return result

getNodesByType2 = getNodesFromType2 # alias

def getNodesFromType2__(node, ntype, result):
    if node[3] == ntype: result.append(node)
    for c in node[2]:
        if c[3] == ntype: result.append(c)
        for d in c[2]:
            if d[3] == ntype: result.append(d)

# Parcours 3 niveaux de recursivite seulement
def getNodesFromType3(node, ntype):
    result = []
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getNodesFromType3__(c, ntype, result)
    else: getNodesFromType3__(node, ntype, result)
    return result

getNodesByType3 = getNodesFromType3 # alias

def getNodesFromType3__(node, ntype, result):
    if node[3] == ntype: result.append(node)
    for c in node[2]:
        if c[3] == ntype: result.append(c)
        for d in c[2]:
            if d[3] == ntype: result.append(d)
            for e in d[2]:
                if e[3] == ntype: result.append(e)

def getByType(t, ntype, recursive=-1):
    """Return a standard node containing the list of matching type nodes as children."""
    if recursive == 1: ret = getNodesFromType1(t, ntype)
    elif recursive == 2: ret = getNodesFromType2(t, ntype)
    elif recursive == 3: ret = getNodesFromType3(t, ntype)
    else: ret = getNodesFromType(t, ntype)
    return [ntype, None, ret, None]

# -- Retourne un seul noeud (no wildcard) - Fast
def getNodeFromType(t, ntype):
    """Return the first matching node with given type."""
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t[isStd:]:
            ret = getNodeFromType(c, ntype)
            if ret is not None: return ret
        return None
    else: return getNodeFromType__(t, ntype)

getNodeByType = getNodeFromType # alias

def getNodeFromType__(node, ntype):
    if node[3] == ntype: return node
    for i in node[2]:
        ret = getNodeFromType(i, ntype)
        if ret is not None: return ret
    return None

# Parcours un seul niveau de recursivite, node doit etre un noeud du pyTree
def getNodeFromType1(node, ntype):
    if node[3] == ntype: return node
    for c in node[2]:
        if c[3] == ntype: return c
    return None

# node doit etre un noeud du pyTree
def getNodeFromType2(node, ntype):
    if node[3] == ntype: return node
    for c in node[2]:
        if c[3] == ntype: return c
        for d in c[2]:
            if d[3] == ntype: return d
    return None

getNodeByType2 = getNodeFromType2 # alias

# node doit etre un noeud du pyTree
def getNodeFromType3(node, ntype):
    if node[3] == ntype: return node
    for c in node[2]:
        if c[3] == ntype: return c
        for d in c[2]:
            if d[3] == ntype: return d
            for e in d[2]:
                if e[3] == ntype: return e
    return None

getNodeByType3 = getNodeFromType3 # alias

# -- Retourne une liste des chemins qui ont le type 'ntype'.
# On demarre le parcours a partir de node. Parcours complet de l'arbre.
# Si pyCGNSLike=True, le chemin ne comporte pas de /CGNSTree
def getPathsFromType(node, ntype, pyCGNSLike=False):
    """Return a list of paths corresponding to a given type."""
    result = []; current = ''
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getPathsFromType__(c, ntype, current, result)
    else:
        getPathsFromType__(node, ntype, current, result)
    return fixPaths__(result, int(pyCGNSLike))

getPathsByType = getPathsFromType # alias

def getPathsFromType__(node, ntype, current, result):
    current += '/'+node[0]
    if node[3] == ntype: result.append(current)
    for c in node[2]: getPathsFromType__(c, ntype, current, result)

def getPathsFromType1(node, ntype, pyCGNSLike=False):
    result = []; current = ''
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getPathsFromType1__(c, ntype, current, result)
    else:
        getPathsFromType1__(node, ntype, current, result)
    return fixPaths__(result, int(pyCGNSLike))

getPathsByType1 = getPathsFromType1 # alias

def getPathsFromType1__(node, ntype, current, result):
    current += '/'+node[0]
    if node[3] == ntype: result.append(current)
    for c in node[2]:
        if c[3] == ntype: result.append(current+'/'+c[0])

def getPathsFromType2(node, ntype, pyCGNSLike=False):
    result = []; current = ''
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getPathsFromType2__(c, ntype, current, result)
    else:
        getPathsFromType2__(node, ntype, current, result)
    return fixPaths__(result, int(pyCGNSLike))

getPathsByType2 = getPathsFromType2 # alias

def getPathsFromType2__(node, ntype, current, result):
    current += '/'+node[0]
    if node[3] == ntype: result.append(current)
    for c in node[2]:
        if c[3] == ntype: result.append(current+'/'+c[0])
        for d in c[2]:
            if d[3] == ntype: result.append(current+'/'+c[0]+'/'+d[0])

# ntype=0: /CGNSTree/Base/Zone
# ntype=1: /Base/Zone (pyCGNSLike=True)
# ntype=2: Base/Zone
def fixPaths__(result, ntype=0):
    c = 0
    if ntype == 1:
        for i in result: 
            result[c] = i.replace('/CGNSTree', ''); c += 1
    elif ntype == 2:
        for i in result:
            result[c] = i.replace('/CGNSTree/', ''); c += 1
    else:
        for i in result: result[c] = i[1:]; c += 1
    return result

# -- Retourne les chemins des zones
def getZonePaths(t, pyCGNSLike=False):
    """Return the list of paths of zones."""
    result = []; current1 = ''; ntype = 'Zone_t'
    current1 += '/'+t[0]
    if t[3] == ntype: result.append(current1)
    for c in t[2]:
        current2 = current1+'/'+c[0]
        if c[3] == ntype: result.append(current2)
        for d in c[2]:
            current3 = current2+'/'+d[0]
            if d[3] == ntype: result.append(current3)
    return fixPaths__(result, int(pyCGNSLike))

# -- Retourne une liste des chemins qui ont le nom 'name'.
# On demarre le parcours a partir de node. Parcours complet de l'arbre.
# Si pyCGNSLike=True, le chemin ne comporte pas de CGNSTree
def getPathsFromName(node, name, pyCGNSLike=False):
    """Return a list of paths corresponding to a given name."""
    result = []; current = ''
    isStd = isStdNode(node)
    if isStd >= 0:
        if ('*' in name)|('?' in name)|('[' in name):
            for c in node[isStd:]: getPathsFromName___(c, name, current, result)
        else:
            for c in node[isStd:]: getPathsFromName__(c, name, current, result)
    else:
        if ('*' in name)|('?' in name)|('[' in name):
            getPathsFromName___(node, name, current, result)
        else:
            getPathsFromName__(node, name, current, result)
    return fixPaths__(result, int(pyCGNSLike))

getPathsByName = getPathsFromName # alias

def getPathsFromName__(node, name, current, result):
    current += '/'+node[0]
    if node[0] == name: result.append(current)
    for c in node[2]: getPathsFromName__(c, name, current, result)

def getPathsFromName___(node, name, current, result):
    current += '/'+node[0]
    if fnmatch.fnmatch(node[0], name): result.append(current)
    for c in node[2]: getPathsFromName___(c, name, current, result)

def getPathsFromName1(node, name, pyCGNSLike=False):
    result = []; current = ''
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getPathsFromName1__(c, name, current, result)
    else:
        getPathsFromName1__(node, name, current, result)
    return fixPaths__(result, int(pyCGNSLike))

getPathsByName1 = getPathsFromName1 # alias

def getPathsFromName1__(node, name, current, result):
    current += '/'+node[0]
    if node[0] == name: result.append(current)
    for c in node[2]:
        if c[0] == name: result.append(c)

def getPathsFromName2(node, name, pyCGNSLike=False):
    result = []; current = ''
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getPathsFromName2__(c, name, current, result)
    else:
        getPathsFromName2__(node, name, current, result)
    return fixPaths__(result, int(pyCGNSLike))

getPathsByName2 = getPathsFromName2 # alias

def getPathsFromName2__(node, name, current, result):
    current += '/'+node[0]
    if node[0] == name: result.append(current)
    for c in node[2]:
        if c[0] == name: result.append(current+'/'+c[0])
        for d in c[2]:
            if d[0] == name: result.append(current+'/'+c[0]+'/'+d[0])

# -- Retourne une liste des chemins qui ont la valeur value.
# On demarre le parcours a partir de node. Parcours complet de l'arbre.
# Si pyCGNSLike=True, le chemin ne comporte pas de /CGNSTree
def getPathsFromValue(node, value, pyCGNSLike=False):
    """Return a list of paths corresponding to a given value."""
    result = []; current = ''
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getPathsFromValue__(c, value, current, result)
    else:
        getPathsFromValue__(node, value, current, result)
    return fixPaths__(result, int(pyCGNSLike))

getPathsByValue = getPathsFromValue # alias

def getPathsFromValue__(node, value, current, result):
    current += '/'+node[0]
    if isValue(node, value): result.append(current)
    for c in node[2]: getPathsFromValue__(c, value, current, result)

# -- Retourne la fin d'un chemin --
def getPathLeaf(path):
    """Return end of path."""
    s = path.split('/')
    if len(s) > 0: return s[-1]
    else: return ''

# -- Retourne l'ancetre de nom --
def getPathAncestor(path, level=1):
    """Return the path ancestor of path."""
    if level == 0: return path
    s = path.split('/'); ls =  len(s); c = 0
    for i in range(level): c += len(s[ls-i-1])+1
    return path[0:-c]

# -- Retourne les noeuds Zone_t --
def getZones(t):
    """Return a list of all Zone_t nodes."""
    result = []; ntype = 'Zone_t'
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t[isStd:]: getNodesFromType2__(c, ntype, result)
    else: getNodesFromType2__(t, ntype, result)
    return result

# -- getZonesPerIteration --
# Si pas d'argument iteration ou time, retourne une liste de zones pour chaque iteration (a plat?)
# Si argument time ou iteration, retourne la liste de zones correspondantes
def getZonesPerIteration(t, iteration=None, time=None):
    """Return zones for given iteration of time."""
    bid = getNodeFromType2(t, 'BaseIterativeData_t')
    if bid is None:
        raise ValueError('getZonesPerIteration: no BaseIterativeData.')

    nbOfZones = getNodeFromName1(bid, 'NumberOfZones')
    zonePtrs = getNodeFromName1(bid, 'ZonePointers')
    if nbOfZones is None:
        raise ValueError('getZonesPerIteration: no NumberOfZones node in BaseIterativedata.')
    if zonePtrs is None:
        raise ValueError('getZonesPerIteration: no ZonePointers node in BaseIterativeData.')

    nbOfZones = nbOfZones[1]
    zonePtrs = zonePtrs[1]
    
    if iteration is not None:
        zoneNames = []
        for nbz in range(nbOfZones[iteration]):
            zoneName = zonePtrs[:, nbz, iteration].tostring().strip()
            zoneNames.append(zoneName)
        return [getNodeFromName2(t,z) for z in zoneNames]

    if time is not None:
        timeValues = getNodeFromName1(bid, 'TimeValues')[1]
        zoneNames = []
        for i in range(nbOfZones.shape[0]):
            if time == timeValues[i]:
                for nbz in range(nbOfZones[i]):
                    zoneName = zonePtrs[:, nbz, i].tostring().strip()
                    zoneNames.append(zoneName)
        return [getNodeFromName2(t,z) for z in zoneNames]

    zonesPerIteration = []
    for i in range(nbOfZones.shape[0]):
        zones = []
        for nbz in range(nbOfZones[i]):
            zoneName = zonePtrs[:, nbz, i].tostring().strip()
            zones.append(getNodeFromName2(t,zoneName))
        zonesPerIteration.append(zones)
    return zonesPerIteration

# -- Retourne les noeuds CGNSBase_t --
def getBases(t):
    """Return a list of all CGNSBase_t nodes."""
    result = []; ntype = 'CGNSBase_t'
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t[isStd:]: getNodesFromType1__(c, ntype, result)
    else: getNodesFromType1__(t, ntype, result)
    return result

# -- Retourne une liste des noeuds qui ont le nom 'name'.
# On demarre le parcours a partir de node.
# Parcours tout l'arbre. Wildcards possibles.
def getNodesFromName(t, name):
    """Return a list of nodes matching given name."""
    result = []
    isStd = isStdNode(t)
    if isStd >= 0:
        if ('*' in name)|('?' in name)|('[' in name):
            for c in t[isStd:]: getNodesFromName___(c, name, result)
        else:
            for c in t[isStd:]: getNodesFromName__(c, name, result)
    else:
        if ('*' in name)|('?' in name)|('[' in name):
            getNodesFromName___(t, name, result)
        else: getNodesFromName__(t, name, result)
    return result

getNodesByName = getNodesFromName # alias

def getNodesFromName__(node, name, result):
    if node[0] == name: result.append(node)
    for c in node[2]: getNodesFromName__(c, name, result)

def getNodesFromName___(node, name, result):
    if fnmatch.fnmatch(node[0], name): result.append(node)
    for c in node[2]: getNodesFromName___(c, name, result)

# Parcours 1 niveau de recursivite seulement (no widcard)
def getNodesFromName1(node, name):
    result = []
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getNodesFromName1__(c, name, result)
    else: getNodesFromName1__(node, name, result)
    return result

getNodesByName1 = getNodesFromName1 # alias

def getNodesFromName1__(node, name, result):
    if node[0] == name: result.append(node)
    for c in node[2]:
        if c[0] == name: result.append(c)

# Parcours 2 niveaux de recursivite seulement
def getNodesFromName2(node, name):
    result = []
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getNodesFromName2__(c, name, result)
    else: getNodesFromName2__(node, name, result)
    return result

getNodesByName2 = getNodesFromName2 # alias

def getNodesFromName2__(node, name, result):
    if node[0] == name: result.append(node)
    for c in node[2]:
        if c[0] == name: result.append(c)
        for d in c[2]:
            if d[0] == name: result.append(d)

# Parcours 3 niveaux de recursivite seulement
def getNodesFromName3(node, name):
    result = []
    isStd = isStdNode(node)
    if isStd >= 0:
        for c in node[isStd:]: getNodesFromName3__(c, name, result)
    else: getNodesFromName3__(node, name, result)
    return result

getNodesByName3 = getNodesFromName3 # alias

def getNodesFromName3__(node, name, result):
    if node[0] == name: result.append(node)
    for c in node[2]:
        if c[0] == name: result.append(c)
        for d in c[2]:
            if d[0] == name: result.append(d)
            for e in d[2]:
                if e[0] == name: result.append(e)

def getByName(t, name, recursive=-1):
    """Return a standard node containing the list of matching name nodes as children."""
    if recursive == 1: ret = getNodesFromName1(t, name)
    elif recursive == 2: ret = getNodesFromName2(t, name)
    elif recursive == 3: ret = getNodesFromName3(t, name)
    else: ret = getNodesFromName(t, name)
    return [name, None, ret, None]

# -- Retourne un seul noeud (no wildcard) - Fast
def getNodeFromName(t, name):
    """Return the first matching node with given name."""
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t[isStd:]:
            ret = getNodeFromName(c, name)
            if ret is not None: return ret
        return None
    else: return getNodeFromName__(t, name)

getNodeByName = getNodeFromName # alias

def getNodeFromName__(node, name):
    if node[0] == name: return node
    for i in node[2]:
        ret = getNodeFromName(i, name)
        if ret is not None: return ret
    return None

# Parcours un seul niveau de recursivite, node doit etre un noeud du pyTree
def getNodeFromName1(node, name):
    if node[0] == name: return node
    for c in node[2]:
        if c[0] == name: return c
    return None

getNodeByName1 = getNodeFromName1 # alias

# node doit etre un noeud du pyTree
def getNodeFromName2(node, name):
    if node[0] == name: return node
    for c in node[2]:
        if c[0] == name: return c
        for d in c[2]:
            if d[0] == name: return d
    return None

getNodeByName2 = getNodeFromName2 # alias

# node doit etre un noeud du pyTree
def getNodeFromName3(node, name):
    if node[0] == name: return node
    for c in node[2]:
        if c[0] == name: return c
        for d in c[2]:
            if d[0] == name: return d
            for e in d[2]:
                if e[0] == name: return e
    return None

getNodeByName3 = getNodeFromName3 # alias

# -- Retourne une liste des noeuds qui ont la valeur 'value'.
# On demarre le parcours a partir de node.
def getNodesFromValue(t, value):
    """Return a list of nodes matching given value."""
    result = []
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t[isStd:]: getNodesFromValue__(c, value, result)
    else: getNodesFromValue__(t, value, result)
    return result

getNodesByValue = getNodesFromValue # alias

def getNodesFromValue__(node, value, result):
    if isValue(node, value): result.append(node)
    for c in node[2]: getNodesFromValue__(c, value, result)

# -- Retourne une liste des noeuds qui ont le nom 'name' et le type 'ntype'
# On demarre le parcours a partir de node.
# Parcours tout l'arbre. Wildcards possibles.
def getNodesFromNameAndType(t, name, ntype):
    """Return a list of nodes matching given name and type."""
    result = []
    isStd = isStdNode(t)
    if isStd >= 0:
        if ('*' in name)|('?' in name)|('[' in name) or ('*' in ntype)|('?' in ntype)|('[' in ntype):
            for c in t[isStd:]: getNodesFromNameAndType___(c, name, ntype, result)
        else:
            for c in t[isStd:]: getNodesFromNameAndType__(c, name, ntype, result)
    else:
        if ('*' in name)|('?' in name)|('[' in name) or ('*' in ntype)|('?' in ntype)|('[' in ntype):
            getNodesFromNameAndType___(t, name, ntype, result)
        else: getNodesFromNameAndType__(t, name, ntype, result)
    return result

def getNodesFromNameAndType__(node, name, ntype, result):
    if node[0] == name and node[3] == ntype: result.append(node)
    for c in node[2]: getNodesFromNameAndType__(c, name, ntype, result)

def getNodesFromNameAndType___(node, name, ntype, result):
    if fnmatch.fnmatch(node[0], name) and fnmatch.fnmatch(node[3], ntype): result.append(node)
    for c in node[2]: getNodesFromNameAndType___(c, name, ntype, result)

# -- Retourne un seul noeud (no wildcard) - Fast
def getNodeFromNameAndType(t, name, ntype):
    """Return the first matching node with given name and type."""
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t[isStd:]:
            ret = getNodeFromNameAndType(c, name, ntype)
            if ret is not None: return ret
        return None
    else: return getNodeFromNameAndType__(t, name, ntype)

def getNodeFromNameAndType__(node, name, ntype):
    if node[0] == name and node[3] == ntype: return node
    for i in node[2]:
        ret = getNodeFromNameAndType(i, name, ntype)
        if ret is not None: return ret
    return None

# -- (p, c) = getParentOfNode(start, node). Retourne :
# Si start est un noeud standard, retourne le noeud parent de node et le no
# dans la liste des fils. p[2][c] = node.
# Si start est une liste de noeuds standards. Si node est directement dans
# cette liste, retourne la liste et le no de node: p[c] = node.
# Sinon, retourne le noeud parent de node: p[2][c] = node.
# De facon generale:
# if (isStdNode(start) == 0 and id(start) == id(p)): p[c] = node
#                                              else: p[2][c] = node
# Remarque: start doit etre au dessus de node
def getParentOfNode(t, node):
    """Return the parent of given node."""
    idNode = id(node)
    isStd = isStdNode(t)
    if isStd >= 0:
        d = isStd
        for c in t[isStd:]:
            if id(c) == idNode: return (t, d)
            d += 1
        r = None; d = 0
        for c in t[isStd:]:
            (r, d) = getParentOfNode__(c, idNode)
            if r is not None: return (r, d)
        return (None, 0)
    else: return getParentOfNode__(t, idNode)

def getParentOfNode__(start, idNode):
    c = 0; r = None; d = 0
    for i in start[2]:
        if idNode == id(i): return (start, c)
        else:
            (r, d) = getParentOfNode__(i, idNode)
            if r is not None: return (r, d)
        c += 1
    return (r, d)

def getParentFromType(start, node, parentType, prev=None):
    """Return thee first parent node matching type."""
    if id(start) == id(node): return prev
    if start[3] == parentType: prev = start
    for n in start[2]:
        ret = getParentFromType(n, node, parentType, prev)
        if ret is not None: return ret
    return None

def getParentsFromType(start, node, parentType, l=[]):
    """Return all parent nodes matching type."""
    l = l[:]
    if id(start) == id(node): return l
    if start[3] == parentType: l.append(start)#; l = l[:]
    for n in start[2]:
        ret = getParentsFromType(n, node, parentType, l)
        if ret != []: return ret
    return []
    
# -- Return the position of node in parent children list (Return -1 if not found)
def getNodePosition(node, parent):
    """Return the position of node in parent children list."""
    for c, n in enumerate(parent[2]):
        if id(n) == id(node): return c
    return -1

# -- Retourne le nom d'un noeud
def getName(node):
    """Return node name."""
    return node[0]

# -- Retourne le type d'un noeud
def getType(node):
    """Return node type."""
    return node[3]

# -- Retourne les enfants d'un noeud
def getChildren(node):
    """Return children list of a node."""
    return node[2]

# -- Retourne la valeur d'un noeud si celui-ci contient un entier, un float
# ou une chaine. Sinon retourne le tableau
def getValue(node):
    """Return the value of a node."""
    n = node[1]
    if isinstance(n, numpy.ndarray):
        if n.dtype.char == 'S': 
            if version_info[0] == 2: return n.tostring()
            else: return n.tostring().decode()
        elif n.dtype.char == 'c':
            if len(n.shape) == 1:
                if version_info[0] == 2: return n.tostring()
                else: return n.tostring().decode() 
            out = []
            for i in range(n.shape[1]):
                if version_info[0] == 2: v = n[:,i].tostring()
                else: v = n[:,i].tostring().decode()
                out.append(v.strip())
            return out
        elif n.dtype == numpy.int32:
            if n.size == 1: return int(n.flat[0])
            else: return n
        elif n.dtype == numpy.float64:
            if n.size == 1: return float(n.flat[0])
            else: return n
        elif n.dtype == numpy.float32:
            if n.size == 1: return float(n.flat[0])
            else: # cast en float64
                b = numpy.empty(n.shape, dtype=numpy.float64)
                b[:] = n[:]
                return b
    else: return node[1]

#==============================================================================
# -- rm Nodes --
#==============================================================================

# -- rmNode (t is modified)
def rmNode(t, node):
    """Remove given node from t."""
    (p, c) = getParentOfNode(t, node)
    if p is not None:
        if isStdNode(t) == 0 and id(p) == id(t): del p[c]
        else: del p[2][c]
    return None

def _rmNode(t, node):
    (p, c) = getParentOfNode(t, node)
    if p is not None:
        if isStdNode(t) == 0 and id(p) == id(t): del p[c]
        else: del p[2][c]
    return None

# -- rmNodeByPath
def rmNodeByPath(t, path):
    """Remove node by specifying its path."""
    tp = copyRef(t)
    _rmNodeByPath(tp, path)
    return tp

rmNodeFromPath = rmNodeByPath # alias

def _rmNodeByPath(t, path):
    if path == '' or path == '/': t = None; return None
    if path[0] == '/': path = path[1:]
    if path[-1] == '/': path = path[:-1]
    if t[0] == path: t = None; return None
    if t[3] == 'CGNSTree_t': p = path.replace(t[0]+'/','')
    else: p = path
    p = p.split('/')
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t: rmNodeByPath__(c, p[1:])
    else: rmNodeByPath__(t, p)
    return None

def rmNodeByPath__(t, path):
    children = list(range(len(t[2])-1,-1,-1))
    for ichild in children:
        if t[2][ichild][0] == path[0]:
            if len(path) == 1: t[2].pop(ichild)
            else: rmNodeByPath__(t[2][ichild], path[1:])
    return None

_rmNodeFromPath = _rmNodeByPath # alias

# -- rmNodesByName
def rmNodesByName(t, name):
    """Remove nodes of given name."""
    tp = copyRef(t)
    _rmNodesByName(tp, name)
    return tp

rmNodesFromName = rmNodesByName # alias

def _rmNodesByName(t, name):
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t: rmNodesByName__(c, name)
    else: rmNodesByName__(t, name)
    return None

_rmNodesFromName = _rmNodesByName # alias

def rmNodesByName__(t, name):
    children = list(range(len(t[2])-1,-1,-1))
    for ichild in children:
        if isName(t[2][ichild], name): t[2].pop(ichild)
        else: rmNodesByName__(t[2][ichild], name)
    return None

# -- rmNodesByType
def rmNodesByType(t, ntype):
    """Remove nodes of given type."""
    tp = copyRef(t)
    _rmNodesByType(tp, ntype)
    return tp

rmNodesFromType = rmNodesByType # alias

def _rmNodesByType(t, ntype):
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t: rmNodesByType__(c, ntype)
    else: rmNodesByType__(t, ntype)
    return None

_rmNodesFromType = _rmNodesByType # alias

def rmNodesByType__(t, ntype):
    nodes = getNodesFromType(t, ntype)
    for n in nodes: _rmNode(t, n)
    return None

# -- rmNodesByNameAndType
def rmNodesByNameAndType(t, name, ntype):
    """Remove nodes of that match given name and given type at the same time."""
    tp = copyRef(t)
    _rmNodesByNameAndType(tp, name, ntype)
    return tp

rmNodesFromPathAndType = rmNodesByNameAndType # alias

def _rmNodesByNameAndType(t, name, ntype):
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t: rmNodesByNameAndType__(c, name, ntype)
    else: rmNodesByNameAndType__(t, name, ntype)
    return None

_rmNodesFromPathAndType = _rmNodesByNameAndType # alias

def rmNodesByNameAndType__(t, name, ntype):
    nodes = getNodesFromType(t, ntype)
    for n in nodes:
        if isName(n, name): _rmNode(t, n)
    return None

# -- rmNodesByValue
def rmNodesByValue(t, value):
    """Remove nodes that patch given value."""
    tp = copyRef(t)
    _rmNodesByValue(tp, value)
    return tp

rmNodesFromValue = rmNodesByValue # alias

def _rmNodesByValue(t, value):
    isStd = isStdNode(t)
    if isStd >= 0:
        for c in t: rmNodesByValue__(c, value)
    else: rmNodesByValue__(t, value)
    return None

def rmNodesByValue__(t, value):
    indchild = list(range(len(t[2])-1,-1,-1))
    for ichild in indchild:
        if isValue(t[2][ichild], value): t[2].pop(ichild)
        else: rmNodesByValue__(t[2][ichild], value)
    return None

_rmNodesFromValue = _rmNodesByValue # alias

#==============================================================================
# -- Node name modifications --
#==============================================================================

# -- Sort nodes by name
def sortByName(t, recursive=True):
    """Sort nodes by their names (alphabetical order)."""
    tp = copyRef(t)
    _sortByName(tp, recursive)
    return tp

def _sortByName(t, recursive=True):
    names = [n[0] for n in t[2]]
    names.sort()
    for name in names:
        node = getNodeFromName1(t, name)
        _rmNode(t, node)
        addChild(t, node, pos=names.index(node[0]))
        if recursive: _sortByName(node, recursive)
    return None

# -- Rename node ou nodes
def renameNode(t, name, newName, subName=None):
    """Rename nodes (propagating changes)."""
    tp = copyRef(t)
    _renameNode(tp, name, newName, subName)
    return tp

def _renameNode(t, name, newName, subName=None):
    if ('*' in name)|('?' in name)|('[' in name):
        _renameNodeRe(t, name, newName, subName)
    else:
        nodes = getNodesFromName(t, name)
        nodesValue = getNodesFromValue(t, name)
        for node in nodes: setName(node, newName)
        for nodeval in nodesValue: setValue(nodeval, newName)
    return None

# name: must be a pattern
# newSubName: morceau de nom de remplacement
# subName: morceau de nom a remplacer
def _renameNodeRe(t, name, newSubName, subName):
    if subName is None: # not flawless
        subName = name.replace('*', '')
        subName = subName.replace('?', '')
    nodes = getNodesFromName(t, name)
    nodesValue = getNodesFromValue(t, name)
    for node in nodes:
        v = getName(node)
        n = v.replace(subName, newSubName)
        setName(node, n)
    for nodeval in nodesValue:
        v = getValue(nodeval)
        n = v.replace(subName, newSubName)
        setValue(nodeval, n)
    return None

# -- Append base name to zone name
def appendBaseName2ZoneName(t, updateRef=True, separator='_', trailing=''):
    """Append base name to all zone names."""
    tp = copyRef(t)
    _appendBaseName2ZoneName(tp, updateRef, separator, trailing)
    return tp

def _appendBaseName2ZoneName(t, updateRef=True, separator='_', trailing=''):
    bases = getBases(t)
    st = '%s'+separator+'%s'+trailing
    for b in bases:
        basename = b[0]
        zones = getNodesFromType1(b, 'Zone_t')
        for z in zones:
            oldname = z[0]
            newname = st%(basename, oldname)
            setName(z, newname)
            if updateRef:
                # Change les references a cette zone dans la base
                nodes = getNodesFromValue(b, oldname)
                for n in nodes: setValue(n, newname)
    return None

# -- Enleve les / et les remplace par \, car ADF et HDF ne supporte pas les /
def _adaptZoneNamesForSlash(t):
  zones = getZones(t)
  for z in zones: z[0] = z[0].replace('/', '\\')
  return None

#==============================================================================
# -- node copy/merge --
#==============================================================================

# -- Duplique un arbre ou un sous-arbre par references
def duptree__(node, parent):
    d = [node[0], node[1], [], node[3]]
    if len(parent) == 4: parent[2].append(d)
    for i in node[2]: duptree__(i, d)
    return d

# -- Copy un arbre en gardant des references sur les numpy
def copyRef(node):
    """Copy a tree sharing node values."""
    ret = isStdNode(node)
    if ret == -1:
        dup = duptree__(node, []); return dup
    elif ret == 0:
        l = list(node); lg = len(l)
        for i in range(lg): l[i] = duptree__(l[i], [])
        return l
    else: return node

# -- Copie un arbre ou un sous-arbre en copiant aussi les numpy.array
def copyTree(node, parent=None, order='F'):
    """Fully copy a tree."""
    if node[1] is not None and isinstance(node[1], numpy.ndarray):
        d = [node[0], node[1].copy(order), [], node[3]]
    else: d = [node[0], node[1], [], node[3]]
    if parent is not None and len(parent) == 4: parent[2].append(d)
    for i in node[2]: copyTree(i, d, order=order)
    return d

# -- Copie un noeud (pas de recursivite)
def copyNode(node):
    """Copy only this node (no recursion)."""
    if node[1] is not None and isinstance(node[1], numpy.ndarray):
        d = [node[0], node[1].copy('F'), node[2], node[3]]
    else: d = [node[0], node[1], node[2], node[3]]
    return d

# -- copyOnly(node) --
# Copie recursivement les noms (si names=True), les types (si types=True)
# avec les criteres de nom=byName et/ou type=byType)
def duptree1__(node, byName, byType, parent):
    value = node[1]
    if byName is not None:
        if fnmatch.fnmatch(node[0], byName):
            value = node[1].copy('F')
    if byType is not None:
        if fnmatch.fnmatch(node[3], byType):
            value = node[1].copy('F')
    if byName is None and byType is None:
        value = node[1].copy('F')
    d = [node[0], value, [], node[3]]
    if len(parent) == 4: parent[2].append(d)
    for i in node[2]: duptree1__(i, byName, byType, d)
    return d

def copyValue(node, byName=None, byType=None):
    """Copy the value of nodes specified by byName or byType string."""
    ret = isStdNode(node)
    if ret == -1:
        dup = duptree1__(node, byName, byType, []); return dup
    elif ret == 0:
        l = list(node); lg = len(l)
        for i in range(lg): l[i] = duptree1__(l[i], byName, byType, [])
        return l
    else: return node

# -- Merge nodes
# Merge une liste d'arbres en un seul. Noeud a noeud, niveaux a niveaux,
# si un noeud existe plusieurs fois, on prend le premier dans la liste.
def merge(A, pathList=None):
    """Merge a list of pyTrees in one."""
    # check
    if pathList:
        for path in pathList:
            if not isinstance(path,str): raise AttributeError("merge: pathList should be a list of strings.")
    tp = copyRef(A[0])
    if pathList is None: pathList = [None for tree in A[1:]]
    for path, tree in zip(pathList, A[1:]):
        if path: _append(tp, tree, path)
        else:
            if tp[0] != tree[0]: raise AttributeError("merge: path is None, but root name are differents.")
            for child in tree[2]: _append(tp, child, tree[0])
    return tp

#==============================================================================
# -- print nodes --
#==============================================================================

def repr__(t, Out=None, DEB='', LAST=False, code=['','','']):
    if Out is None: Out = []; First = True
    else: First = False
    out = "['%s%s%s',"%(code[0],t[0],code[2])
    if isinstance(t[1], str):
        out += "'%s',["%t[1]
    elif isinstance(t[1], (int,float,list)):
        out += str(t[1])+",["
    elif isinstance(t[1],numpy.ndarray):
        if t[1].dtype==numpy.array('a').dtype and t[1].size < 32:
            out += "array('%s',dtype='%s'),["%(t[1].tostring(),t[1].dtype)
        elif t[1].size == 1 :
            out += "array(%s,dtype='%s'),["%(t[1].tolist(),t[1].dtype)
        else:
            Flag = 'C'
            if t[1].flags.fortran: Flag='F'
            out += "array(shape=%s,dtype='%s',order='%s'),["%(t[1].shape,t[1].dtype,Flag)
    elif t[1] is None: out += 'None,['
    else: out += '%s value is not printable,['%(type(t[1]))
    if len(t[2]) < 2: out += "%d son],'%s%s%s']"%(len(t[2]),code[1],t[3],code[2])
    else: out += "%d sons],'%s%s%s']"%(len(t[2]),code[1],t[3],code[2])
    if DEB != '': Out.append(DEB+'_'+out+'\n')
    else: Out.append(out+'\n')
    if LAST: DEB = DEB[:-1]+' '
    outDEB = DEB+3*' '+'|'
    for ichild in range(len(t[2])):
        if ichild == (len(t[2])-1):
            repr__(t[2][ichild], Out, outDEB, LAST=True, code=code)
        else:
            repr__(t[2][ichild], Out, outDEB, code=code)
    if First: return Out

def printTree(node, file=None, stdOut=None, editor=None, color=False):
    """Pretty print a pyTree node."""
    if color: code = ['\033[34m', '\033[32m', '\033[0m']
    else: code = ['', '', '']
    stdNode = isStdNode(node)
    if stdNode == -2: return
    elif stdNode == -1: rep = ''.join(repr__(node, code=code))
    else:
        out = []
        for b in node: out += repr__(b)
        rep = ''.join(out)

    if stdOut: stdOut.write(rep)
    if editor:
        from tempfile import mkstemp
        import os
        fd, tmpfl = mkstemp('.py')
        os.write(fd, rep)
        os.close(fd)
        os.system('%s %s ;rm -f %s' %(editor, tmpfl, tmpfl))
    if file:
        fi = open(file, 'w')
        fi.write(rep)
        fi.close()
    if file is None and editor is None and stdOut is None:
        import sys
        sys.stdout.write(rep)

# Mesure de la taille de a en octets
def getSizeOf__(a, s):
    s += len(a[0])
    s += len(a[3])
    r = a[1]
    if r is not None:
        if isinstance(r, numpy.ndarray):
            if r.dtype == numpy.int32: s += r.size*4
            else: s += r.size*8
    for i in a[2]:
        s = getSizeOf__(i, s)
    return s

def getSizeOf(a):
    """Return the size of a in octets."""
    s = 0
    if isStdNode(a) == 0:
        for i in a: 
            sl = 0
            s += getSizeOf__(i, sl)
    else: s = getSizeOf__(a, s)
    return s

#==============================================================================
# -- Conversion zones, bases, listes de zones <-> tree + noms + ranges --
#==============================================================================

# -- Converti un noeud zone, liste de zones, base, liste de bases,
# tree en pyTree (avec les memes adresses)
# Retourne (t, ntype) avec ntype=1 (zone), 2 (liste de zones), 3 (tree),
# 4 (base), 5 (liste de bases)
def node2PyTree(node):
    ntype = typeOfNode(node)
    if ntype == 1: # zone
        t = newCGNSTree()
        dim = getZoneDim(node)[4]
        base = createBaseNode('Base', dim)
        t[2].append(base)
        t[2][1][2].append(node)
    elif ntype == 2: # une liste de zones
        t = newCGNSTree()
        dim = getZoneDim(node[0])[4]
        base = createBaseNode('Base', dim)
        t[2].append(base)
        t[2][1][2] += node
    elif ntype == 4: # base
        t = newCGNSTree()
        t[2] += [node]
    elif ntype == 5: # liste de bases
        t = newCGNSTree()
        t[2] += node
    else: t = node
    return t, ntype

# -- Converti un pyTree en noeud de type fourni (avec les memes adresses)
# IN: type=1 (zone), 2 (liste de zones), 3 (tree), 4 (base), 5 (liste de bases)
def pyTree2Node(t, type):
    if type == 1: # zone: renvoit une zone ou une liste de zones
        node = t[2][1][2]
        if len(node) == 1: node = node[0]
    elif type == 2: # liste de zones: renvoit une liste de zones
        node = t[2][1][2]
    elif type == 4: # base: renvoit une base
        node = t[2][1]
    elif type == 5: # liste de bases: renvoit une liste de bases
        node = t[2][1:]
    else: node = t
    return node

# -- EltName2EltNo
# Convertit un nom CGNS d'elt en no CGNS d'elt et son nombre de noeuds associes
def eltName2EltNo(name):
    eltno = 0; nnodes = 0
    if name == 'NODE':
        eltno = 2; nnodes = 1
    elif name[0:3] == 'BAR':
        if len(name) == 3: nnodes = 2
        else: nnodes = int(name[4:])
        if nnodes == 2: eltno = 3
        elif nnodes == 3: eltno = 4
    elif name[0:3] == 'TRI':
        if len(name) == 3: nnodes = 3
        else: nnodes = int(name[4:])
        if nnodes == 3: eltno = 5
        elif nnodes == 6: eltno = 6
    elif name[0:4] == 'QUAD':
        if len(name) == 4: nnodes = 4
        else: nnodes = int(name[5:])
        if nnodes == 4: eltno = 7
        elif nnodes == 8: eltno = 8
        elif nnodes == 9: eltno = 9
    elif name[0:5] == 'TETRA':
        if len(name) == 5: nnodes = 4
        else: nnodes = int(name[6:])
        if nnodes == 4: eltno = 10
        elif nnodes == 10: eltno = 11
    elif name[0:4] == 'PYRA':
        if len(name) == 4: nnodes = 5
        else: nnodes = int(name[5:])
        if nnodes == 5: eltno = 12
        elif nnodes == 14: eltno = 13
    elif name[0:5] == 'PENTA':
        if len(name) == 5: nnodes = 6
        else: nnodes = int(name[6:])
        if nnodes == 6: eltno = 14
        elif nnodes == 15: eltno = 15
        elif nnodes == 18: eltno = 16
    elif name[0:4] == 'HEXA':
        if len(name) == 4: nnodes = 8
        else: nnodes = int(name[5:])
        if nnodes == 8: eltno = 17
        elif nnodes == 20: eltno = 18
        elif nnodes == 27: eltno = 19
    elif name == 'MIXED':
        print('Warning: eltName2EltNo: MIXED elements not supported.')
        eltno = 20; nnodes = -1
    elif name == 'NGON' or name == 'NGON_n': eltno = 22; nnodes = 1
    elif name == 'NFACE' or name == 'NFACE_n': eltno = 23; nnodes = 1
    return eltno, nnodes

# -- EltNo2EltName
# Convertit un numero CGNS d'elt en nom CGNS d'elt et son nombre de noeuds
# associes
def eltNo2EltName(eltno):
    name = 'UNKNOWN'; nnodes = 0
    if eltno == 2: name = 'NODE'; nnodes = 1
    elif eltno == 3: name = 'BAR'; nnodes = 2
    elif eltno == 4: name = 'BAR_3'; nnodes = 3
    elif eltno == 5: name = 'TRI'; nnodes = 3
    elif eltno == 6: name = 'TRI_6'; nnodes = 6
    elif eltno == 7: name = 'QUAD'; nnodes = 4
    elif eltno == 8: name = 'QUAD_8'; nnodes = 8
    elif eltno == 9: name = 'QUAD_9'; nnodes = 9
    elif eltno == 10: name = 'TETRA'; nnodes = 4
    elif eltno == 11: name = 'TETRA_10'; nnodes = 10
    elif eltno == 12: name = 'PYRA'; nnodes = 5
    elif eltno == 13: name = 'PYRA_14'; nnodes = 14
    elif eltno == 14: name = 'PENTA'; nnodes = 6
    elif eltno == 15: name = 'PENTA_15'; nnodes = 15
    elif eltno == 16: name = 'PENTA_18'; nnodes = 18
    elif eltno == 17: name = 'HEXA'; nnodes = 8
    elif eltno == 18: name = 'HEXA_20'; nnodes = 20
    elif eltno == 19: name = 'HEXA_27'; nnodes = 27
    elif eltno == 20:
        print('Warning: eltNo2EltName: MIXED elements not supported.')
        name = 'MIXED'; nnodes = -1
    elif eltno == 22: name = 'NGON'; nnodes = 1
    elif eltno == 23: name = 'NFACE'; nnodes = 1
    return name, nnodes

# Donne la dimension d'un element a partir de son no
def dimFromEltNo(eltno):
    if eltno == 2: return 0 # NODE
    elif eltno >= 3 and eltno <= 4: return 1 # BAR
    elif eltno >= 5 and eltno <= 6: return 2 # TRI
    elif eltno >= 7 and eltno <= 9: return 2 # QUAD
    else: return 3

# -- Convertit un PointRange (pyTree) en indices de fenetres (Converter)
def range2Window(r):
    if r[0,0] > r[0,1]: imin = r[0,1]; imax = r[0,0]
    else: imin = r[0,0]; imax = r[0,1]
    if r[1,0] > r[1,1]: jmin = r[1,1]; jmax = r[1,0]
    else: jmin = r[1,0]; jmax = r[1,1]
    if r.shape == (1,1): return [imin,imax,1,1,1,1]
    elif r.shape == (2,2): return [imin,imax,jmin,jmax,1,1]
    if r[2,0] > r[2,1]: kmin = r[2,1]; kmax = r[2,0]
    else: kmin = r[2,0]; kmax = r[2,1]
    return [imin, imax, jmin, jmax, kmin, kmax]

# -- Convertit une fenetre [imin,imax,jmin,jmax,kmin,kmax] en PointRange pyTree
def window2Range(win):
    r = numpy.empty((3,2), numpy.int32, order='Fortran')
    r[0,0] = win[0]; r[0,1] = win[1]
    r[1,0] = win[2]; r[1,1] = win[3]
    r[2,0] = win[4]; r[2,1] = win[5]
    return r

# -- PointList (face) to range (struct zones)
def pointList2Windows(PL, ni, nj, nk):
    # find face
    ind = PL[0]
    # find i,j
    # regrouper dans des vector[j]
    # pour chaque j, trouve imin, imax
    # Regroupe en j si meme imin, imax


# -- ClearList: supprime les [] d'une liste
def clearList(list):
  t = []
  for i in list:
    if i != []: t.append(i)
  return t

# -- Change a var name to a CGNS name
def getCGNSName(v):
    if v in name2CGNS: return name2CGNS[v]
    else: return v

#==============================================================================
# -- Conversion Arrays / Nodes --
#==============================================================================

# -- Array2PyTreeDim
# Retourne un numpy contenant les dimensions de la zone comme celui utilise
# dans le noeud zone de l'arbre.
# Les arrays structures sont tous nixnjxnk.
# Dans l'arbre python, il faut mettre les grilles 2D sous forme nixnj
# et les grilles 1D sous forme ni.
# ATTENTION: l'array doit etre en noeuds!
def array2PyTreeDim(a):
    if len(a) == 5: # structure
        ni = a[2]; nj = a[3]; nk = a[4]
        if ni == 1:
            if nj == 1:
                d = numpy.empty((1,3), numpy.int32, order='Fortran')
                d[0,0] = nk
                d[0,1] = nk-1
                d[0,2] = 0
            elif nk == 1:
                d = numpy.empty((1,3), numpy.int32, order='Fortran')
                d[0,0] = nj
                d[0,1] = nj-1
                d[0,2] = 0
            else:
                d = numpy.empty((2,3), numpy.int32, order='Fortran')
                d[0,0] = nj;   d[1,0] = nk
                d[0,1] = nj-1; d[1,1] = nk-1
                d[0,2] = 0;    d[1,2] = 0
        elif nj == 1:
            if nk == 1:
                d = numpy.empty((1,3), numpy.int32, order='Fortran')
                d[0,0] = ni
                d[0,1] = ni-1
                d[0,2] = 0
            else:
                d = numpy.empty((2,3), numpy.int32, order='Fortran')
                d[0,0] = ni;   d[1,0] = nk
                d[0,1] = ni-1; d[1,1] = nk-1
                d[0,2] = 0;    d[1,2] = 0
        elif nk == 1:
            d = numpy.empty((2,3), numpy.int32, order='Fortran')
            d[0,0] = ni;   d[1,0] = nj
            d[0,1] = ni-1; d[1,1] = nj-1
            d[0,2] = 0;    d[1,2] = 0
        else:
            d = numpy.empty((3,3), numpy.int32, order='Fortran')
            d[0,0] = ni;   d[1,0] = nj;   d[2,0] = nk
            d[0,1] = ni-1; d[1,1] = nj-1; d[2,1] = nk-1
            d[0,2] = 0;    d[1,2] = 0;    d[2,2] = 0
    elif len(a) == 4: # non structure
        if isinstance(a[2], list): # Array2
            if a[3] == 'NGON':
                nelts = a[2][3].size
                d = numpy.empty((1,3), numpy.int32, order='Fortran')
                d[0,0] = a[1][0].size; d[0,1] = nelts; d[0,2] = 0
            else:
                d = numpy.empty((1,3), numpy.int32, order='Fortran')
                d[0,0] = a[1][0].size; d[0,1] = a[2][0].shape[0]; d[0,2] = 0 
        else:   # Array1
            if a[3] == 'NGON':
                sizeFN = a[2][0,1]; nelts = a[2][0,2+sizeFN]
                d = numpy.empty((1,3), numpy.int32, order='Fortran')
                d[0,0] = a[1].shape[1]; d[0,1] = nelts; d[0,2] = 0
            else:
                d = numpy.empty((1,3), numpy.int32, order='Fortran')
                d[0,0] = a[1].shape[1]; d[0,1] = a[2].shape[1]; d[0,2] = 0
    else:
        d = None
    return d

# -- Cree un noeud de donnees (DataArray_t) a partir d'un array field nfld --
# cellDim = zone cell dimension
def createDataNode(name, array, nfld, cellDim=3):
  name = getCGNSName(name)
  if len(array) == 5: # structure
      ni = array[2]; nj = array[3]; nk = array[4]
      if isinstance(array[1], list): a = array[1][nfld] # array 2
      else: a = array[1][nfld,:] # array 1
      if cellDim == 3:
          a = numpy.reshape(a, (ni, nj, nk), order='Fortran')
      elif cellDim == 2:
          if ni == 1 or ni == 0:
              a = numpy.reshape(a, (nj, nk), order='Fortran')
          elif nj == 1 or nj == 0:
              a = numpy.reshape(a, (ni, nk), order='Fortran')
          else:
              a = numpy.reshape(a, (ni, nj), order='Fortran')
      else: # cellDim = 1
          if ni == 1 or ni == 0:
              if nj == 1 or nj == 0:
                  a = numpy.reshape(a, (nk), order='Fortran')
              elif nk == 1 or nk == 0:
                  a = numpy.reshape(a, (nj), order='Fortran')
          elif nj == 1 or nj == 0:
              if nk == 1 or nk == 0:
                  a = numpy.reshape(a, (ni), order='Fortran')
  else: # non structure
    if isinstance(array[1], list):
        a = array[1][nfld] # Array2
        a = numpy.reshape(a, (a.size), order='Fortran')
    else: a = array[1][nfld,:] # Array1
  node = [name, a, [], 'DataArray_t']
  return node

# -- Create a zone node from 2 arrays --
# name: nom de la zone
# array: champs en noeuds
# array2: champs en centres
# Les champs sont stockes dans __FlowSolutionNodes__, ...
def createZoneNode(name, array, array2=[],
                   GridCoordinates=__GridCoordinates__,
                   FlowSolutionNodes=__FlowSolutionNodes__,
                   FlowSolutionCenters=__FlowSolutionCenters__):
  zone = [name, array2PyTreeDim(array), [], 'Zone_t']
  info = zone[2]
  cellDim = zone[1].shape[0]

  # Zone type
  if len(array) == 5: type = 'Structured'
  else: type = 'Unstructured'
  v = numpy.fromstring(type, 'c')
  info.append(['ZoneType', v, [], 'ZoneType_t'])

  # Coordonnees
  px = KCore.isNamePresent(array, 'CoordinateX')
  if px == -1: px = KCore.isNamePresent(array, 'x')
  if px == -1: px = KCore.isNamePresent(array, 'X')
  py = KCore.isNamePresent(array, 'CoordinateY')
  if py == -1: py = KCore.isNamePresent(array, 'y')
  if py == -1: py = KCore.isNamePresent(array, 'Y')
  pz = KCore.isNamePresent(array, 'CoordinateZ')
  if pz == -1: pz = KCore.isNamePresent(array, 'z')
  if pz == -1: pz = KCore.isNamePresent(array, 'Z')
  if px != -1 or py != -1 or pz != -1:
    info.append([GridCoordinates, None, [], 'GridCoordinates_t'])
    if px != -1:
      node = createDataNode('CoordinateX', array, px, cellDim)
      info[1][2].append(node)
    if py != -1:
      node = createDataNode('CoordinateY', array, py, cellDim)
      info[1][2].append(node)
    if pz != -1:
      node = createDataNode('CoordinateZ', array, pz, cellDim)
      info[1][2].append(node)

  # Connectivite
  if len(array) == 4: # non structure
      etype,stype = eltName2EltNo(array[3])
      i = numpy.empty((2), numpy.int32); i[0] = etype; i[1] = 0
      if etype == 22: # Faces->Nodes and Elements->Faces connectivities (NGON array)
          if isinstance(array[2], list): # Array2
            setElementConnectivity2(zone, array)
          else: # Array1
            setElementConnectivity(zone, array)
      else:  # Elements -> Nodes connectivities
          info.append(['GridElements', i, [], 'Elements_t'])
          info2 = info[len(info)-1][2]
          i = numpy.empty((2), numpy.int32); i[0] = 1
          if isinstance(array[2], list): # Array2
             i[1] = array[2][0].shape[0]
             info2.append(['ElementRange', i, [], 'IndexRange_t']) 
             sizeConnectivity = array[2][0].size
             connect = array[2][0]
             connect = numpy.reshape(connect, (sizeConnectivity))
             info2.append(['ElementConnectivity', connect, [], 'DataArray_t'])  
          else: # Array1
             i[1] = array[2].shape[1]
             sizeConnectivity = array[2].size
             info2.append(['ElementRange', i, [], 'IndexRange_t'])
             connect = numpy.transpose(numpy.copy(array[2]))
             connect = numpy.reshape(connect, (sizeConnectivity))
             info2.append(['ElementConnectivity', connect, [], 'DataArray_t'])

  # field
  if isinstance(array[1], list): nvar = len(array[1])
  else: nvar = array[1].shape[0]
  createFlow = False
  for i in range(nvar):
    if i != px and i != py and i != pz: createFlow = True
  if createFlow:
    vars = array[0].split(',')
    info.append([FlowSolutionNodes, None, [], 'FlowSolution_t'])
    info = info[len(info)-1]
    for i in range(nvar):
      if i != px and i != py and i != pz:
        node = createDataNode(vars[i], array, i, cellDim)
        info[2].append(node)
  if array2 != []:
      if isinstance(array2[1], list): nvar = len(array2[1])
      else: nvar = array2[1].shape[0]
      vars = array2[0].split(',')
      info.append([FlowSolutionCenters, None, [], 'FlowSolution_t'])
      info = info[len(info)-1]
      v = numpy.fromstring('CellCenter', 'c')
      info[2].append(['GridLocation', v, [], 'GridLocation_t'])
      for i in range(nvar):
          node = createDataNode(vars[i], array2, i, cellDim)
          info[2].append(node)
  return zone

# -- convert a data node to an array (array1)
# IN: node: noeud de l'arbre a transformer en array
# IN: dim: comme issu de getZoneDim
# IN: connects: liste des noeuds Elements_t de la zone correspondante
# IN: loc: localisation de node si connue (-1: inconnu, 0: noeuds, 1: centres)
# Retourne [loc, array], loc='nodes', 'centers', 'unknown'
# array=un array derive de node si possible, sinon None
def convertDataNode2Array(node, dim, connects, loc=-1):
    gtype = dim[0]
    array = None
    ar = node[1]
    if ar is None: return ['unknown', None]
    if ar.dtype != numpy.float64: # mauvais type de noeud, copie
        ar2 = numpy.empty(ar.shape, numpy.float64)
        ar2[:] = ar
        ar = ar2

    if isinstance(gtype, numpy.ndarray): gtype = gtype.tostring()
    if gtype == 'Structured':
        ni = dim[1]; nj = dim[2]; nk = dim[3]
        ni1 = max(ni-1,1); nj1 = max(nj-1,1); nk1 = max(nk-1,1)
        if ni*nj*nk == ar.size: # champ en noeuds
            a = numpy.reshape(ar, (ni*nj*nk), order='Fortran')
            a = numpy.reshape(a, (1,ni*nj*nk))
            array = [node[0], a, ni, nj, nk]
            return ['nodes', array]
        elif ni1*nj1*nk1 == ar.size: # champ en centres
            a = numpy.reshape(ar, (ni1*nj1*nk1), order='Fortran')
            a = numpy.reshape(a, (1,ni1*nj1*nk1))
            array = [node[0], a, ni1, nj1, nk1]
            return ['centers', array]
        else:
            # Le tableau n'est pas coherent avec les noeuds ou les centres
            # de la zone
            size = ar.shape; lsize = len(size)
            if lsize == 1: ni = size[0]; nj = 1; nk = 1
            elif lsize == 2: ni = size[0]; nj = size[1]; nk = 1
            elif lsize == 3: ni = size[0]; nj = size[1]; nk = size[2]
            a = numpy.reshape(ar, (1,ni*nj*nk), order='Fortran')
            array = [node[0], a, ni, nj, nk]
            print("Warning: convertDataNode2Array: incoherency zone/array (%s)."%node[0])
            return ['unknown', array]
    elif gtype == 'Unstructured':
        cr, ettype = adaptConnect__(connects, dim)
        lconnects = len(connects)
        s = ar.size
        a = numpy.reshape(ar, (1, s), order='Fortran')
        locout = 'nodes'
        if dim[1] != dim[2]: # on peut decider
            if s == dim[2]: ettype += '*'; locout = 'centers'
            elif s != dim[1]:
                print("Warning: convertDataNode2Array: incoherency zone/array (%s)."%node[0])
        else: # force + no check
            if loc == 1: ettype += '*'; locout = 'centers'

        array = [node[0], a, cr, ettype]
        return [locout, array]

# Pour array2
def convertDataNode2Array2(node, dim, connects, loc=-1):
    gtype = dim[0]
    array = None
    ar = node[1]
    if ar is None: return ['unknown', None]
    if ar.dtype != numpy.float64: # mauvais type de noeud, copie
        ar2 = numpy.empty(ar.shape, numpy.float64)
        ar2[:] = ar
        ar = ar2

    if isinstance(gtype, numpy.ndarray): gtype = gtype.tostring()
    if gtype == 'Structured':
        ni = dim[1]; nj = dim[2]; nk = dim[3]
        ni1 = max(ni-1,1); nj1 = max(nj-1,1); nk1 = max(nk-1,1)
        if ni*nj*nk == ar.size: # champ en noeuds
            array = [node[0], [ar], ni, nj, nk]
            return ['nodes', array]
        elif ni1*nj1*nk1 == ar.size: # champ en centres
            array = [node[0], [ar], ni1, nj1, nk1]
            return ['centers', array]
        else:
            # Le tableau n'est pas coherent avec les noeuds ou les centres
            # de la zone
            size = ar.shape; lsize = len(size)
            if lsize == 1: ni = size[0]; nj = 1; nk = 1
            elif lsize == 2: ni = size[0]; nj = size[1]; nk = 1
            elif lsize == 3: ni = size[0]; nj = size[1]; nk = size[2]
            array = [node[0], [ar], ni, nj, nk]
            print("Warning: convertDataNode2Array: incoherency zone/array (%s)."%node[0])
            return ['unknown', array]

    # non structure
    iBE = -1; iBEMultiple = -1; iNGon = -1; iNFace = -1; i = 0
    for c in connects:
        ctype = c[1][0]
        if ctype == 22: iNGon = i
        elif ctype == 23: iNFace = i
        else:
            if iBE == -1: iBE = i # first
            else: iBEMultiple = 1
        i += 1
    if iNGon != -1 and iNFace != -1: # NGON
        eltType = 22
        connect1 = connects[iNGon]; connect2 = connects[iNFace]
    elif iBE != -1: # une connectivite BE existe -> on la prend
        connect = connects[iBE]; eltType = connect[1][0]
        if iBEMultiple == 1:
            print('Warning: convertDataNode2Array: different connectivities in a single zone is not possible. Only first connectivity is kept.')
    else:
        #raise ValueError("convertDataNode2Array: no valid connectivity found.")
        #print('Warning: convertDataNode2Array: no valid connectivity found (using NODE).')
        cr = numpy.empty((1,0), dtype=numpy.int32); ettype='NODE'

    ettype, stype = eltNo2EltName(eltType)

    if eltType == 22: # NGON
        info = connect1
        e = getNodeFromName(connect1, 'ElementConnectivity')
        f = getNodeFromName(connect2, 'ElementConnectivity')
        g = getNodeFromName(connect1, 'FaceIndex')
        if g is None:
            no = getNodeFromName1(connect1, 'ElementRange')[1]
            nfaces = no[1]-no[0]+1
            n = converter.adaptNGon2Index(e[1], nfaces)
            _createUniqueChild(connect1, 'FaceIndex', 'DataArray_t', value=n)
            g = getNodeFromName(connect1, 'FaceIndex')
        h = getNodeFromName(connect1, 'ElementIndex')
        if h is None: 
            no = getNodeFromName1(connect2, 'ElementRange')[1]
            nelts = no[1]-no[0]+1
            n = converter.adaptNFace2Index(f[1], nelts)
            _createUniqueChild(connect2, 'ElementIndex', 'DataArray_t', value=n)
            h = getNodeFromName(connect2, 'ElementIndex')
        cr = [e[1], f[1], g[1], h[1]]
    else:
        e = getNodeFromName(connect, 'ElementConnectivity')
        if e[1] is not None:
            size = e[1].size//stype
            b = e[1].reshape((size,stype)) # reshape needed
            cr = [b]
        else: cr = [None] # nodes
    locout = 'nodes'
    s = ar.size
    if dim[1] != dim[2]: # on peut decider
        if s == dim[2]: ettype += '*'; locout = 'centers'
        elif s != dim[1]:
            print("Warning: convertDataNode2Array: incoherency zone/array (%s)."%node[0])
    else: # force + no check
        if loc == 1: ettype += '*'; locout = 'centers'

    array = [node[0], [ar], cr, ettype]
    return [locout, array]

# -- convert a data node list to an array
# IN: nodes: liste de noeuds de l'arbre DataArray_t a transformer en array
# IN: dim: comme issu de getZoneDim
# IN: connects: liste des noeuds Elements_t de la zone correspondante
# (connectivites volumiques)
# IN: loc: localisation de nodes (-1: inconnu, 0: nodes, 1: centers)
# Retourne un array,
# array=un array derive de nodes si possible, sinon None
def convertDataNodes2Array(nodes, dim, connects, loc=-1):
    gtype = dim[0]
    nfld = len(nodes)
    if nfld == 0: return []
    ar = nodes[0][1]
    if ar is None: return []

    s = ar.size
    a = numpy.empty((nfld,s), numpy.float64)
    vars = nodes[0][0]
    for n in nodes[1:]: vars += ','+n[0]
    c = 0
    for n in nodes:
        if n[1] is None: a[c,:] = 0.
        elif n[1].dtype != numpy.float64: # mauvais type de noeud
            n2 = numpy.empty(n[1].shape, numpy.float64)
            n2[:] = n[1]
            n[1] = n2
            converter.cpyValueByField(a, n[1], s, c)
        else:
            converter.cpyValueByField(a, n[1], s, c)
        c+= 1

    if gtype == 'Structured':
        ni = dim[1]; nj = dim[2]; nk = dim[3]
        ni1 = max(ni-1,1); nj1 = max(nj-1,1); nk1 = max(nk-1,1)
        if ni*nj*nk == s:
            return [vars, a, ni, nj, nk]
        elif ni1*nj1*nk1 == s:
            return [vars, a, ni1, nj1, nk1]
        else:
            # Le tableau n'est pas coherent avec les noeuds ou les centres
            # de la zone
            print("Warning: convertDataNodes2Array: incoherency zone (%dx%dx%d)/array (%s,%d)."%(ni,nj,nk,nodes[0][0],s))
            size = ar.shape; lsize = len(size)
            if lsize == 1: ni = size[0]; nj = 1; nk = 1
            elif lsize == 2: ni = size[0]; nj = size[1]; nk = 1
            elif lsize == 3: ni = size[0]; nj = size[1]; nk = size[2]
            return [vars, a, ni, nj, nk]

    # unstructured
    cr, ettype = adaptConnect__(connects, dim)

    if dim[1] != dim[2]: # on peut decider
        if s == dim[2]: ettype += '*'
        elif s != dim[1]:
            print("Warning: convertDataNodes2Array: incoherency zone (%d,%d)/array (%s,%d)."%(dim[1],dim[2],nodes[0][0],s))
    else: # force + no check
        if loc == 1: ettype += '*'
    return [vars, a, cr, ettype]

# Pour array2
def convertDataNodes2Array2(nodes, dim, connects, loc=-1):
    gtype = dim[0]
    nfld = len(nodes)
    if nfld == 0: return []
    ar = nodes[0][1]
    if ar is None: return []
    s = ar.size

    vars = nodes[0][0]
    for n in nodes[1:]: vars += ','+n[0]
    field = []
    for n in nodes: field.append(n[1])

    if gtype == 'Structured':
        ni = dim[1]; nj = dim[2]; nk = dim[3]
        ni1 = max(ni-1,1); nj1 = max(nj-1,1); nk1 = max(nk-1,1)
        if ni*nj*nk == s:
            return [vars, field, ni, nj, nk]
        elif ni1*nj1*nk1 == s:
            return [vars, field, ni1, nj1, nk1]
        else:
            # Le tableau n'est pas coherent avec les noeuds ou les centres
            # de la zone
            size = ar.shape; lsize = len(size)
            if lsize == 1: ni = size[0]; nj = 1; nk = 1
            elif lsize == 2: ni = size[0]; nj = size[1]; nk = 1
            elif lsize == 3: ni = size[0]; nj = size[1]; nk = size[2]
            print("Warning: convertDataNodes2Array: incoherency zone/array.")
            return [vars, field, ni, nj, nk]

    # unstructured
    iBE = -1; iBEMultiple = -1; iNGon = -1; iNFace = -1; i = 0
    for c in connects:
        ctype = c[1][0]
        if ctype == 22: iNGon = i
        elif ctype == 23: iNFace = i
        else:
            if iBE == -1: iBE = i # first
            else: iBEMultiple = 1
        i += 1
    if iNGon != -1 and iNFace != -1: # NGON
        eltType = 22
        connect1 = connects[iNGon]; connect2 = connects[iNFace]
    elif iBE != -1: # une connectivite BE existe -> on la prend
        connect = connects[iBE]; eltType = connect[1][0]
        if iBEMultiple == 1:
            print('Warning: convertDataNode2Array: different connectivities in a single zone is not possible. Only first connectivity is kept.')
    else:
        #raise ValueError("convertDataNode2Array: no valid connectivity found.")
        #print('Warning: convertDataNode2Array: no valid connectivity found (using NODE).')
        cr = numpy.empty((1,0), dtype=numpy.int32); ettype='NODE'; eltType = 0

    ettype, stype = eltNo2EltName(eltType)

    if eltType == 22: # NGON
        info = connect1
        e = getNodeFromName(connect1, 'ElementConnectivity')
        f = getNodeFromName(connect2, 'ElementConnectivity')
        g = getNodeFromName(connect1, 'FaceIndex')
        if g is None:
            node = getNodeFromName1(connect1, 'ElementRange')[1]
            nfaces = node[1]-node[0]+1
            n = converter.adaptNGon2Index(e[1], nfaces)
            _createUniqueChild(connect1, 'FaceIndex', 'DataArray_t', value=n)
            g = getNodeFromName(connect1, 'FaceIndex')
        h = getNodeFromName(connect1, 'ElementIndex')
        if h is None: 
            node = getNodeFromName1(connect2, 'ElementRange')[1]
            nelts = node[1]-node[0]+1
            n = converter.adaptNFace2Index(f[1], nelts)
            _createUniqueChild(connect2, 'ElementIndex', 'DataArray_t', value=n)
            h = getNodeFromName(connect2, 'ElementIndex')
        cr = [e[1], f[1], g[1], h[1]]
    elif eltType != 0: # all elements except NODE
        e = getNodeFromName(connect, 'ElementConnectivity')
        if e[1] is not None:
            size = e[1].size//stype
            b = e[1].reshape((size,stype)) # reshape needed
            cr = [b]
        else: cr = [] # ou None?

    # tag *
    if dim[1] != dim[2]: # on peut decider
        if s == dim[2]: ettype += '*'
        elif s != dim[1]:
            print("Warning: convertDataNodes2Array: incoherency zone/array.")
    else: # force + no check
        if loc == 1: ettype += '*'
    return [vars, field, cr, ettype]

def _groupByFamily(t, familyChilds=None, unique=False):
    """Group all nodes of each base which have the same family name but not link to a family node."""
    if familyChilds is None: familyChilds=[]
    for b in getBases(t):
        # Name of Nodes Family_t
        familiesNodeNames = [fam[0] for fam in getNodesFromType(b,'Family_t')]
        # All node BC_t not family specified
        BCNotSpec = []
        for bc in getNodesFromType3(b, 'BC_t'):
            if not isValue(bc, 'FamilySpecified'):
                childname = getNodesFromName(bc, 'FamilyName')                
                if childname != []: BCNotSpec.append(bc)
        for bc in BCNotSpec:
            FamNameNode = getNodeFromType(bc, 'FamilyName_t')
            FamilyName= getValue(FamNameNode)
            solverBC = getNodeFromName(bc, '.Solver#BC')
            if not solverBC: solverBCChilds = []
            else: solverBCChilds = solverBC[2]

            # if a family node already exist
            if FamilyName in FamiliesNodeNames:
                setValue(bc, 'FamilySpecified')
            else:
                FamiliesNodeNames.append(FamilyName)
                FamBC=createNode('FamilyBC','FamilyBC_t',Internal.getValue(bc))                
                children = FamilyChilds[:]+solverBCChilds[:]# faut il faire une deepcopy ? dans XTree c est le cas
                if not unique: children.append(FamBC)
                else:
                    e = -1
                    for child in children:
                        if isName(child,getName(FamBC)) and isType(child,'Family_t'):
                            setValue(child,getValue(FamBC))
                            e = 0
                            break
                    if e==-1: children.append(FamBC)
                        
                #
                nodep=createNode(FamilyName,'Family_t', value=None,children=children,parent=b)
                setValue(bc,'FamilySpecified')
                if solverBC: _rmNode(bc,solverBC)
    return None

def groupByFamily(t,familyChilds=None, unique=False):
    tp = copyRef(t)
    _groupByFamily(tp)
    return tp

#==============================================================================
# -- Get zone dims --
#==============================================================================

# -- retourne les dimensions et le type de la grille
# pour un noeud de Zone donne --
# Retourne ['Structured', ni, nj, nk, celldim] pour une grille structuree
# Retourne ['Unstructured', np, ne, EltsName, celldim] pour une grille
# non structuree
# Avec: celldim: 1 pour les i-array ou les BAR-array, 2 pour les i,j-array
# ou les TRI, QUAD-array, 3 pour les autres
# EltsName: nom des elements (BAR, NODE, ..., NGON ou MULTI)
def getZoneDim(zone):
  """Return dimension information from a Zone_t node."""
  if len(zone) < 4: raise TypeError("getZoneDim: not a zone node.")
  if zone[3] != 'Zone_t': raise TypeError("getZoneDim: '%s' is not a zone node."%zone[0])
  dims = zone[1]
  # zone type
  info = zone[2]
  for i in info:
    # Noeud zone type
    if i[3] == 'ZoneType_t':
      cellDim = 3
      gtype = getValue(i)
      if gtype == 'Structured':
        if dims.shape[0] == 1:
          ni = int(dims[0,0]); nj = 1; nk = 1; cellDim = 1
        elif dims.shape[0] == 2:
          ni = int(dims[0,0]); nj = int(dims[1,0]); nk = 1; cellDim = 2
        else:
          ni = int(dims[0,0]); nj = int(dims[1,0]); nk = int(dims[2,0])
        return [gtype, ni, nj, nk, cellDim]
      elif gtype == 'Unstructured':
        np = int(dims[0,0]); ne = int(dims[0,1])
        c = getElementNodes(zone)
        lc = len(c)
        if lc == 0: return [gtype, np, ne, 'UNKNOWN', 3]
        if lc >= 2:
            # Y a t-il NGON et NFACE? Si oui, renvoie NGON meme si il
            # y a des BEs
            NGON = -1; NFACE = -1; NGONp = None
            for i in c:
                if i[1][0] == 22: NGON = 1; NGONp = i
                elif i[1][0] == 23: NFACE = 1
            if NGON == 1 and NFACE == 1:
                data = getNodeFromName1(NGONp, 'ElementConnectivity')
                datar = data[1]
                if datar is not None and datar.size > 0:
                    if datar[0] == 1: cellDim = 1
                    elif datar[0] == 2: cellDim = 2
                return [gtype, np, ne, 'NGON', cellDim]
            else: return [gtype, np, ne, 'MULTIPLE', 3]
        eltName,stype = eltNo2EltName(c[0][1][0])
        if lc == 2 and eltName != 'NGON' and eltName != 'NFACE':
            return [gtype, np, ne, 'MULTIPLE', 3]
        if eltName == 'NODE': cellDim = 0
        elif eltName == 'BAR': cellDim = 1
        elif eltName == 'TRI' or eltName == 'QUAD': cellDim = 2
        elif eltName == 'NGON':
            data = getNodeFromName1(c[0], 'ElementConnectivity')
            datar = data[1]
            if datar is not None and datar.size > 0:
                if datar[0] == 1: cellDim = 1
                elif datar[0] == 2: cellDim = 2
        elif eltName == 'NFACE':
            if len(c) > 1:
                data = getNodeFromName1(c[1], 'ElementConnectivity')
                datar = data[1]
                if datar is not None and datar.size > 0:
                    if datar[0] == 1: cellDim = 1
                    elif datar[0] == 2: cellDim = 2
                eltName = 'NGON'
        return [gtype, np, ne, eltName, cellDim]
      else:
        raise TypeError("getZoneDim: cannot find zone type for zone '%s'."%zone[0])
      break
  raise TypeError("getZoneDim: cannot find zone type for zone '%s'."%zone[0])

# -- getZoneType --
# Retourne 1 si la zone est structuree
# Retourne 2 si la zone est non structuree
# Retourne 0 sinon.
# Cette routine est plus rapide que getZoneDim
def getZoneType(zone):
    info = zone[2]
    for i in info:
        if i[3] == 'ZoneType_t':
            gtype = i[1][0] # OK si numpy ou string
            if gtype == 'S' or gtype == b'S': return 1
            else: return 2
    return 0

#==============================================================================
# -- Ghost cells management --
#===============================================================================

# -- Add ghost cells in a pyTree or a zone
# Returns a new zone with ghost cells
# IN: t: top tree
# IN: b: tree/basis/zone to modify (defined in t)
# IN: d: number of ghost cells to add
# IN: modified: container names or the list of variables to be modified
# If modified != all: only the provided field is filled with ghost cells
#                     -> flow field and zone dimensions are not consistent
#
# adaptBCs=0: BCs are not modified
#         =1: Physical BCs/BCMatch/BCOverlap ranges are modified
#             BCMatch ranges are not relevant geometrically (graphically)
# fillCorner: method to fill edges and corners
#         =1: edges and corners are filled (grid coordinates+flow solution)
#             according to the grid connectivity -> geometrically, the corners
#             and edges can be wrong
#         =0: neighbouring vectors are extrapolated to build edge cells, no
#             flow field is filled
def addGhostCells(t, b, d, adaptBCs=0, modified=[], fillCorner=1):
    """Add ghost cells to a pyTree.
    Usage: addGhostCells(t, b, d, modified, fillCorner)"""
    from . import GhostCells
    return GhostCells.addGhostCells(t, b, d, adaptBCs, modified, fillCorner)
def _addGhostCells(t, b, d, adaptBCs=0, modified=[], fillCorner=1):
    """Add ghost cells to a pyTree.
    Usage: addGhostCells(t, b, d, modified, fillCorner)"""
    from . import GhostCells
    GhostCells._addGhostCells(t, b, d, adaptBCs, modified, fillCorner)
    return None

# -- Remove ghost Cells on GridConnectivity and all fields of FlowSolution
# and FlowSolution#Centers
# Return a new zone
def rmGhostCells(t, b, d, adaptBCs=0, modified=[]):
    """Remove ghost cells to a pyTree.
    Usage: rmGhostCells(t, b, d, modified)"""
    from . import GhostCells
    return GhostCells.rmGhostCells(t, b, d, adaptBCs, modified)
def _rmGhostCells(t, b, d, adaptBCs=0, modified=[]):
    """Remove ghost cells to a pyTree.
    Usage: rmGhostCells(t, b, d, modified)"""
    from . import GhostCells
    return GhostCells._rmGhostCells(t, b, d, adaptBCs, modified)

# -- Remove rind cells
def rmRindCells(t, d, modified=[]):
    from . import GhostCells
    return GhostCells.rmRindCells(t, d, modified)

#==============================================================================
# -- Check --
#==============================================================================

# -- checkPyTree
def checkPyTree(t, level=-20):
    """Check pyTree conformity.
    Usage: checkPyTree(t, level)"""
    from . import Check
    return Check.checkPyTree(t, level)

# -- correctPyTree
def correctPyTree(t, level=-20):
    """Correct a pyTree.
    Usage: correctPyTree(t, level)"""
    from . import Check
    return Check.correctPyTree(t, level)

def _correctPyTree(t, level=-20):
    """Correct a pyTree.
    Usage: correctPyTree(t, level)"""
    from . import Check
    Check._correctPyTree(t, level)
    return None

# -- checkMultigrid
def checkMultigrid(t, level=1, nbMinCoarseB=5, nbMinCoarseW=3):
    """Check if zones, BC ranges and grid connectivity ranges are
    compatible with multigrid of given level.
    Usage: checkMultigrid(t, level, nbMinCoarseB, nbMinCoarseW)"""
    from . import Check
    return Check.checkMultigrid(t, level, nbMinCoarseB, nbMinCoarseW)

# -- checkSize
def checkSize(t, sizeMax=100000000):
    """Check if the number of points of a zone exceeds the sizeMax.
    Usage: checkSize(t,sizeMax)"""
    from . import Check
    return Check.checkSize(t,sizeMax)

# -- correctBCElementNodes
def _correctBCElementNodes(t):
    """Correct element nodes to tag them as BC."""
    from . import Check
    return Check._correctBCElementNodes(t)

#==============================================================================
# -- BC management --
#==============================================================================

# -- Add one layer to BCs. dir=1,2,3 (strctured zones)
def addOneLayer2BC(t, dir, N=1):
    a = copyRef(t)
    _addOneLayer2BC(a, dir, N)
    return a

def _addOneLayer2BC(a, dir, N=1):
    nodes = getZones(a)
    for z in nodes:
        # ZoneBC
        bnds = getNodesFromType2(z, 'BC_t')
        for bnd in bnds:
            wins = getNodesFromName1(bnd, 'PointRange')
            for w in wins:
                (parent, d) = getParentOfNode(bnd, w)
                win = range2Window(w[1])
                win[-1]+=N
                parent[2][d][1] = window2Range(win)

        # Connectivite match/nearmatch/overlap
        connect = getNodesFromType2(z, 'ZoneGridConnectivity_t')
        for cn in connect:
            wins = getNodesFromName2(cn, 'PointRange')
            transf =  getNodesFromName(cn, 'Transform')
            winsopp = getNodesFromName(cn, 'PointRangeDonor')
            for w in wins:
                (parent, d) = getParentOfNode(cn, w)
                win = range2Window(w[1])
                win[-1]+=N
                parent[2][d][1] = window2Range(win)
            nom = 0
            for w in winsopp: # passe seult pour les matchs/nearmatch
                if dir == 3 and len(transf[nom][1]) == 2:
                    n = numpy.empty((3), numpy.int32)
                    n[0] = transf[nom][1][0]; n[1] = transf[nom][1][1]; n[2] = 3
                    transf[nom][1] = n; diropp = 3
                else: diropp = abs(transf[nom][1][dir-1])
                (parent, d) = getParentOfNode(cn, w)
                win = range2Window(w[1])
                win[-1]+=N
                parent[2][d][1] = window2Range(win)
                nom += 1
    return None

# -- for each base, group all BCs of same BCType in a family named FamilyName
def groupBCByBCType(t, btype='BCWall', name='FamWall'):
    """Tag all BCs of given type with family named FamilyName."""
    a = copyRef(t)
    _groupBCByBCType(a, btype, name)
    return a

def _groupBCByBCType(t, btype='BCWall', name='FamWall'):
    for base in getBases(t):
        found = False
        for bc in getNodesFromType3(base,'BC_t'):
            if isValue(bc, btype):
                found = True
                setValue(bc, 'FamilySpecified')
                _createUniqueChild(bc, 'FamilyName', 'FamilyName_t', value=name)
        if found:
            bcvalue = createNode('FamilyBC', 'FamilyBC_t', value=btype)
            _createUniqueChild(base, name, 'Family_t', children=[bcvalue])
    return None

#==============================================================================
# -- ElementConnectivity management --
#==============================================================================

# -- adaptConnect
# Cette fonction choisit la connectivite volumique a passer en array
# avec la regle suivante:
# - si NGON et NFACE existe -> NGON
# - si existe une BE -> premiere BE
# IN: connects: liste de noeuds 'Elements_t' (connectivites volumiques)
# IN: dim: getZoneDim
# OUT: (numpy connect, type elts 'NODES' par ex)
def adaptConnect__(connects, dim):
    iBE = -1; iBEMultiple = -1; iNGon = -1; iNFace = -1; i = 0
    for c in connects:
        ctype = c[1][0]
        if ctype == 22: iNGon = i
        elif ctype == 23: iNFace = i
        else:
            if iBE == -1: iBE = i # first
            else: iBEMultiple = 1
        i += 1
    if iNGon != -1 and iNFace != -1: # NGON
        eltType = 22
        connect1 = connects[iNGon]; connect2 = connects[iNFace]
    elif iBE != -1: # une connectivite BE existe -> on la prend
        connect = connects[iBE]; eltType = connect[1][0]
        if iBEMultiple == 1:
            print('Warning: convertDataNode2Array: different connectivities in a single zone is not possible. Only first connectivity is kept.')
    else:
        #raise ValueError("convertDataNode2Array: no valid connectivity found.")
        #print('Warning: convertDataNode2Array: no valid connectivity found (using NODE).')
        return (numpy.empty((1,0), dtype=numpy.int32), 'NODE')

    np = dim[1]; ne = dim[2]
    ettype, stype = eltNo2EltName(eltType)

    if eltType == 22: # NGON
        info = connect1[2]
        c1 = None; nfaces = 0
        for i in info:
          if i[0] == 'ElementConnectivity': c1 = i[1]
          if i[0] == 'ElementRange': nfaces = i[1][1]-i[1][0]+1
        info = connect2[2]
        c2 = None; nelts = 0
        for i in info:
          if i[0] == 'ElementConnectivity': c2 = i[1]
          if i[0] == 'ElementRange': nelts = i[1][1]-i[1][0]+1
        st1 = c1.size; st2 = c2.size
        st = st1 + st2 + 4
        cr = numpy.empty((1, st), numpy.int32)
        converter.cpyConnectP2ConnectA(cr, c1, c2, stype, ne, nfaces, nelts)
    else: # Autres types
        info = connect[2]
        cr = None
        for i in info:
          if i[0] == 'ElementConnectivity':
            c = i[1]
            if c is None: # a NODE
                cr = numpy.empty((stype, 0), numpy.int32)
            else:
                nelts = c.size // stype # elts = ne sauf si MULTIPLE
                cr = numpy.empty((stype, nelts), numpy.int32)
                c2 = None
                converter.cpyConnectP2ConnectA(cr, c, c2, stype, nelts, -1, -1)
    return cr, ettype

# -- setElementConnectivity
# Remplit le noeud GridElements et ElementConnectivity de z a partir
# des donnees de array (non structure)
def setElementConnectivity(z, array):
  etype, stype = eltName2EltNo(array[3])
  GENodes = getElementNodes(z)
  i = numpy.empty((2), numpy.int32); i[0] = etype; i[1] = 0
  if GENodes == []:
      if etype != 22 and etype != 23: # Elements->Nodes connectivities
          z[2].append(['GridElements', i, [], 'Elements_t'])
          info = z[2][len(z[2])-1]
          i = numpy.empty((2), numpy.int32); i[0] = 1; i[1] = array[2].shape[1]
          info[2].append(['ElementRange', i, [], 'IndexRange_t'])
          connect = numpy.empty((array[2].size), numpy.int32)
          converter.cpyConnectA2ConnectP(array[2], connect,
                                         stype, array[2].shape[1])
          info[2].append(['ElementConnectivity', connect, [], 'DataArray_t'])
          _updateElementRange(z)
      else: # Faces->Nodes and Elements->Faces connectivities (NGON or NFACE)
          info = z[2]
          connect = numpy.transpose(numpy.copy(array[2]))
          # Creation du noeud NGON_n: connectivite Faces->Noeuds
          info.append(['NGonElements', i, [], 'Elements_t'])
          info2 = info[len(info)-1][2]
          # Size of ElementRange : nb de faces
          nfaces = array[2][0][0]
          i = numpy.empty((2), numpy.int32); i[0] = 1; i[1] = nfaces
          info2.append(['ElementRange', i, [], 'IndexRange_t'])
          # Tableau de connectivite Face/Noeuds
          sizeNGON = array[2][0][1]
          cFN = connect[2:sizeNGON+2]
          sizeFN = cFN.size
          cFN = numpy.reshape(cFN, (sizeFN))
          info2.append(['ElementConnectivity', cFN, [], 'DataArray_t'])
          _updateElementRange(z)
          # Creation du noeud NFACE_n: connectivite Elements->Faces
          etype,stype = eltName2EltNo('NFACE')
          i2 = numpy.empty((2), numpy.int32); i2[0] = etype; i2[1] = 0
          info.append(['NFaceElements', i2, [], 'Elements_t'])
          info2 = info[len(info)-1][2]
          # Size of ElementRange
          nelts = array[2][0][sizeNGON+2]
          i = numpy.empty((2), numpy.int32)
          i[0] = nfaces+1; i[1] = nfaces+nelts
          info2.append(['ElementRange', i, [], 'IndexRange_t'])
          # Tableau de connectivite Elements/Faces
          sizeNFACE = array[2][0][sizeNGON+3]
          cEF = connect[sizeNGON+4:sizeNGON+4+sizeNFACE]
          sizeEF = cEF.size; cEF = numpy.reshape(cEF, (sizeEF))
          info2.append(['ElementConnectivity', cEF, [], 'DataArray_t'])
          _updateElementRange(z)
  else: # la connectivite existe deja
      if etype != 22 and etype != 23: # Elements->Nodes connectivities
          GENodes[0][1] = i
          nodeE = getNodeFromName2(z, 'ElementRange')
          i = numpy.empty((2), numpy.int32); i[0] = 1; i[1] = array[2].shape[1]
          nodeE[1] = i
          connect = numpy.empty((array[2].size), numpy.int32)
          converter.cpyConnectA2ConnectP(array[2], connect,
                                         stype, array[2].shape[1])
          nodeE = getNodeFromName2(z, 'ElementConnectivity')
          nodeE[1] = connect
          _updateElementRange(z)
          for n in GENodes[1:]: # delete NFace ou NGon (si il existe)
              if n[1][0] == 23 or n[1][0] == 22:
                  (p, r) = getParentOfNode(z, n)
                  del p[2][r]
      else: # Faces->Nodes and Elements->Faces connectivities (NGON or NFACE)
          GEl = getNodesFromType1(z, 'Elements_t')
          for GE in GEl:
              if GE[1][1] == 0:
                  (p, r) = getParentOfNode(z, GE)
                  del p[2][r]
          info = z[2]
          # Creation du noeud NGON_n: connectivite Faces->Noeuds
          info.append(['NGonElements', i, [], 'Elements_t'])
          info2 = info[len(info)-1][2]
          # Size of ElementRange : nb de faces
          nfaces = array[2][0][0]
          i = numpy.empty((2), numpy.int32); i[0] = 1; i[1] = nfaces
          info2.append(['ElementRange', i, [], 'IndexRange_t'])
          # Tableau de connectivite Face/Noeuds
          sizeNGON = array[2][0][1]
          cFN = array[2][0][2:sizeNGON+2]
          sizeFN = cFN.size
          cFN = numpy.reshape(cFN, (sizeFN))
          info2.append(['ElementConnectivity', cFN, [], 'DataArray_t'])
          # Creation du noeud NFACE_n : connectivite Elements->Faces
          etype, stype = eltName2EltNo('NFACE')
          i2 = numpy.empty((2), numpy.int32); i2[0] = etype; i2[1] = 0
          info.append(['NFaceElements', i2, [], 'Elements_t'])
          info2 = info[len(info)-1][2]
          # Size of ElementRange
          nelts = array[2][0][sizeNGON+2]
          i = numpy.empty((2), numpy.int32)
          i[0] = nfaces+1; i[1] = nfaces+nelts
          info2.append(['ElementRange', i, [], 'IndexRange_t'])
          # Tableau de connectivite Elements/Faces
          sizeNFACE = array[2][0][sizeNGON+3]
          cEF = array[2][0][sizeNGON+4:sizeNGON+4+sizeNFACE]
          sizeEF = cEF.size; cEF = numpy.reshape(cEF, (sizeEF))
          info2.append(['ElementConnectivity', cEF, [], 'DataArray_t'])
          _updateElementRange(z)

# Pour array - api2
def setElementConnectivity2(z, array):
  etype, stype = eltName2EltNo(array[3])
  GENodes = getElementNodes(z)
  i = numpy.empty((2), numpy.int32); i[0] = etype; i[1] = 0
  if GENodes == []:
      if etype != 22 and etype != 23: # Elements->Nodes connectivities
          z[2].append(['GridElements', i, [], 'Elements_t'])
          info = z[2][len(z[2])-1]
          i = numpy.empty((2), numpy.int32); i[0] = 1; i[1] = array[2][0].shape[1]
          info[2].append(['ElementRange', i, [], 'IndexRange_t'])
          info[2].append(['ElementConnectivity', array[2][0], [], 'DataArray_t'])
          _updateElementRange(z)
      else: # Faces->Nodes and Elements->Faces connectivities (NGON or NFACE)
          info = z[2]
          # Creation du noeud NGON_n: connectivite Faces->Noeuds
          info.append(['NGonElements', i, [], 'Elements_t'])
          info2 = info[len(info)-1][2]
          # Size of ElementRange : nb de faces
          nfaces = array[2][2].size
          i = numpy.empty((2), numpy.int32); i[0] = 1; i[1] = nfaces
          info2.append(['ElementRange', i, [], 'IndexRange_t'])
          # Tableau de connectivite Face/Noeuds
          info2.append(['ElementConnectivity', array[2][0], [], 'DataArray_t'])
          # Tableau FaceIndex (PH)
          info2.append(['FaceIndex', array[2][2], [], 'DataArray_t'])
          _updateElementRange(z)
          # Creation du noeud NFACE_n: connectivite Elements->Faces
          etype,stype = eltName2EltNo('NFACE')
          i2 = numpy.empty((2), numpy.int32); i2[0] = etype; i2[1] = 0
          info.append(['NFaceElements', i2, [], 'Elements_t'])
          info2 = info[len(info)-1][2]
          # Size of ElementRange
          nelts = array[2][3].size
          i = numpy.empty((2), numpy.int32)
          i[0] = nfaces+1; i[1] = nfaces+nelts
          info2.append(['ElementRange', i, [], 'IndexRange_t'])
          # Tableau de connectivite Elements/Faces
          info2.append(['ElementConnectivity', array[2][1], [], 'DataArray_t'])
          # Tableau ElementIndex (PH)
          info2.append(['ElementIndex', array[2][3], [], 'DataArray_t'])
          _updateElementRange(z)
  else: # la connectivite existe deja
      if etype != 22 and etype != 23: # Elements->Nodes connectivities
          GENodes[0][1] = i
          nodeE = getNodeFromName2(z, 'ElementRange')
          i = numpy.empty((2), numpy.int32); i[0] = 1; i[1] = array[2][0].shape[1]
          nodeE[1] = i
          nodeE = getNodeFromName2(z, 'ElementConnectivity')
          nodeE[1] = array[2][0]
          _updateElementRange(z)
          for n in GENodes[1:]: # delete NFace ou NGon (si il existe)
              if n[1][0] == 23 or n[1][0] == 22:
                  (p, r) = getParentOfNode(z, n)
                  del p[2][r]
      else: # Faces->Nodes and Elements->Faces connectivities (NGON or NFACE)
          GEl = getNodesFromType1(z, 'Elements_t')
          for GE in GEl:
              if GE[1][1] == 0:
                  (p, r) = getParentOfNode(z, GE)
                  del p[2][r]
          info = z[2]
          # Creation du noeud NGON_n: connectivite Faces->Noeuds
          info.append(['NGonElements', i, [], 'Elements_t'])
          info2 = info[len(info)-1][2]
          # Size of ElementRange : nb de faces
          nfaces = array[2][2].size
          i = numpy.empty((2), numpy.int32); i[0] = 1; i[1] = nfaces
          info2.append(['ElementRange', i, [], 'IndexRange_t'])
          # Tableau de connectivite Face/Noeuds
          info2.append(['ElementConnectivity', array[2][0], [], 'DataArray_t'])
          # Tableau FaceIndex (PG)
          createUniqueChild(info2, 'FaceIndex', 'DataArray_t', array[2][2])
          # Creation du noeud NFACE_n : connectivite Elements->Faces
          etype, stype = eltName2EltNo('NFACE')
          i2 = numpy.empty((2), numpy.int32); i2[0] = etype; i2[1] = 0
          info.append(['NFaceElements', i2, [], 'Elements_t'])
          info2 = info[len(info)-1][2]
          # Size of ElementRange
          nelts = array[2][3].size
          i = numpy.empty((2), numpy.int32)
          i[0] = nfaces+1; i[1] = nfaces+nelts
          info2.append(['ElementRange', i, [], 'IndexRange_t'])
          # Tableau de connectivite Elements/Faces
          info2.append(['ElementConnectivity', array[2][1], [], 'DataArray_t'])
          # Tableau ElementIndex (PH)
          createUniqueChild(info2, 'ElementIndex', 'DataArray_t', array[2][3])
          _updateElementRange(z)

# -- Retourne une liste des noeuds Elements_t volumiques d'une zone
# Retourne [] si il n'y en a pas.
def getElementNodes(z):
    GEl = getNodesFromType1(z, 'Elements_t')
    #if GEl == []: return []
    #dimElt = [dimFromEltNo(i[1][0]) for i in GEl]
    #maxdim = max(dimElt)
    #c = 0; out = []
    #for GE in GEl:
    #    if dimElt[c] == maxdim: # connectivite volumique
    #        out.append(GE)
    #    c += 1
    out = []
    for GE in GEl:
        if GE[1][1] == 0: out.append(GE)
    return out

# -- Retourne une liste des noeuds Elements_t de boundary d'une zone
# Retourne [] si il n'y en a pas.
def getElementBoundaryNodes(z):
    GEl = getNodesFromType1(z, 'Elements_t')
    #if GEl == []: return []
    #dimElt = [dimFromEltNo(i[1][0]) for i in GEl]
    #maxdim = max(dimElt)
    #c = 0; out = []
    #for GE in GEl:
    #    if dimElt[c] < maxdim: # connectivite BC
    #        out.append(GE)
    #    c += 1
    out = []
    for GE in GEl:
        if GE[1][1] > 0: out.append(GE)
    return out

# -- Update la numerotation des ElementRanges dans l'ordre des connectivites
def _updateElementRange(z):
    GEl = getNodesFromType1(z, 'Elements_t')
    c = 0
    for GE in GEl:
        r = getNodeFromName1(GE, 'ElementRange')
        r[1] = numpy.copy(r[1]); a = r[1]
        size = a[1]-a[0]+1
        a[0] = c+1; c += size; a[1] = c
    return None

# -- getElementRange
# Retourne le range d'une connectivite d'une zone z
# IN: name: nom de la connectivite
# ou IN: type: type de la connectivite (premier rencontre)
# ou IN: number: no de la connectivite (first 0)
# Si une seule connectivite, la retourne
# Si pas trouve, retourne None
def getElementRange(z, name=None, type=None, number=None):
    dims = getZoneDim(z)
    if dims[0] == 'Unstructured':
        elts = getNodesFromType1(z, 'Elements_t')
        if len(elts) == 1:
            r = getNodeFromName1(elts[0], 'ElementRange')
            if r is not None: return [r[1][0], r[1][1]]
        if name is not None:
            for e in elts:
                if e[0] == name:
                    r = getNodeFromName1(e, 'ElementRange')
                    if r is not None: return [r[1][0], r[1][1]]
        elif number is not None:
            if (number < 0 or number > len(elts)-1): return None
            e = elts[number]
            r = getNodeFromName1(e, 'ElementRange')
            if r is not None: return [r[1][0], r[1][1]]
        elif type is not None:
            noType = eltName2EltNo(type)[0]
            for e in elts:
                if getValue(e)[0] == noType:
                    r = getNodeFromName1(e, 'ElementRange')
                    if r is not None: return [r[1][0], r[1][1]]
        return None
    else: # structure
        ne = dims[1]*dims[2]*dims[3]
        return [1, ne]

#==============================================================================
# -- Adapter --
#==============================================================================

# -- Adapte une connectivite avec ParentElement en connectivite NFACE
# remove = True: detruit la connectivite ParentElements
def adaptPE2NFace(t, remove=True):
    tp = copyRef(t)
    _adaptPE2NFace(tp, remove)
    return tp

def _adaptPE2NFace(t, remove=True):
    zones = getZones(t)
    for z in zones:
        parentElt = getNodeFromName2(z, 'ParentElements')
        if parentElt is not None:
            cFE = parentElt[1]
            cNFace, nelts = converter.adaptPE2NFace(cFE)
            p = createUniqueChild(z, 'NFaceElements', 'Elements_t',
                                  value=[23,0])
            createUniqueChild(p, 'ElementRange', 'IndexRange_t',
                              value=[1,nelts])
            createUniqueChild(p, 'ElementConnectivity',
                              'DataArray_t', value=cNFace)
            if remove: _rmNodesByName(z, 'ParentElements')
    return None

# -- Adapte une connectivite NFACE en connectivite ParentElement
# remove = True: detruit la connectivite NFace
# methodPE = 0 : methode geometrique pour generer le ParentElement (pour un maillage relativement regulier, sans cellules concaves).
# methodPE = 1 : methode topologique (pour un maillage quelconque).
def adaptNFace2PE(t, remove=True, methodPE=0):
    tp = copyRef(t)
    _adaptNFace2PE(tp, remove, methodPE)
    return tp

def _adaptNFace2PE(t, remove=True, methodPE=0):
    zones = getZones(t)
    for z in zones:
        npts = 0; nelts = 0; cNFace = None; NGON = None; noNFace = 0
        cNGon = None
        c = 0
        # Coordonnees
        XN = getNodeFromName2(z, 'CoordinateX')[1]
        YN = getNodeFromName2(z, 'CoordinateY')[1]
        ZN = getNodeFromName2(z, 'CoordinateZ')[1]

        # Connectivites NFace/NGon
        for e in z[2]:
            if e[3] == 'Elements_t':
                val = getValue(e)[0]
                if val == 22: # NGon
                    NGON = e
                    cNGon = getNodeFromName1(e, 'ElementConnectivity')[1]
                    node = getNodeFromName1(e, 'ElementRange')[1]
                    nfaces = node[1]-node[0]+1
                elif val == 23: # NFace
                    node = getNodeFromName1(e, 'ElementRange')[1]
                    nelts = node[1]-node[0]+1
                    noNFace = c
                    cNFace = getNodeFromName1(e, 'ElementConnectivity')[1]
            c += 1

        if (cNFace is not None and NGON is not None and cNGon is not None):
            cFE = converter.adaptNFace2PE(cNFace, cNGon, XN, YN, ZN, nelts, nfaces, methodPE)
            createUniqueChild(NGON, 'ParentElements', 'DataArray_t', value=cFE)
        if remove: del z[2][noNFace]
    return None

# -- Adapte NGON en FaceIndex
def _adaptNGon2Index(t):
    zones = getZones(t)
    for z in zones:
        for e in z[2]:
            if e[3] == 'Elements_t':
                val = getValue(e)[0]
                if val == 22: # NGon
                    NGON = e
                    cNGon = getNodeFromName1(e, 'ElementConnectivity')[1]
                    node = getNodeFromName1(e, 'ElementRange')[1]
                    nfaces = node[1]-node[0]+1
                    pos = converter.adaptNGon2Index(cNGon, nfaces)
                    createUniqueChild(NGON, 'FaceIndex', 'DataArray_t', value=pos)
    return None

# -- Adapte NFACE en ElementIndex
def _adaptNFace2Index(t):
    zones = getZones(t)
    for z in zones:
        for e in z[2]:
            if e[3] == 'Elements_t':
                val = getValue(e)[0]
                if val == 23: # NFace
                    NFACE = e
                    node = getNodeFromName1(e, 'ElementRange')[1]
                    nelts = node[1]-node[0]+1
                    cNFace = getNodeFromName1(e, 'ElementConnectivity')[1]
                    pos = converter.adaptNFace2Index(cNFace, nelts)
                    createUniqueChild(NFACE, 'ElementIndex', 'DataArray_t', value=pos)
    return None

# -- Adapte une condition aux limites definies par faces en conditions aux
# limites definies avec une BCC
def adaptBCFace2BCC(t, remove=True):
    tp = copyRef(t)
    _adaptBCFace2BCC(tp, remove)
    return tp

def _adaptBCFace2BCC(t, remove=True):
    zones = getZones(t)
    for z in zones:
        dims = getZoneDim(z)
        if dims[0] == 'Unstructured' and dims[3] != 'NGON':
            connect = getElementNodes(z)[0]
            connect = getNodeFromName(connect, 'ElementConnectivity')
            bnds = getNodesFromType2(z, 'BC_t')
            for b in bnds:
                a1 = getNodeFromName1(b, 'PointList')
                a2 = getNodeFromName1(b, 'GridLocation')
                if (a1 is not None and a2 is not None and getValue(a2) == 'FaceCenter'):

                    (BAR,TRI,QUAD) = converter.adaptBCFace2BCC(a1[1], connect[1], dims[3])
                    if BAR.size > 0:
                        n = createNode('Bnd_BAR', 'Elements_t', value=[2,BAR.size])
                        _addChild(z, n)
                        createChild(n, 'ElementRange', 'IndexRange_t', value=[1,BAR.size])
                        createChild(n, 'ElementConnectivity', 'DataArray_t', value=BAR)
                    if TRI.size > 0:
                        n = createNode('Bnd_TRI', 'Elements_t', value=[3,TRI.size])
                        _addChild(z, n)
                        createChild(n, 'ElementRange', 'IndexRange_t', value=[1,TRI.size])
                        createChild(n, 'ElementConnectivity', 'DataArray_t', value=TRI)
                    if QUAD.size > 0:
                        n = createNode('Bnd_QUAD', 'Elements_t', value=[4,QUAD.size])
                        _addChild(z, n)
                        createChild(n, 'ElementRange', 'IndexRange_t', value=[1,QUAD.size])
                        createChild(n, 'ElementConnectivity', 'DataArray_t', value=QUAD)
                    _updateElementRange(z)
                _createUniqueChild(b, 'ElementRange', 'IndexRange_t', value=[1,2])
    return None

# -- Adapte une condition aux limites definie par une BCC en
# condition aux limites definies par des faces
def adaptBCC2BCFace(t, remove=True):
    tp = copyRef(t)
    _adaptBCC2BCFace(tp, remove)
    return tp

# A finir
def _adaptBCC2BCFace(t, remove=True):
    zones = getZones(t)
    for z in zones:
        dims = getZoneDim(z)
        # get volumic connectivity
        cEV = getElementNodes(z)[0]
        # get boundary connectivities
        BCCs = getElementBoundaryNodes(z)
        # zone BCs
        zoneBCs = getNodesFromType1(z, 'ZoneBC_t')
        for bc in zoneBCs:
            # find BCC
            ra = getNodeFromType1(bc, 'IndexRange_t')
            bcc = None
            if ra is not None:
                for e in BCCs:
                    r = getNodeFromName1(e, 'ElementRange')
                    if (r is not None and r[0] == ra[0] and r[1] == ra[1]):
                        bcc = getNodeFromName1(e, 'ElementConnectivity'); break # found
            if bcc is not None:
                faces = converter.adaptBCC2BCFace(bcc[1], cEV[1], dims[3])
                #createChild(bc, 'PointList', 'IndexRange_t', value=faces)
            if remove: _rmNodesByName(z, bcc[0])
    return None

# -- adaptZoneBCEltRange2EltList
# Transforme des ElementRange des zoneBC en ElementList
# Uniquement pour les zones non-structurees
def adaptZoneBCEltRange2EltList(t):
    tp = copyRef(t)
    _adaptZoneBCEltRange2EltList(tp)
    return tp

def _adaptZoneBCEltRange2EltList(t):
    bcs = []
    zones = getZones(t)
    if zones == []: bcs = getNodesFromType(t, 'BC_t')
    else:
        for z in zones:
            if getZoneType(z) == 2: bcs += getNodesFromType(z, 'BC_t')

    for b in bcs:
        r = getNodeFromType1(b, 'IndexRange_t')
        v = getValue(r)
        val = numpy.arange(v[0][0],v[0][1]+1,dtype=numpy.int32)
        val = val.reshape((1,val.size))
        setValue(r, val)
        setType(r, 'IndexArray_t')
        setName(r, 'ElementList')
    return None

# input dict, t
def getElementRangeDict(t, d):
    zones = getZones(t)
    for z in zones:
        if z[0] not in d: d[z[0]] = {}
        elts = getNodesFromType1(z, 'Elements_t')
        for e in elts:
            n = getNodeFromName1(e, 'ElementRange')[1]
            if n is not None: d[z[0]][e[0]] = (n[0],n[1])
            else: d[z[0]][e[0]] = (1,1)

# input dict, index of a Elements of z
# trouve la connectivite de zname referencee par ind
def referencedElement(ind, zname, d):
    if zname not in d: return None
    p = d[zname]
    for k in p:
        if ind >= p[k][0] and ind <= p[k][1]: return k 
    return None

# -- fixNGon
# Fixing (temporaire)
# Si zone NGON + BE:
# remove=True: supprime les BE si il y a des NGONs dans une zone
# remove=False: met NGON en premier
# Si zone BE et breakBE=True
# Break la zone en plusieurs BE a simple element
# Si zone MIX et convertMIX=True
# Converti la zone MIXED en NGON
# Si addNFace=True, ajoute le NFace
def fixNGon(t, remove=False, breakBE=True, convertMIXED=True, addNFace=True):
    tp = copyRef(t)
    _fixNGon(tp, remove, breakBE, convertMIXED)
    return tp

def _fixNGon(t, remove=False, breakBE=True, convertMIXED=True, addNFace=True):
    zones = getZones(t)

    dictOfZTypes = {} # dictionnaire des types de zone (0: struct, 1: non struct)
    for z in zones:
        stype = getZoneType(z)
        dictOfZTypes[z[0]] = stype

    shift0 = {}
    getElementRangeDict(t, shift0)
    for z in zones:
        stype = dictOfZTypes[z[0]]
        if stype == 2: # Unstructured
            # Check NGON, NFACE, BE (first BE met)
            NGON = -1; NFACE = -1; BE = -1; NFACEORIG = True
            sons = z[2]; no = 0
            for s in sons:
                if s[3] == 'Elements_t':
                    ztype = s[1][0]
                    if ztype == 22: NGON = no
                    elif ztype == 23: NFACE = no
                    elif BE == -1: BE = no
                no += 1

            # si NFACE, abs face index
            if NFACE != -1:
                c = getNodeFromName1(sons[NFACE], 'ElementConnectivity')
                if c[1] is not None: c[1] = numpy.absolute(c[1])

            # Si NGON, ajoute NFace si manquant
            if addNFace == True and NFACE == -1 and NGON >= 0:
                # Essaie de transformer la connectivite ParentElements en NFACE
                ngon = sons[NGON]
                parentElt = getNodeFromName1(ngon, 'ParentElements')
                if parentElt is not None: # parent element est present
                    cFE = parentElt[1]
                    sh = cFE.shape
                    if len(sh) == 1: # Bug elsA
                        cFE = cFE.reshape((sh[0]/2,2), order='F'); parentElt[1] = cFE
                    cNFace, nelts = converter.adaptPE2NFace(cFE)
                    p = createUniqueChild(z, 'NFaceElements', 'Elements_t',
                                          value=[23,0])
                    createUniqueChild(p, 'ElementRange', 'IndexRange_t',
                                      value=[1,nelts])
                    createUniqueChild(p, 'ElementConnectivity',
                                      'DataArray_t', value=cNFace)
                    NFACE = len(sons)-1; NFACEORIG = False
                    if remove: _rmNode(z, parentElt)
                else:
                    print('Warning: cannot create NFACE. ParentElements node is not present.')
                _updateElementRange(z)

            # Reorder: NGON en premier, NFACE en deuxieme si BE
            if BE > 0 and NGON >= 0:
                ngon = sons[NGON]
                
                # Supprime ou reordonne les connectivites BE
                if remove: # supprime les BE
                    connects = getElementNodes(z)
                    for c in connects:
                        if c[1][0] != 22 and c[1][0] != 23:
                            _rmNode(z, c)
                    connects = getElementBoundaryNodes(z)
                    for c in connects: _rmNode(z, c)
                else: # swap
                    if NGON != -1 and NGON > BE: # swap BE and NGON
                        # swap un BE et NGON
                        temp = sons[BE]
                        sons[BE] = sons[NGON]; sons[NGON] = temp
                    if NFACE != -1 and NFACE != BE+1:
                        temp = sons[BE+1]
                        sons[BE+1] = sons[NFACE]; sons[NFACE] = temp

                _updateElementRange(z)

            # Break connectivity sur les zones BE multiple (mais pas NGON)
            if NGON == -1 and breakBE:
                connects = getElementNodes(z)
                if len(connects) > 1: # multiple volume connectivity
                    from . import PyTree; import Generator.PyTree as G
                    zones = PyTree.breakConnectivity(z)
                    zones = G.close(zones)
                    (p,c) = getParentOfNode(t, z)
                    if p is not None: 
                        p[2][c] = zones[0]; p[2] += zones[1:]
                        for z in zones: dictOfZTypes[z[0]] = 2
            # Convert MIXED to NGON
            if NGON == -1 and convertMIXED:
                connects = getElementNodes(z)
                if connects != []:
                    if connects[0][1][0] == 20:
                        from . import PyTree
                        PyTree._convertArray2NGon(z) # ineffective (2.2)
                                                 
    # Remet les BCs d'aplomb (en fonction de la renumerotation shift0->shift1)
    shift1 = {}
    getElementRangeDict(t, shift1)
    for z in zones:
        stype = dictOfZTypes[z[0]]
        if stype == 2: # Unstructured
            zoneBC = getNodeFromType1(z, 'ZoneBC_t')
            if zoneBC is not None: BCs = getNodesFromType1(zoneBC, 'BC_t')
            else: BCs = []
            for b in BCs:
               loc = getNodeFromName1(b, 'GridLocation')
               if loc is None:
                  createChild(b, 'GridLocation', 'GridLocation_t', 'FaceCenter')
               pl = getNodesFromType1(b, 'IndexArray_t')
               for p in pl:
                    if p[0] == 'ElementList': p[0] = 'PointList' # forced
                    if p[0] == 'PointList':
                        pln = p[1]
                        ind = pln.ravel('k')[0]
                        ref = referencedElement(ind, z[0], shift0)
                        if ref is None: 
                            print('Warning: Cannot find', ind, b[0], z[0])
                        else:
                            shiftn = shift1[z[0]][ref][0]-shift0[z[0]][ref][0]
                            pln[:] += shiftn
                            p[1] = pln.reshape((1, pln.size))
            zoneGC = getNodeFromType1(z, 'ZoneGridConnectivity_t')
            if zoneGC is not None:
                joins = getNodesFromType2(z, 'GridConnectivity_t')
            else: joins = []
            for b in joins:
                loc = getNodeFromName1(b, 'GridLocation')
                if loc is None:
                    createChild(b, 'GridLocation', 'GridLocation_t', 'FaceCenter')
                pl = getNodesFromType1(b, 'IndexArray_t')
                for p in pl:
                    if p[0] == 'ElementList': p[0] = 'PointList' # forced
                    if p[0] == 'PointList':
                        pln = p[1]
                        ind = pln.ravel('k')[0]
                        ref = referencedElement(ind, z[0], shift0)
                        if ref is None: 
                            print('Warning: Cannot find', ind, b[0], z[0])
                        else:
                            shiftn = shift1[z[0]][ref][0]-shift0[z[0]][ref][0]
                            pln[:] += shiftn
                            p[1] = pln.reshape((1, pln.size))
                    if p[0] == 'PointListDonor':
                        zdonor = getValue(b)
                        pln = p[1]
                        ind = pln.ravel('k')[0]
                        ref = referencedElement(ind, zdonor, shift0)
                        if ref is None:
                            print('Warning: cannot find', ind, b[0], z[0])
                        else:
                            shiftn = shift1[zdonor][ref][0]-shift0[zdonor][ref][0]
                            pln[:] += shiftn
                            p[1] = pln.reshape((1, pln.size))

        else: # Structured Zone
            # > For hybrid join we need to offset PointListDonor
            zoneGC = getNodeFromType1(z, 'ZoneGridConnectivity_t')
            if zoneGC is not None:
                joins = getNodesFromType2(z, 'GridConnectivity1to1_t')
            else: joins = []
            for b in joins:
                pld = getNodeFromName1(b, 'PointListDonor')
                if pld is not None:
                    zdonorname = getValue(b)
                    if zdonorname is None:
                       zdonorname = getNodeFromName(b, 'ZoneDonorName')
                       zdonorname = getValue(zdonorname)
                    dtype = dictOfZTypes.get(zdonorname, 0)
                    if dtype == 2: # hybrid join
                        pln = pld[1]
                        ind = pln.ravel('k')[0]
                        ref = referencedElement(ind, zdonorname, shift0)
                        if ref is None:
                            print('Warning: Cannot find', ind, b[0], z[0])
                        else:
                            shiftn = shift1[zdonorname][ref][0]-shift0[zdonorname][ref][0]
                            pln[:] += shiftn
                            pld[1] = pln.reshape((1, pln.size))
    return None

#==============================================================================
# Remet des elements a un NGon
#==============================================================================
# methodPE = 0 : methode geometrique pour generer le ParentElement (pour un maillage relativement regulier, sans cellules concaves).
# methodPE = 1 : methode topologique (pour un maillage quelconque).
def _unfixNGon(t, methodPE=0):
    from . import PyTree; import Transform.PyTree as T
    zones = getZones(t)
    for z in zones:
        _adaptNFace2PE(z, remove=False, methodPE=0)
        B = T.breakElements(z)
        for b in B:
            PyTree._mergeConnectivity(z, b, boundary=0)
    return None

#==============================================================================
# Create :elsAHybrid node, add ParentElements
# method=0: avec tri TRI,QUAD+interior, exterior
# method=1: avec tri interior, exterior
#==============================================================================
def createElsaHybrid(t, method=0, axe2D=0, methodPE=0):
    from . import elsAProfile
    tp = elsAProfile.createElsaHybrid(t, method, axe2D, methodPE)
    return tp

def _createElsaHybrid(t, method=0, axe2D=0, methodPE=0):
    from . import elsAProfile
    elsAProfile._createElsaHybrid(t, method, axe2D, methodPE)
    return None

#==============================================================================
# -- Purely internal (undocumented) --
#==============================================================================

# -- Convert an array with (i,j,k) numerotation to an array with 1D-index
# numerotation
# IN: array with (i,j,k) numerotation
# OUT: array with 1D-index numerotation
def convertIJKArray21DArray(*thetuple):
  if len(thetuple) == 4:
    (a,im,jm,km) = thetuple
    size = a.shape[1]
    b = numpy.empty((1,size), numpy.int32)
    for l in range(size):
      i = a[0][l]-1; j = a[1][l]-1; k = a[2][l]-1
      ind = adrNode1__(i,j,k,im,jm,km,0)
      b[0][l] = ind
    return b
  elif len(thetuple) == 3:
    (a,im,jm) = thetuple
    size = a.shape[1]
    b = numpy.empty((1,size), numpy.int32)
    for l in range(size):
      b[0][l] = a[0][l]-1; j = a[1][l]-1
      ind = adrNode2__(i,j,im,jm,0)
    return b
  elif len(thetuple) == 2:
    (a,im) = thetuple
    size = a.shape[1]
    b = numpy.empty((1,size), numpy.int32)
    for l in range(size):
      b[0][l] = a[l]-1
    return b
  else:
    raise TypeError("convertIJKArray2IndexArray: bad number of arguments.")
    return 0

# -- Return 1D-indice corresponding to (i,j) in 2D and (i,j,k) in 3D
def adrNode1__(i,j,k,im,jm,km,d):
    return i+d + (j+d)*(im+2*d)+(k+d)*(im+2*d)*(jm+2*d)
def adrNode2__(i,j,im,jm,d):
    return i+d + (j+d)*(im+2*d)

# -- Gather indices as structured patch
def gatherInStructPatch2D__(listIndices, indirWin, niw, njw, dirf, niZ, njZ, nkZ):
    restart = 1; res = []
    noBlk = -1
    while restart == 1:
        # Recherche du premier pt identifie
        noind = 0; indstart=-1; indend=-1
        for indH in listIndices:
            if indH != -1: indstart = noind; noBlk = indirWin[indH-1]; break
            noind += 1
        if indstart == -1: return []

        jstart = indstart//niw; istart = indstart-jstart*niw
        iend = niw-1; jend = njw-1;
        for i in range(istart+1,niw):
            ind = i+jstart*niw
            indH = listIndices[ind]
            if indH == -1 or indirWin[indH-1] != noBlk:break
            else: iend=i; jend=jstart; indend=ind
        for j in range(jstart+1,njw):
            ok = 1
            for i in range(istart,niw):
                ind = i+j*niw
                indH = listIndices[ind]
                if indH == -1 or indirWin[indH-1] != noBlk:
                    if i == istart: ok = 0; break
                    else: iend=min(i-1,iend); jend=j
                else: indend=ind
            if ok == 0: break
            else: jend=j
        # update listIndices
        for j in range(jstart,jend+1):
            for i in range(istart,iend+1):
                ind = i+j*niw
                listIndices[ind] = -1;
        res.append([istart,iend,jstart,jend,noBlk])
        if max(listIndices) == -1: restart=0

    out = []; noblks = []
    for r in res:
        istart = r[0]; iend = r[1]; jstart = r[2]; jend = r[3]
        if dirf==1: range0 = [1,1,istart+1,iend+1,jstart+1,jend+1]
        elif dirf==2: range0 = [niZ,niZ,istart+1,iend+1,jstart+1,jend+1]
        elif dirf==3: range0 = [istart+1,iend+1,1,1,jstart+1,jend+1]
        elif dirf==4: range0 = [istart+1,iend+1,njZ,njZ,jstart+1,jend+1]
        elif dirf==5: range0 = [istart+1,iend+1,jstart+1,jend+1,1,1]
        else: range0 = [istart+1,iend+1,jstart+1,jend+1,nkZ,nkZ]
        out.append(range0); noblks.append(r[4])
    return out, noblks

def getRotationAngleValueInDegrees(RotationAngleNode):
    if RotationAngleNode is None: return None
    angleUnitNode = getNodeFromType1(RotationAngleNode, 'DimensionalUnits_t')
    if angleUnitNode is not None:
        angleUnit = getValue(angleUnitNode)[4]
        if angleUnit == 'Radian': alpha = __RAD2DEG__
        elif angleUnit == 'Degree': alpha = 1.
        else:
            print('Warning: getRotationAngleValueInDegrees: unknown angle unit=%s.'%angleUnit)
            print('Warning: getRotationAngleInDegrees: AngleUnits must be Radian or Degree. Assuming Degree')
            alpha = 1.
    else: 
        alpha=1.#__RAD2DEG__
        print('Warning: getRotationAngleInDegrees: no angle units defined in RotationAngle node: assuming angle unit is Degree.')

    anglev = getValue(RotationAngleNode)
    if anglev.shape[0] != 3: raise ValueError("getRotationAngleInDegrees: RotationAngle value must be a numpy of size 3.")
    anglev = [v*alpha for v in anglev]
    return anglev

#=============================================================================
# Get the periodic connectivity information (rotation or translation)
# IN: gcnode: GridConnectivity node
# OUT : [xc,yc,zc,axisX,axisY,axisZ,angleDeg],[translX,translY,translZ]
#       avec angleDeg en degres ! 
#=============================================================================
def getPeriodicInfo__(gcnode):
    cont = getNodeFromName1(gcnode, 'GridConnectivityProperty')
    if cont is None: return [],[]
    per = getNodeFromName2(gcnode, 'Periodic')
    if per is None: return [],[]
    translVect = []; rotationData = []
    transl = getNodeFromName2(cont, 'Translation')
    tx = 0.; ty = 0.; tz = 0. # translation
    isTranslated = 0; isRotated = 0
    if transl is not None:
        transl = getValue(transl)
        tx = transl[0]; ty = transl[1]; tz = transl[2]
        if tx != 0. or ty != 0. or tz != 0.: translVect=[tx,ty,tz]

    angle=0.; xc=0.; yc=0.; zc=0.; vx=0.; vy=0.; vz=0. # rotation
    center = getNodeFromName2(cont, 'RotationCenter')
    anglev = getNodeFromName2(cont, 'RotationAngle')
    if center is not None and anglev is not None:
        (xc,yc,zc) = getValue(center)
        anglev = getRotationAngleValueInDegrees(anglev)
        if   anglev[0] != 0: angle = anglev[0]; vx = 1.; isRotated = 1
        elif anglev[1] != 0: angle = anglev[1]; vy = 1.; isRotated = 1
        elif anglev[2] != 0: angle = anglev[2]; vz = 1.; isRotated = 1
        if angle != 0.: rotationData = [xc,yc,zc,vx, vy, vz,angle]

    return rotationData, translVect
#=================================================================================
# Merge BCDataSets: special bcdataset (etc/FFD) are prefered, other are destroyed
#=================================================================================
def _mergeBCDataSets__(z, bcNode):
    dataSets = getNodesFromType1(bcNode, 'BCDataSet_t')
    if len(dataSets)==1: return None

    # parmi les BCDataSet on privilegie les specifiques etc/FFD dans cet ordre 
    path = '%s#ZoneBC/%s'%(__FlowSolutionCenters__, bcNode[0])
    dataSetByPath = getNodeFromPath(z, path)
    if dataSetByPath is not None: dataSetByPathName=dataSetByPath[0]
    else: dataSetByPathName = None
    dataSetLocs = []
    nod = -1; no = 0
    for d in dataSets:
        loc = getNodeFromName1(d,'GridLocation')
        if loc is not None:
            if getValue(loc)=='FaceCenter': dataSetLocs.append(1)
            else : dataSetLocs.append(0)
        else: dataSetLocs.append(0)
        dname = d[0]
        if dname == 'BCDataSet#EndOfRun': nod = no; break
        elif dname==dataSetByPathName and dataSetByPathName is not None: nod = no; break
        elif dname=='FFD72SurfaceSolution': nod = no; break
        no+=1
    if nod == -1: nod = 0# other cases : first bc data set is kept

    no = 0; locd = dataSetLocs[nod]
    bcdatas = getNodeFromType2(dataSets[nod], 'DataArray_t')
    cont, c = getParentOfNode(z, bcdatas)
    for d in dataSets:
        if no != nod:
            if dataSetLocs[nod]==locd: 
                parent,nop = getParentOfNode(z,d)
                datas = getNodesFromType2(d, 'DataArray_t')
                for data0 in datas:
                    _createUniqueChild(cont, data0[0], 'DataArray_t', value=data0[1])
                del parent[2][nop]
            else:
                print('Warning: BCDataSet location of %s is different from %s. Not merged.'%(dataSets[no], dataSets[nod]))
        no += 1
    return None

#==============================================================================
# Acces generalise au BCDataSet
# Soit dans la zoneBC: ZoneBC/BC/BCName/BCDataSet/BCNeumann/fields
#                   ou ZoneBC/BC/BCName/BCDataSet#EndOfRun/BCNeumann/fields
# Soit dans FlowSolution#Centers#ZoneBC/BCName/Neumann/fields
# Soit dans ZoneBC/FFD72SurfaceSolution/NeumannData/fields
# IN: bc: BC_t node
# IN: z: zone node
# OUT: data node list (liste de DataArray_t)
#==============================================================================
def getBCDataSet(z, bcNode, withLoc=False):
    datas = []; ploc = 'Vertex'
    # Try from BCDataSet#EndOfRun (etc)
    dataSet = getNodeFromName1(bcNode, 'BCDataSet#EndOfRun')
    if dataSet is not None:
        #print('found new etc style dataSet')
        datas = getNodesFromType2(dataSet, 'DataArray_t')
        if withLoc:
            l = 'Vertex'
            l = getNodeFromType1(dataSet, 'GridLocation_t')
            if l is not None: ploc = getValue(l)
            return datas,ploc
        else: return datas

    # Try from old style etc
    path = '%s#ZoneBC/%s'%(__FlowSolutionCenters__, bcNode[0])
    node = getNodeFromPath(z, path)
    if node is not None:
        #print('found old style etc Flow#ZoneBC')
        datas = getNodesFromType2(node, 'DataArray_t')
        if withLoc:
            l = 'Vertex'
            l = getNodeFromType1(node, 'GridLocation_t')
            if l is not None: ploc = getValue(l)
            return datas,ploc
        else: return datas

    # Try from FFD extraction
    node = getNodeFromName(bcNode, 'FFD72SurfaceSolution')
    if node is not None:
        #print('found FFD72')
        datas = getNodesFromType2(node, 'DataArray_t')
        if withLoc:
            l = 'Vertex'
            l = getNodeFromType1(node, 'GridLocation_t')
            if l is not None: ploc = getValue(l)
            return datas,ploc
        else: return datas

    # Try from other BCDataSet
    dataSet = getNodeFromType1(bcNode, 'BCDataSet_t')
    if dataSet is not None:
        #print('Found dataSet')
        datas = getNodesFromType2(dataSet, 'DataArray_t')
        if withLoc:
            l = 'Vertex'
            l = getNodeFromType1(dataSet, 'GridLocation_t')
            if l is not None: ploc = getValue(l)
            return datas,ploc
        else: return datas
    if withLoc: return None
    else: return datas
 
#==============================================================================
# Retourne une liste des faces de la BC (node)
# si la grille est structuree, retourne un indicage de faces
# IN: z: zone node
# IN: bcNode: BC node
# OUT: facelist node
#==============================================================================
def getBCFaceNode(z, bcNode):
    dims = getZoneDim(z)
    if dims[0] == 'Unstructured': return getNodeFromName1(bcNode, 'PointList')

    r = getNodeFromName1(bcNode, 'PointRange') # structure maintenant
    ni = dims[1]; nj = dims[2]; nk = dims[3]
    wins = range2Window(r[1])
    listIndices = converter.range2PointList(wins[0], wins[1], wins[2], wins[3], wins[4], wins[5], 
                                            ni, nj, nk)
    return createNode(__FACELIST__, 'DataArray_t', value=listIndices)
    
# Return a container per variable available in a BCDataSet 
# FlowSolution and BCDataSet must be at the same location Vertex/Centers 
# return list of type: [[numpyFS, 'CellCenter', [[BCRange1, numpyBCDS1],[BCRange2, numpyBCDS2]]]
def  getBCDataSetContainers(name, z):
    containers = []
    dims = getZoneDim(z)
    if dims[0] == 'Unstructured': 
        print('Internal: getBCDataSetContainers not yet implemented for unstructured zones.')
        return None

    if name == __GridCoordinates__: return None

    # Extract BCDataSet info: loc, bc range, fields
    allBCDatas=[]; allRanges=[]; allLocDS=[]
    for bcNode in getNodesFromType3(z, "BC_t"):
        res = getBCDataSet(z, bcNode, withLoc=True) # numpy of fields at bcs + loc
        if res is not None:
            datas,locds = res
            w = getNodeFromName2(bcNode, 'PointRange')
            w = range2Window(w[1]) # imin,imax,jmin,jmax,kmin,kmax
            allBCDatas.append(datas)# all the fields
            allRanges.append(w)
            allLocDS.append(locds)
    if allBCDatas == []: return None
    nbcdataset = len(allRanges)

    # Case 1 : get all the fields at a given location
    if (name  == __FlowSolutionNodes__) or (name == __FlowSolutionCenters__):
        loc = 'Vertex'
        f = getNodeFromName2(z, name)
        if f is not None:
            # location of flow solution
            l = getNodeFromType1(f, 'GridLocation_t')
            if l is not None: loc = getValue(l)
            else:
                if name == __FlowSolutionNodes__: loc = 'Vertex'
                else: loc = 'CellCenter'

            if loc == 'Vertex':  matchingLoc = 'Vertex'
            else: matchingLoc = 'FaceCenter'   

            for f0 in f[2]: # for all the fields stored in flow solution
                if f0[3] == 'DataArray_t':
                    fname = f0[0] # variable name to be extracted from bcdataset
                    dataSetL = [] # [[BCRange1,fieldBC],[BCRange2,fieldBC2],[...],...]
                    for nobcdata in range(nbcdataset):
                        if allLocDS[nobcdata]==matchingLoc: # OK: locations are consistent
                            datas = getNodesFromName1(allBCDatas[nobcdata],fname) 
                            if datas != []:
                                dataSetL+=[[allRanges[nobcdata], datas[0][1]]]
                    if dataSetL!=[]:# variable exists both in flow solution and in bcdataset
                        dataFS = [fname, loc, dataSetL]
                        containers.append(dataFS)
    else: 
        if not isinstance(name, list): name = [name]
        for v in name:
            varname = v.split(':')
            loc = 'Vertex'; container = __FlowSolutionNodes__
            if len(varname) == 2: # center value
                if varname[0] == 'centers': loc = 'CellCenter'; container = __FlowSolutionCenters__
                varname = varname[1]
            else: varname = varname[0] # node value

            # search vars in nodes/centers container
            containerNode = getNodeFromName1(z, container)
            if containerNode is None: return None
            else: # location of flow solution
                l = 'Vertex'
                l = getNodeFromType1(containerNode, 'GridLocation_t')
                if l is not None: ploc = getValue(l)
                else:
                    if name == __FlowSolutionNodes__: ploc = 'Vertex'
                    else: ploc = 'CellCenter'
            if loc == ploc: # variable is consistent with the container FlowSolution
                if loc == 'Vertex':  matchingLoc = 'Vertex'
                else: matchingLoc = 'FaceCenter' 
                f0 = getNodeFromName1(containerNode,varname)
                if f0 is not None:
                    if f0[3] == 'DataArray_t':
                        dataSetL = [] # [[BCRange1,fieldBC],[BCRange2,fieldBC2],[...],...]
                        for nobcdata in range(nbcdataset):
                            if allLocDS[nobcdata]==matchingLoc: # OK: locations are consistent
                                datas = getNodesFromName1(allBCDatas[nobcdata],varname) 
                                if datas != []:
                                    dataSetL+=[[allRanges[nobcdata], datas[0][1]]]
                        if dataSetL != []: # variable exists both in flow solution and in bcdataset
                            dataFS = [varname, loc, dataSetL]
                            containers.append(dataFS)
    if containers == []: return None
    else: return containers
    
#==============================================================================
# Automatically set the containers to the first found
# IN: t: pyTree ou base ou zone
# OUT: __GridCoordinates__, __FlowSolutionCenters__, ...
#==============================================================================
def autoSetContainers(t):
    global __FlowSolutionCenters__, __FlowSolutionNodes__
    zones = getZones(t)
    if len(zones) == 0: return
    z = zones[0]
    node = getNodeFromType1(z, 'GridCoordinates_t')
    if node is not None:
        __GridCoordinates__ = node[0]
    nodes = getNodesFromType1(z, 'FlowSolution_t')
    foundNode = 0; foundCenter = 0
    for n in nodes:
        loc = getNodeFromType1(n, 'GridLocation_t')
        data = getNodeFromType1(n, 'DataArray_t')
        if data is not None and data[1] is not None:
            if loc is not None and getValue(loc) == 'CellCenter':
                if foundCenter == 0: __FlowSolutionCenters__ = n[0]
                foundCenter += 1
            else:
                if foundNode == 0: __FlowSolutionNodes__ = n[0]
                foundNode += 1
    if __FlowSolutionNodes__ == __FlowSolutionCenters__:
        if foundNode == 0: __FlowSolutionNodes__ = ''
        elif foundCenter == 0: __FlowSolutionCenters__ = ''
        else: print('Warning: FlowSolutionNodes and FlowSolutionCenters have the same name (%s).'%__FlowSolutionNodes__)
    if foundNode > 1: print('Warning: multiple FlowSolutionNodes containers found (selected: %s)'%__FlowSolutionNodes__)
    if foundCenter > 1: print('Warning: multiple FlowSolutionCenters containers found (selected: %s)'%__FlowSolutionCenters__)
    return None

#==============================================================================
# IN: varName ('Density', 'centers:Density', 'nodes:Density', 
# {centers:Density]...
# OUT: varName, loc=0 (nodes), 1 (centers) 
#==============================================================================
def fixVarName(var):
    var = var.replace('{', '')
    var = var.replace('}', '')
    varS = var.split(':')
    if len(varS) == 1: return var, 0
    if varS[0] == 'centers': return varS[1], 1
    return varS[1], 0
