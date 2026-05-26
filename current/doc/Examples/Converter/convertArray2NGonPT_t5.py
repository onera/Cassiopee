# - convertArray2NGon(pyTree): topologic -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

# Simple box, Elements, with BCs
N = 10

def _addBCsBySubzone(t, nxblocks=1, is3D=True, quadFaces=True):
    Nz = N if is3D else 1
    Lz = 1. if is3D else 0.
    cartFunc = G.cartHexa if quadFaces else G.cartTetra
    subzone = cartFunc((0., 0., 0.), (0., 1., Lz), (1, N, Nz))
    C._addBC2Zone(t, 'inlet', 'BCInflow', subzone=subzone)
    if is3D: subzone = cartFunc((0., 0., 0.), (1., 1., 0.), (nxblocks*(N - 1) + 1, N, 1))
    else: subzone = cartFunc((0., 0., 0.), (1., 0., 0.), (nxblocks*(N - 1) + 1, 1, 1))
    C._addBC2Zone(t, 'wall', 'BCWallInviscid', subzone=subzone)
    if is3D: subzone = cartFunc((0., 0., Nz - 1.), (1., 1., 0.), (nxblocks*(N - 1) + 1, N, 1))
    else: subzone = cartFunc((0., N - 1., 0.), (1., 0., 0.), (nxblocks*(N - 1) + 1, 1, 1))
    C._addBC2Zone(t, 'farfield', 'BCFarfield', subzone=subzone)
    subzone = cartFunc(((N - 1.) * nxblocks, 0., 0.), (0., 1., Lz), (1, N, Nz))
    C._addBC2Zone(t, 'outlet', 'BCOutflow', subzone=subzone)
    return None

def _addData(t):
    C._initVars(t, "Density=1.05")
    C._initVars(t, "centers:Pressure=1.")

def _addBCDataAtFaceCenters(t):
    b = Internal.getBCNodesFromType(t, bndType='BCWall*')[0]
    er = Internal.getNodeFromName(b, Internal.__ELEMENTRANGE__)
    er = Internal.getValue(er)[0]
    nfaces = er[1] - er[0] + 1
    d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                              gridLocation='FaceCenter', parent=b)
    d = Internal.newBCData('BCNeumann', parent=d)
    Internal.newDataArray('Density', value=nfaces*[1.05], parent=d)
    Internal.newDataArray('MomentumX', value=nfaces*[1.], parent=d)
    Internal.newDataArray('MomentumY', value=nfaces*[0.], parent=d)
    Internal.newDataArray('MomentumZ', value=nfaces*[0.], parent=d)
    Internal.newDataArray('EnergyStagnationDensity', value=nfaces*[1.], parent=d)

def _convertStruct2NGon(method="geometric", addDataSets=False, api=1):
    C.clearAllNames()
    t = G.cart((0., 0., 0.), (1., 1., 1.), (N, N, N))
    C._addBC2Zone(t, 'inlet', 'BCInflow', wrange='imin')
    C._addBC2Zone(t, 'wall', 'BCWallInviscid', wrange='imax')
    C._addBC2Zone(t, 'farfield', 'BCFarfield', wrange='jmin')
    C._addBC2Zone(t, 'outlet', 'BCOutflow', wrange='jmax')
    _addData(t)
    if addDataSets: _addBCDataAtFaceCenters(t)
    tng = C.convertArray2NGon(t, recoverBC=True, method=method, api=api)
    test.testT(tng, 1)

def _convertQuad2NGon(method="geometric", addDataSets=False, api=1):
    C.clearAllNames()
    t = G.cartHexa((0., 0., 0.), (1., 1., 0.), (N, N, 1))
    _addBCsBySubzone(t, nxblocks=1, quadFaces=True, is3D=False)
    _addData(t)
    # _addBCDataAtNodes(t)
    if addDataSets: _addBCDataAtFaceCenters(t)
    tng = C.convertArray2NGon(t, recoverBC=True, method=method, api=api)
    test.testT(tng, 2)

def _convertHexa2NGon(method="geometric", addDataSets=False, api=1):
    C.clearAllNames()
    t = G.cartHexa((0., 0., 0.), (1., 1., 1.), (N, N, N))
    _addBCsBySubzone(t, nxblocks=1, quadFaces=True)
    _addData(t)
    if addDataSets: _addBCDataAtFaceCenters(t)
    tng = C.convertArray2NGon(t, recoverBC=True, method=method, api=api)
    test.testT(tng, 3)

def _convertME2NGon(method="geometric", addDataSets=False):
    C.clearAllNames()
    a = G.cartPyra((0., 0., 0.), (1., 1., 1.), (N, N, N))
    b = G.cartHexa((N-1., 0., 0.), (1., 1., 1.), (N, N, N))
    t = C.mergeConnectivity(a, b, boundary=0)
    _addBCsBySubzone(t, nxblocks=2, quadFaces=True)
    _addData(t)
    if addDataSets: _addBCDataAtFaceCenters(t)
    tng = C.convertArray2NGon(t, recoverBC=True, method=method, api=3)
    test.testT(tng, 4)

# 3D STRUCT
_convertStruct2NGon(method="topologic", addDataSets=False, api=3)

# 2D BE
_convertQuad2NGon(method="topologic", addDataSets=False, api=3)

# 3D BE
_convertHexa2NGon(method="topologic", addDataSets=False, api=1)

# 3D ME: pyra - hexa
_convertME2NGon(method="topologic", addDataSets=False)
