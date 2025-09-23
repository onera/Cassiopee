"""Transformation of pyTrees.
"""
#
# Python Interface to make basic transformations on pyTrees
#
from . import Transform
from . import transform
import numpy
__version__ = Transform.__version__

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Converter
except:
    raise ImportError("Transform.PyTree: requires Converter module.")

def collapse(a):
    """Collapse the smallest edge of each element for TRI arrays. Return a BAR.
    Usage: collapse(a)"""
    return C.TZGC1(a, 'nodes', True, Transform.collapse)

def _collapse(a):
    """Collapse the smallest edge of each element for TRI arrays. Return a BAR."""
    return C._TZGC1(a, 'nodes', True, Transform.collapse)

def cart2Cyl(t, center, axis, depth=0, thetaShift=0):
    """Transform a mesh in Cartesian coordinates to cylindrical coordinates.
    Usage: cart2Cyl(t, center, axis, depth, thetaShift)"""
    return C.TZGC1(t, 'nodes', True, Transform.cart2Cyl, center, axis, depth, thetaShift)

def _cart2Cyl(t, center, axis, depth=0, thetaShift=0):
    """Transform a mesh in Cartesian coordinates to cylindrical coordinates."""
    for z in Internal.getZones(t):
        transform._cart2CylZ(z, center, axis, depth, thetaShift,
                             Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__)
    return None

def cyl2Cart(t, center, axis):
    """Transform a mesh in Cylindrical coordinates to Cartesian coordinates.
    Usage: cyl2Cart(t, center, axis)"""
    return C.TZGC1(t, 'nodes', True, Transform.cyl2Cart, center, axis)

def _cyl2Cart(t, center, axis):
    """Transform a mesh in Cylindrical coordinates to Cartesian coordinates."""
    for z in Internal.getZones(t):
        transform._cyl2CartZ(z, center, axis, Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__)
    return None

def translate(t, transvect):
    """Translate a zone.
    Usage: translate(z, (v1,v2,v3))"""
    return C.TZGC3(t, 'nodes', False, Transform.translate, transvect)

def _translate(t, transvect):
    """Translate a zone."""
    C.__TZGC3(t, Transform._translate, transvect)
    return None

# Really in place - on coordinates only
def _rotate(t, center, arg1, arg2=None, vectors=[]):
    """Rotate a zone."""
    vectorsN = []; vectorsC = []
    for vect in vectors:
        if len(vect) == 3:
            loc = 0; vectname=[]
            for nov in range(3):
                spl = vect[nov].split(':')
                if len(spl) == 2:
                    vectname.append(spl[1])
                    if spl[0] == 'centers': loc += 1
                    else: loc += 4
                else: vectname.append(spl[0]); loc += 4
            if loc == 3: vectorsC += [vectname]
            elif loc == 12: vectorsN += [vectname]
    C.__TZA3(t, 'nodes', Transform._rotate, center, arg1, arg2, vectorsN)
    C.__TZA3(t, 'centers', Transform._rotate, center, arg1, arg2, vectorsC)
    return None

def rotate(t, center, arg1, arg2=None, vectors=[]):
    """Rotate a zone."""
    vectorsN = []; vectorsC = []
    for vect in vectors:
        if len(vect) == 3:
            loc = 0; vectname=[]
            for nov in range(3):
                spl = vect[nov].split(':')
                if len(spl) == 2:
                    vectname.append(spl[1])
                    if spl[0] == 'centers': loc += 1
                    else: loc += 4
                else: vectname.append(spl[0]); loc += 4
            if loc == 3: vectorsC += [vectname]
            elif loc == 12: vectorsN += [vectname]
    tp = C.TZA3(t, 'nodes', 'nodes', True, Transform.rotate, center, arg1, arg2, vectorsN)
    C.__TZA3(tp, 'centers', Transform._rotate, center, arg1, arg2, vectorsC)
    return tp

def homothety(a, center, alpha):
    """Make for a mesh defined by an array an homothety of center Xc and
    of factor alpha.
    Usage: homothety(a, (xc,yc,zc), alpha)"""
    return C.TZGC3(a, 'nodes', False, Transform.homothety, center, alpha)

def _homothety(a, center, alpha):
    """Make for a mesh defined by an array an homothety of center Xc and of factor alpha."""
    return C.__TZGC3(a, Transform._homothety, center, alpha)

def contract(a, center, dir1, dir2, alpha):
    """Contract a mesh around a plane defined by (center, dir1, dir2) and of factor alpha.
    Usage: contract(a, (xc,yc,zc), dir1, dir2, alpha)"""
    return C.TZGC3(a, 'nodes', False, Transform.contract, center, dir1, dir2, alpha)

def _contract(a, center, dir1, dir2, alpha):
    """Contract a mesh around a plane defined by (center, dir1, dir2) and of factor alpha.
    Usage: contract(a, (xc,yc,zc), dir1, dir2, alpha)"""
    return C.__TZGC3(a, Transform._contract, center, dir1, dir2, alpha)

def scale(a, factor=1., X=None):
    """Scale a mesh of given factor."""
    if X is None:
        try:
            import Generator.PyTree as G
            X = G.barycenter(a)
        except: pass
    return C.TZGC3(a, 'nodes', False, Transform.scale, factor, X)

def _scale(a, factor=1., X=None):
    """Scale a mesh of given factor."""
    if X is None:
        try:
            import Generator.PyTree as G
            X = G.barycenter(a)
        except: pass
    return C.__TZGC3(a, Transform._scale, factor, X)

def symetrize(a, point, vector1, vector2):
    """Make a symetry of mesh from plane passing by point and of director vector: vector1 and vector2.
    Usage: symetrize(a, (xc,yc,zc), (v1x,v1y,v1z), (v2x,v2y,v2z))"""
    return C.TZGC3(a, 'nodes', False, Transform.symetrize, point, vector1, vector2)

def _symetrize(a, point, vector1, vector2):
    """Make a symetry of mesh from plane passing by point and of director vector: vector1 and vector2.
    Usage: symetrize(a, (xc,yc,zc), (v1x,v1y,v1z), (v2x,v2y,v2z))"""
    return C.__TZGC3(a, Transform._symetrize, point, vector1, vector2)

def perturbate(a, radius, dim=3):
    """Perturbate a mesh randomly of radius
    Usage: perturbate(a, radius, dim)"""
    return C.TZA1(a, 'nodes', 'nodes', False, Transform.perturbate, radius, dim)

def _perturbate(a, radius, dim=3):
    """Perturbate a mesh randomly of radius
    Usage: perturbate(a, radius, dim)"""
    return C._TZA1(a, 'nodes', 'nodes', False, Transform.perturbate, radius, dim)

def smoothField(t, eps=0.1, niter=1, type=0, varNames=[]):
    """Smooth given fields."""
    return C.TZA3(t, 'nodes', 'nodes', False, Transform.smoothField, eps, niter, type, varNames)

def _smoothField(t, eps=0.1, niter=1, type=0, varNames=[]):
    """Smooth given fields."""
    varNodes = []; varCenters = []
    for v in varNames:
        s = v.split(':')
        if len(s) == 2:
            if s[0] == 'centers': varCenters.append(s[1])
            else: varNodes.append(s[0])
        else: varNodes.append(v)
    if varNodes != []: C.__TZA3(t, 'nodes', Transform._smoothField, eps, niter, type, varNodes)
    if varCenters != []: C.__TZA3(t, 'centers', Transform._smoothField, eps, niter, type, varCenters)
    return None

def smooth(t, eps=0.5, niter=4, type=0, fixedConstraints=[],
           projConstraints=[], delta=1., point=(0,0,0), radius=-1.):
    """Smooth a mesh with a Laplacian.
    Usage: smooth(t, eps, niter, type, fixedConstraints, projConstraints, delta)"""
    tp = Internal.copyRef(t)
    _smooth(tp, eps, niter, type, fixedConstraints,
            projConstraints, delta, point, radius)
    return tp

def _smooth(t, eps=0.5, niter=4, type=0, fixedConstraints=[],
            projConstraints=[], delta=1., point=(0,0,0), radius=-1.):
    """Smooth a mesh with a Laplacian.
    Usage: smooth(t, eps, niter, type, fixedConstraints, projConstraints, delta)"""
    if fixedConstraints != []:
        c = []
        for z in fixedConstraints:
            c.append(C.getFields(Internal.__GridCoordinates__, z)[0])
        fixedConstraints = c
    if projConstraints != []:
        c = []
        for z in projConstraints:
            c.append(C.getFields(Internal.__GridCoordinates__, z)[0])
        projConstraints = c
    zones = Internal.getZones(t)
    coords = []
    for z in zones: coords.append(C.getFields(Internal.__GridCoordinates__, z)[0])
    coordsp = Transform.smooth(coords, eps, niter, type,
                               fixedConstraints, projConstraints, delta,
                               point, radius)
    C.setFields(coordsp, t, 'nodes', writeDim=False)
    return None

def deform(t, vector=['dx','dy','dz']):
    """Deform surface by moving surface of the vector (dx, dy, dz).
    Usage: deform(t, vector=['dx','dy','dz'])"""
    tp = Internal.copyRef(t)
    _deform(tp, vector)
    return tp

def _deform(t, vector=['dx','dy','dz']):
    """Deform surface by moving surface of the vector (dx, dy, dz).
    Usage: deform(t, vector=['dx','dy','dz'])"""
    if len(vector) != 3: raise ValueError("deform: 3 variables are required.")
    return C._TZA1(t, 'nodes', 'nodes', False, Transform.deform, vector)

def deformNormals(t, alpha, niter=1):
    """Deform a a surface of alpha times the surface normals.
    alpha is a string denoting the alpha variable in a.
    Usage: deformNormals(t, alpha, niter)"""
    tp = Internal.copyRef(t)
    _deformNormals(tp, alpha, niter)
    return tp

def _deformNormals(t, alpha, niter=1):
    """Deform a a surface of alpha times the surface normals."""
    alphat = C.getField(alpha, t)
    res = C.TLAGC(t, Transform.deformNormals, alphat, niter, None, None)
    C.setFields(res[0], t, 'nodes')
    return None

def deformPoint(a, xyz, dxdydz, depth, width):
    """Deform mesh by moving point (x,y,z) of a vector (dx, dy, dz).
    Usage: deformPoint(a, (x,y,z), (dx,dy,dz), width, depth)"""
    return C.TZGC1(a, 'nodes', True, Transform.deformPoint, xyz, dxdydz, depth, width)

def _deformPoint(a, xyz, dxdydz, depth, width):
    """Deform mesh by moving point (x,y,z) of a vector (dx, dy, dz)."""
    return C._TZGC1(a, 'nodes', False, Transform.deformPoint, xyz, dxdydz, depth, width)

def deformMesh(a, surfDelta, beta=4., type='nearest'):
    """Deform a mesh a wrt surfDelta defining surface grids and deformation vector on it.
    Usage: deformMesh(a, surfDelta, beta, type)"""
    info = C.getAllFields(surfDelta, 'nodes')
    return C.TZA1(a, 'nodes', 'nodes', True, Transform.deformMesh, info, beta, type)

def _deformMesh(a, surfDelta, beta=4., type='nearest'):
    """Deform a mesh a wrt surfDelta defining surface grids and deformation vector on it."""
    info = C.getAllFields(surfDelta, 'nodes')
    return C._TZA1(a, 'nodes', 'nodes', True, Transform.deformMesh, info, beta, type)

def controlPoints(a, N=(2,2,2)):
    """Create control points for free form."""
    import Generator.PyTree as G
    bb = G.bbox(a)
    b = G.cart((bb[0],bb[1],bb[2]),(bb[3]-bb[0],bb[4]-bb[1],bb[5]-bb[2]),N)
    C._addVars(b, ['dx','dy','dz'])
    return b

def freeForm(a, controlPoints):
    """Coupute free form deformation vector."""
    b = Internal.copyRef(a)
    _freeForm(b, controlPoints)
    return b

def _freeForm(a, controlPoints):
    """Coupute free form deformation vector."""
    array = C.getAllFields(a, 'nodes')
    control = C.getAllFields(controlPoints, 'nodes')
    Transform.freeForm(array, control)
    return None

def join(t, t2=None, tol=1.e-10):
    """Join two zones in one or join a list of zones in one.
    Usage: join(t,t2) or join(t)"""
    nodes = Internal.getZones(t)
    if len(nodes) == 0: return []
    allBCInfos = C.extractBCInfo(t)

    if t2 is not None:
        nodes += Internal.getZones(t2)
        allBCInfos += C.extractBCInfo(t2)
    Internal._orderFlowSolution(nodes, loc='both')

    fieldn = C.getAllFields(nodes, 'nodes')
    fieldc = []
    for f in C.getAllFields(nodes, 'centers'):
        if f != []: fieldc.append(f)
    res = Transform.join(fieldn, arrayc=fieldc, tol=tol)
    if not isinstance(res[0], list): # join sans les centres
        z = C.convertArrays2ZoneNode('join', [res])
    else: # avec les centres res=[an,ac]
        z = C.convertArrays2ZoneNode('join', [res[0]])
        z = C.setFields([res[1]], z, 'centers')
    if Internal.getZoneType(z) == 1: # structured
        z = C.identifyBC(z, allBCInfos, tol)
    return z

def merge(t, sizeMax=1000000000, dir=0, tol=1.e-10, alphaRef=180., mergeBCs=False):
    """Merge a list of matching zones.
    Usage: merge(t, sizeMax, dir, tol, alphaRef)"""
    Internal._orderFlowSolution(t, loc='both')
    allBCInfos = C.extractBCInfo(t)
    fieldn = C.getAllFields(t, 'nodes')
    fieldc = []
    for f in C.getAllFields(t, 'centers'):
        if f != []: fieldc.append(f)
    res = Transform.merge(fieldn, Ac=fieldc, sizeMax=sizeMax,
                          dir=dir, tol=tol, alphaRef=alphaRef)
    zones = []
    if fieldc == []: # merge sans les centres
        for i in res: zones += [C.convertArrays2ZoneNode('zone',[i])]
    else:
        nzones = len(res[0])
        for noi in range(nzones):
            z = C.convertArrays2ZoneNode('zone',[res[0][noi]])
            z = C.setFields([res[1][noi]], z, 'centers')
            zones += [z]
    zones = C.identifyBC(zones, allBCInfos, tol)

    # Merge BCs by type (keeping families) and connectMatch
    if mergeBCs:
        import Connector.PyTree as X
        zones = X.connectMatch(zones)
        for z in zones:
            bct = {}; fct = {}
            zbcs = Internal.getNodesFromType1(z, 'ZoneBC_t')
            for zbc in zbcs:
                bcs = Internal.getNodesFromType1(zbc, 'BC_t')
                for bc in bcs:
                    btype = Internal.getValue(bc)
                    if btype == 'FamilySpecified':
                        fn = Internal.getNodeFromType1(bc, 'FamilyName_t')
                        btype = Internal.getValue(fn)
                        fct[btype] = bc[0].split('.')[0]
                    else:
                        bct[btype] = bc[0].split('.')[0]
            # Rm classical BC and fill
            for k in bct:
                Internal._rmNodesFromValue(z, k)
                C._fillEmptyBCWith(z, bct[k], k)
            # Rm Family BC and fill
            for k in fct:
                zbcs = Internal.getNodesFromType1(z, 'ZoneBC_t')
                for zbc in zbcs:
                    out = []
                    for i in zbc[2]:
                        if Internal.getValue(i) == 'FamilySpecified':
                            fn = Internal.getNodeFromType1(i, 'FamilyName_t')
                            btype = Internal.getValue(fn)
                            if btype != k: out.append(i)
                        else: out.append(i)
                    zbc[2] = out
                C._fillEmptyBCWith(z, fct[k], 'FamilySpecified:'+k)
    return zones

def mergeCart(t, sizeMax=1000000000, tol=1.e-10):
    """Merge Cartesian grids defined as a list of zones using Rigby algorithm.
    Usage: mergeCart(t, sizeMax, tol)"""
    Internal._orderFlowSolution(t, loc='both')
    allBCInfos = C.extractBCInfo(t)
    A = C.getAllFields(t, 'nodes', api=3)
    A = Transform.mergeCart(A, sizeMax, tol)
    for noz in range(len(A)):
        A[noz] = C.convertArrays2ZoneNode('Cart',[A[noz]])
    A = C.identifyBC(A, allBCInfos, tol)
    return A

def patch(t1, t2, position=None, nodes=None, order=None):
    """Patch mesh2 defined by t2 in mesh1 defined by t1
    at position (i,j,k).
    Usage: patch(t1, t2, (i,j,k))"""
    tp1 = Internal.copyRef(t1)
    _patch(tp1, t2, position=position, nodes=nodes, order=order)
    return tp1

def _patch(t1, t2, position=None, nodes=None, order=None):
    """Patch mesh2 defined by t2 in mesh1 defined by t1 at position (i,j,k)."""
    zones1 = Internal.getZones(t1)
    zones2 = Internal.getZones(t2)
    for z1,z2 in zip(zones1, zones2):
        a2 = C.getAllFields(z2, 'nodes')[0]
        C._TZA1(z1, 'nodes', 'nodes', True, Transform.patch, a2, position, nodes, order)
    return None

#===============
# oneovernBC
#===============
def _oneovernBC__(t, N):
    import math
    C._rmBCOfType(t, 'BCMatch'); C._rmBCOfType(t, 'BCNearMatch')
    nodes = Internal.getZones(t)
    for z in nodes:
        wins = Internal.getNodesFromType(z, 'BC_t')
        for w in wins:
            w0 = Internal.getNodeFromName1(w, 'PointRange')
            (parent, d) = Internal.getParentOfNode(z, w0)
            w0 = w0[1]
            i1 = w0[0,0]; j1 = w0[1,0]; k1 = w0[2,0]
            i2 = w0[0,1]; j2 = w0[1,1]; k2 = w0[2,1]
            addi1=0; addi2=0; addj1=0; addj2=0; addk1=0; addk2=0
            i1n = math.floor(1.*(i1-1)/N[0])+1
            i2n = math.floor(1.*(i2-1)/N[0])+1
            j1n = math.floor(1.*(j1-1)/N[1])+1
            j2n = math.floor(1.*(j2-1)/N[1])+1
            k1n = math.floor(1.*(k1-1)/N[2])+1
            k2n = math.floor(1.*(k2-1)/N[2])+1
            if i1 - i1n*N[0] != 1-N[0]: addi1 = 1
            if i2 - i2n*N[0] != 1-N[0]: addi2 = 1
            if j1 - j1n*N[1] != 1-N[1]: addj1 = 1
            if j2 - j2n*N[1] != 1-N[1]: addj2 = 1
            if k1 - k1n*N[2] != 1-N[2]: addk1 = 1
            if k2 - k2n*N[2] != 1-N[2]: addk2 = 1
            i1n+=addi1;i2n+=addi2;j1n+=addj1;j2n+=addj2;k1n+=addk1;k2n+=addk2
            range0 = [int(i1n),int(i2n),int(j1n),int(j2n),int(k1n),int(k2n)]
            r2 = Internal.window2Range(range0)
            parent[2][d][1] = r2

        connect = Internal.getNodesFromType(z, 'ZoneGridConnectivity_t')
        for cn in connect:
            wins = Internal.getNodesFromName2(cn, 'PointRange')
            for w in wins:
                (parent, d) = Internal.getParentOfNode(cn, w)
                w0 = w[1]
                i1 = w0[0,0]; j1 = w0[1,0]; k1 = w0[2,0]
                i2 = w0[0,1]; j2 = w0[1,1]; k2 = w0[2,1]
                addi1=0; addi2=0; addj1=0; addj2=0;addk1=0; addk2=0
                i1n = math.floor(1.*(i1-1)/N[0])+1
                i2n = math.floor(1.*(i2-1)/N[0])+1
                j1n = math.floor(1.*(j1-1)/N[1])+1
                j2n = math.floor(1.*(j2-1)/N[1])+1
                k1n = math.floor(1.*(k1-1)/N[2])+1
                k2n = math.floor(1.*(k2-1)/N[2])+1
                if i1 - i1n*N[0] != 1-N[0]: addi1 = 1
                if i2 - i2n*N[0] != 1-N[0]: addi2 = 1
                if j1 - j1n*N[1] != 1-N[1]: addj1 = 1
                if j2 - j2n*N[1] != 1-N[1]: addj2 = 1
                if k1 - k1n*N[2] != 1-N[2]: addk1 = 1
                if k2 - k2n*N[2] != 1-N[2]: addk2 = 1
                i1n+=addi1;i2n+=addi2;j1n+=addj1;j2n+=addj2;k1n+=addk1;k2n+=addk2
                range0 = [int(i1n),int(i2n),int(j1n),int(j2n),int(k1n),int(k2n)]
                r2 = Internal.window2Range(range0)
                parent[2][d][1] = r2
    return t

#==============================================================================
# oneovern tree
#==============================================================================
def oneovern(t, N):
    """Take one over N points from mesh.
    Usage: oneovern(t, (Ni,Nj,Nk))"""
    tp = Internal.copyRef(t)
    _oneovern(tp, N)
    return tp

def _oneovern(t, N):
    """Take one over N points from mesh."""
    C._TZA1(t, 'nodes', 'nodes', True, Transform.oneovern, N, 1)
    C._TZA1(t, 'centers', 'centers', False, Transform.oneovern, N, 0)
    _oneovernBC__(t, N)
    return None

#==============================================================================
# Subzone des indices d'une fenetre
# IN: imin, imax, jmin, jmax, kmin, kmax: indices de subzone
# IN: w: window a subzoner
# OUT: [i1, i2, j1, j2, k1, k2]: indice de la window subzonee
#==============================================================================
def getBCRange__(w, imin, imax, jmin, jmax, kmin, kmax):
    i1 = w[0,0]; j1 = w[1,0]; k1 = w[2,0]
    i2 = w[0,1]; j2 = w[1,1]; k2 = w[2,1]
    io1 = i1; io2 = i2; jo1 = j1; jo2 = j2; ko1 = k1; ko2 = k2
    i1=min(io1,io2); i2=max(io1,io2);j1=min(jo1,jo2);j2=max(jo1,jo2);k1=min(ko1,ko2);k2=max(ko1,ko2)
    if imin >= i1 and imax <= i2: io1 = 1; io2 = imax-imin+1
    elif imin >= i1 and imin <= i2: io1 = 1; io2 = i2-imin+1
    elif imax >= i1 and imax <= i2: io2 = imax
    if imin < i1: io1 = i1-imin+1; io2 = io2-imin + 1

    if jmin >= j1 and jmax <= j2: jo1 = 1; jo2 = jmax-jmin+1
    elif jmin >= j1 and jmin <= j2: jo1 = 1; jo2 = j2-jmin+1
    elif jmax >= j1 and jmax <= j2: jo2 = jmax
    if jmin < j1: jo1 = j1-jmin+1; jo2 = jo2-jmin+1

    if kmin >= k1 and kmax <= k2: ko1 = 1; ko2 = kmax-kmin+1
    elif kmin >= k1 and kmin <= k2: ko1 = 1; ko2 = k2-kmin+1
    elif kmax >= k1 and kmax <= k2: ko2 = kmax
    if kmin < k1: ko1 = k1-kmin+1; ko2 = ko2-kmin+1
    if i1 == i2:
        if i1 >= imax: io1 = imax-imin+1; io2 = imax-imin+1
        elif i1 < imin: io1 = 1; io2 = 1
    if j1 == j2:
        if j1 >= jmax: jo1 = jmax-jmin+1; jo2 = jmax-jmin+1
        elif j1 < jmin: jo1 = 1; jo2 = 1
    if k1 == k2:
        if k1 >= kmax: ko1 = kmax-kmin+1; ko2 = kmax-kmin+1
        elif k1 < kmin: ko1 = 1; ko2 = 1
    return [min(io1,io2), max(io1,io2), min(jo1,jo2), max(jo1,jo2), min(ko1,ko2), max(ko1,ko2)]

# Retourne 0 si la fenetre initiale est interieure a la zone subzonee
def isWindowInSubzone__(w, dim, imin, imax, jmin, jmax, kmin, kmax,
                        ni0, nj0, nk0):
    io1 = w[0,0]; jo1 = w[1,0]; ko1 = w[2,0]
    io2 = w[0,1]; jo2 = w[1,1]; ko2 = w[2,1]
    i1 = min(io1,io2); i2 = max(io1,io2)
    j1 = min(jo1,jo2); j2 = max(jo1,jo2)
    k1 = min(ko1,ko2); k2 = max(ko1,ko2)
    isout = 0
    if dim == 2:
        if ((i1 == i2 and i1 == 1 and imin > 1) or \
            (i1 == i2 and i1 == ni0 and imax < ni0) or \
            (j1 == j2 and j1 == 1 and jmin > 1) or \
                (j1 == j2 and j1 == nj0 and jmax < nj0)):
            isout = 1
    else:
        if ((i1 == i2 and i1 == 1 and imin > 1) or \
            (i1 == i2 and i1 == ni0 and imax < ni0) or \
            (j1 == j2 and j1 == 1 and jmin > 1) or \
            (j1 == j2 and j1 == nj0 and jmax < nj0) or \
            (k1 == k2 and k1 == 1 and kmin > 1) or \
                (k1 == k2 and k1 == nk0 and kmax < nk0)):
            isout = 1

    if imax < i1 or imin > i2 or jmax < j1 or jmin > j2 or kmax < k1 or kmin > k2: isout = 1
    return isout

# subzone les BCs de z de l'arbre t
def subzoneBCStruct__(t, z, dim, imin, imax, jmin, jmax, kmin, kmax,
                      dim0, nip, njp, nkp, ni0, nj0, nk0):
    wins = Internal.getNodesFromType2(z, 'BC_t')
    for w in wins:
        pr = Internal.getNodeFromName1(w, 'PointRange')
        isout = isWindowInSubzone__(pr[1], dim,
                                    imin, imax, jmin, jmax, kmin, kmax,
                                    ni0, nj0, nk0)

        (parent, d) = Internal.getParentOfNode(w, pr)
        if isout == 1:
            (parentw, dw) = Internal.getParentOfNode(z, w)
            del parentw[2][dw]
        else:
            range0 = getBCRange__(pr[1], imin, imax, jmin, jmax, kmin, kmax)
            notvalid = 0
            if dim == 3:
                if range0[0] == range0[1] and range0[2] == range0[3]: notvalid = 1
                elif range0[0] == range0[1] and range0[4] == range0[5]: notvalid = 1
                elif range0[2] == range0[3] and range0[4] == range0[5]: notvalid = 1
            elif dim == 2:
                if range0[0] == range0[1] and range0[2] == range0[3]: notvalid = 1

            if notvalid == 1:
                (parentw, dw) = Internal.getParentOfNode(z, w)
                del parentw[2][dw]
            else:
                bndType = w[1]
                if dim == 2 and dim0 == 3:
                    if imin == imax and (imin == 1 or imin == nip): # cas
                        range0 = [range0[2], range0[3], range0[4], range0[5],1,1]
                    elif jmin == jmax and (jmin == 1 or jmin == njp):
                        range0 = [range0[0], range0[1], range0[4], range0[5],1,1]
                    elif kmin == kmax and (kmin == 1 or kmin == nkp):
                        range0 = [range0[0], range0[1], range0[2], range0[3],1,1]
                r2 = Internal.window2Range(range0)
                parent[2][d][1] = r2
                w[0] = C.getBCName(w[0])

# Subzone la GridConnectivity de type BCOverlap
def subzoneGCStruct__(z, dim, imin, imax, jmin, jmax, kmin, kmax, \
                      dim0, nip, njp, nkp, ni0, nj0, nk0):
    nodes = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    if nodes == []: return z
    ranges = []; ddDnrs = []; bcnames = []
    for i in nodes:
        pr = Internal.getNodeFromName1(i, 'PointRange')
        if pr is None: break
        isout = isWindowInSubzone__(pr[1], dim,
                                    imin, imax, jmin, jmax, kmin, kmax,
                                    ni0, nj0, nk0)

        r = Internal.getNodeFromType1(i, 'GridConnectivityType_t')
        #(parent, d) = Internal.getParentOfNode(z, i)
        if r is not None:
            val = Internal.getValue(r)
            if val == 'Overset':
                if isout == 0:
                    DDDnrName=None
                    for UDN in Internal.getNodesFromType1(i,'UserDefinedData_t'):
                        if Internal.getNodeFromName1(UDN,'doubly_defined') is not None:
                            DDDnrName=Internal.getValue(i)
                            break

                    w0 = pr[1]
                    range0 = getBCRange__(w0, imin, imax, jmin, jmax, kmin, kmax)
                    notvalid = 0
                    if dim == 3:
                        if range0[0] == range0[1] and range0[2] == range0[3]: notvalid = 1
                        elif range0[0] == range0[1] and range0[4] == range0[5]: notvalid = 1
                        elif range0[2] == range0[3] and range0[4] == range0[5]: notvalid = 1
                    elif dim == 2:
                        if range0[0] == range0[1] and range0[2] == range0[3]: notvalid = 1
                    if notvalid == 0:
                        if dim == 2 and dim0 == 3:
                            if imin == imax and (imin == 1 or imin == nip): # cas
                                range0 = [range0[2], range0[3], range0[4], range0[5],1,1]
                            elif jmin == jmax and (jmin == 1 or jmin == njp):
                                range0 = [range0[0], range0[1], range0[4], range0[5],1,1]
                            elif kmin == kmax and (kmin == 1 or kmin == nkp):
                                range0 = [range0[0], range0[1], range0[2], range0[3],1,1]
                        ranges.append(range0)
                        ddDnrs.append(DDDnrName)
                        bcnames.append(Internal.getName(i))
    C._rmBCOfType(z, 'BCOverlap')
    for nor in range(len(ranges)):
        if ddDnrs[nor] is None:
            C._addBC2Zone(z, C.getBCName('overlap'), 'BCOverlap', ranges[nor])
        else: # doubly defined
            C._addBC2Zone(z, C.getBCName(bcnames[nor]), 'BCOverlap', wrange=ranges[nor],\
                          zoneDonor=[ddDnrs[nor]], rangeDonor='doubly_defined')
    return z

def subzone(t, minIndex, maxIndex=None, type=None):
    """Take a subzone of mesh.
    Usage: subzone(t, (imin,jmin,kmin), (imax,jmax,kmax))"""
    if maxIndex is None: return subzoneUnstruct__(t, minIndex, type)
    else: return subzoneStruct__(t, minIndex, maxIndex)

def subzoneUnstruct__(t, indices, type):
    tp = Internal.copyRef(t)
    nodes = Internal.getZones(tp)
    for z in nodes:
        dimz = Internal.getZoneDim(z)
        fc = C.getFields(Internal.__GridCoordinates__, z, api=1)[0]
        fa = C.getFields(Internal.__FlowSolutionNodes__, z, api=1)[0]
        fb = C.getFields(Internal.__FlowSolutionCenters__, z, api=1)[0]
        if fa != []: fc = Converter.addVars([fc, fa])
        if fb == []: # no flow sol at centers
            nodes = Transform.subzone(fc, indices, type=type)
            C.setFields([nodes], z, 'nodes')
        else:
            if dimz[0] == 'Structured': # faceList as global indices of structured interfaces
                if type == 'faces':
                    [nodes, centers] = Transform.transform.subzoneStructIntBoth(fc, fb, indices)
                    C.setFields([nodes], z, 'nodes')
                    C.setFields([centers], z, 'centers')
                else:
                    raise TypeError("subzone with list (loc='both'): not yet implemented for structured zones.")
            else:
                if type == 'faces':
                    [nodes, centers] = Transform.transform.subzoneFacesBoth(fc, fb, indices)
                    C.setFields([nodes], z, 'nodes')
                    C.setFields([centers], z, 'centers')
                elif type == 'elements':
                    [nodes, centers] = Transform.transform.subzoneElementsBoth(fc, fb, indices)
                    C.setFields([nodes], z, 'nodes')
                    C.setFields([centers], z, 'centers')
                else:
                    [nodes, centers] = Transform.transform.subzoneUnstructBoth(fc, fb, indices)
                    C.setFields([nodes], z, 'nodes')
                    C.setFields([centers], z, 'centers')
        z[0] = C.getZoneName(z[0])
    return tp

def subzoneStruct__(t, minIndex, maxIndex):
    # indices pour les centres
    imin = minIndex[0]; jmin = minIndex[1]; kmin = minIndex[2]
    imax,jmax,kmax = [max(1, val-1) if val > -1 else val for val in maxIndex]
    if imin == maxIndex[0] and imin != 1: imin = imax
    if jmin == maxIndex[1] and jmin != 1: jmin = jmax
    if kmin == maxIndex[2] and kmin != 1: kmin = kmax
    t2 = C.TZA(t, 'both', 'both',
               Transform.subzone, Transform.subzone,
               minIndex, maxIndex,
               (imin, jmin, kmin), (imax,jmax,kmax))
    C._rmBCOfType(t2, 'BCMatch'); C._rmBCOfType(t2, 'BCNearMatch')
    imin = minIndex[0]; imax = maxIndex[0]
    jmin = minIndex[1]; jmax = maxIndex[1]
    kmin = minIndex[2]; kmax = maxIndex[2]
    nodes = Internal.getZones(t)
    nodes2 = Internal.getZones(t2)
    noz = 0
    for z2 in nodes2:
        z = nodes[noz]; dimt0 = Internal.getZoneDim(z)
        if imax < 0: imax = dimt0[1]+imax+1
        if jmax < 0: jmax = dimt0[2]+jmax+1
        if kmax < 0: kmax = dimt0[3]+kmax+1
        (parent, nb) = Internal.getParentOfNode(t2, z2)
        dimt = Internal.getZoneDim(z2)
        ni0 = dimt[1]; nj0 = dimt[2]; nk0 = dimt[3]; dim = dimt[4]
        z2[0] = C.getZoneName(z2[0]) # ?? pourquoi changer le nom?
        subzoneBCStruct__(t2, z2, dim, imin, imax, jmin, jmax, kmin, kmax, \
                          dimt0[4], dimt0[1], dimt0[2], dimt0[3], ni0, nj0, nk0)
        z2 = subzoneGCStruct__(z2, dim, imin, imax, jmin, jmax, kmin, kmax, \
                               dimt0[4], dimt0[1], dimt0[2], dimt0[3], ni0, nj0, nk0)
        noz += 1
        if parent is not None:
            if Internal.isStdNode(t2) == 0: parent[nb] = z2
            else: parent[2][nb] = z2
        else: t2 = z2
    return t2

# intersection de fenetres
# win1 = fenetre par rapport a un bloc B
# win2 = fenetre par rapport a un bloc B
# si ret=0: retourne fenetre commune par rapport a B
# si ret=1: retourne fenetre commune par rapport a win1
def intersectWins__(win1, win2, ret=0):
    fi1,fi2,fj1,fj2,fk1,fk2 = win1
    si1,si2,sj1,sj2,sk1,sk2 = win2
    if fi2 < si1: return None
    if si2 < fi1: return None
    if fj2 < sj1: return None
    if sj2 < fj1: return None
    if fk2 < sk1: return None
    if sk2 < fk1: return None
    oi1 = max(fi1,si1)
    oi2 = min(fi2,si2)
    oj1 = max(fj1,sj1)
    oj2 = min(fj2,sj2)
    ok1 = max(fk1,sk1)
    ok2 = min(fk2,sk2)
    dim1 = 2
    if fk1 != fk2: dim1 = 3
    dim2 = 2
    if sk1 != sk2: dim2 = 3
    dim = max(dim1, dim2)
    if dim == 2: # avoid degenerated windows in 2D
        if oi1 == oi2 and oj1 == oj2 and ok1 == ok2: return None
    if dim == 3: # avoid degenerated windows in 3D
        if oi1 == oi2 and oj1 == oj2: return None
        if oi1 == oi2 and ok1 == ok2: return None
        if oj1 == oj2 and ok1 == ok2: return None

    if ret == 0: return [oi1,oi2,oj1,oj2,ok1,ok2]
    else: return [oi1-fi1+1,oi2-fi1+1,oj1-fj1+1,oj2-fj1+1,ok1-fk1+1,ok2-fk1+1]

# composition de fenetres
# win1 est une fenetre par rapport a win2
# win2 est une fenetre par rapport a B
# Retourne la fenetre 1 par rapport a B
def composeWins__(win1, win2):
    fi1,fi2,fj1,fj2,fk1,fk2 = win1
    si1,si2,sj1,sj2,sk1,sk2 = win2
    return [si1+fi1-1,si1+fi2-1,sj1+fj1-1,sj1+fj2-1,sk1+fk1-1,sk1+fk2-1]

# Transforme un trirac en shift matrix
def shiftMatrix__(trirac):
    m = numpy.zeros((9), dtype=numpy.float64)
    s = trirac.size
    for c in range(s):
        if trirac[c] == 1: m[0+3*c] = 1
        elif trirac[c] == -1: m[0+3*c] = -1
        elif trirac[c] == 2: m[1+3*c] = 1
        elif trirac[c] == -2: m[1+3*c] = -1
        elif trirac[c] == 3: m[2+3*c] = 1
        elif trirac[c] == -3: m[2+3*c] = -1
    if s == 2: m[2+3*2] = 1
    return m

def triracopp__(trirac):
    m = numpy.zeros((3,3),dtype=Internal.E_NpyInt)
    for no in range(3):
        signt = (trirac[no]>0)*1+(-1*(trirac[no]<0))
        i1 = abs(trirac[no])-1
        m[i1,no] = signt
    mopp = numpy.linalg.inv(m)
    triracopp=[1,2,3]
    for i in range(3):
        for no in range(3):
            if mopp[no,i] != 0:
                triracopp[i] = int(mopp[no,i]*(no+1))
    return triracopp

# Retourne un indice sur le donneur
# win: fenetre sur B (raccord match)
# winDonor: fenetre correspondante a win1 sur opp(B)
# trirac: trirac de B vers donor
# i,j,k: indices sur B
# retourne: indices de ijk sur opp(B)
def donorIndex__(win, winDonor, trirac, tupleIJK):
    m = shiftMatrix__(trirac)
    (i,j,k) = tupleIJK
    fi1,fi2,fj1,fj2,fk1,fk2 = win
    si1,si2,sj1,sj2,sk1,sk2 = winDonor

    x0 = si1; y0 = sj1; z0 = sk1
    if m[0+3*0]<-0.5 or m[0+3*1]<-0.5 or m[0+3*2]<-0.5: x0 = si2
    if m[1+3*0]<-0.5 or m[1+3*1]<-0.5 or m[1+3*2]<-0.5: y0 = sj2
    if m[2+3*0]<-0.5 or m[2+3*1]<-0.5 or m[2+3*2]<-0.5: z0 = sk2

    ip = x0+m[0+3*0]*(i-fi1)+m[0+3*1]*(j-fj1)+m[0+3*2]*(k-fk1)
    jp = y0+m[1+3*0]*(i-fi1)+m[1+3*1]*(j-fj1)+m[1+3*2]*(k-fk1)
    kp = z0+m[2+3*0]*(i-fi1)+m[2+3*1]*(j-fj1)+m[2+3*2]*(k-fk1)

    return (int(ip),int(jp),int(kp))

#=====================================================================
# Retourne les raccords matchs d'une zone z
# Si donorName est specifie, retourne les raccords matchs correspondants
# uniquement a ce donneur
#=====================================================================
def getBCMatchs__(z, donorName=None):
    bcs = []
    gc = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')

    if donorName is None: # return all BCMatchs
        if gc is not None:
            bcs = Internal.getNodesFromType1(gc, 'GridConnectivity1to1_t')
        return bcs

    # Return BCMatchs with given donorName
    if gc is not None:
        l = Internal.getNodesFromType1(gc, 'GridConnectivity1to1_t')
        for bc in l:
            oppBlock = Internal.getValue(bc)
            if oppBlock == donorName: bcs.append(bc)
    return bcs

#=====================================================================
# Retourne les datas d'un BCMatch : oppBlockName, range, donorRange, trf
# et infos periodiques si existantes
#=====================================================================
def getBCMatchData__(bc):
    oppBlock = Internal.getValue(bc)
    rnge = Internal.getNodeFromName1(bc, 'PointRange')
    rnge = Internal.range2Window(rnge[1])
    donor = Internal.getNodeFromName1(bc, 'PointRangeDonor')
    donor = Internal.range2Window(donor[1])
    trf = Internal.getNodeFromName1(bc, 'Transform')
    trf = Internal.getValue(trf)
    if trf.size == 2: # met toujours le transform en 3D
        trf2 = numpy.empty(3, dtype=Internal.E_NpyInt)
        trf2[0:2] = trf[0:2]; trf2[2] = 3
        trf = trf2
    # Get periodic data if any
    periodic = None
    gcp = Internal.getNodeFromType1(bc, 'GridConnectivityProperty_t')
    if gcp is not None:
        pr = Internal.getNodeFromType1(gcp, 'Periodic_t')
        if pr is not None:
            rotCenter = Internal.getNodeFromName1(pr, 'RotationCenter')
            translation = Internal.getNodeFromName1(pr, 'Translation')
            rotationAngle = Internal.getNodeFromName1(pr, 'RotationAngle')
            if rotCenter is not None: rotCenter = Internal.getValue(rotCenter)
            else: rotCenter = numpy.zeros((3), numpy.float64)
            if translation is not None: translation = Internal.getValue(translation)
            else: translation = numpy.zeros((3), numpy.float64)
            if rotationAngle is not None: rotationAngle = Internal.getValue(rotationAngle)
            else: rotationAngle = numpy.zeros((3), numpy.float64)
            periodic = [translation, rotCenter, rotationAngle]
    return (oppBlock, rnge, donor, trf, periodic)

# Cree le match interne entre z1 et z2 pour un split en dir
def _createInternalBCMatch(z1, z2, dir):
    if dir == 1:
        C._addBC2Zone(z1, 'match', 'BCMatch', 'imax', z2, 'imin', [1,2,3])
        C._addBC2Zone(z2, 'match', 'BCMatch', 'imin', z1, 'imax', [1,2,3])
    elif dir == 2:
        C._addBC2Zone(z1, 'match', 'BCMatch', 'jmax', z2, 'jmin', [1,2,3])
        C._addBC2Zone(z2, 'match', 'BCMatch', 'jmin', z1, 'jmax', [1,2,3])
    else:
        C._addBC2Zone(z1, 'match', 'BCMatch', 'kmax', z2, 'kmin', [1,2,3])
        C._addBC2Zone(z2, 'match', 'BCMatch', 'kmin', z1, 'kmax', [1,2,3])
    return None

# Delete dans t les BCs match referencant zname
def _deleteBCMatchRef(t, zname):
    zones = Internal.getZones(t)
    for zp in zones:
        bcs = getBCMatchs__(zp, zname)
        for b in bcs:
            Internal._rmNode(zp, b)
    return None

def _replaceZoneWithSplit(t, zname, z1, z2):
    zones = Internal.getZones(t)
    for z in zones:
        if z[0] == zname:
            p,c = Internal.getParentOfNode(t, z)
            p[2][c] = z1
            if len(z2) == 4 and z2[3] == 'Zone_t': p[2].append(z2)
            else: p[2] += z2
    return None

def getWinSize(w):
    """Return the size of a window."""
    s1 = w[1]-w[0]+1
    s2 = w[3]-w[2]+1
    s3 = w[5]-w[4]+1
    return s1*s2*s3

def getWinDir(w):
    """Return the direction of a window."""
    (oi1,oi2,oj1,oj2,ok1,ok2) = w
    dim = 3
    if oi1 != oi2 and ok1 == 1 and ok2 == 1: dim = 2
    if oj1 != oj2 and ok1 == 1 and ok2 == 1: dim = 2
    if dim == 2:
        if oi1 == oi2: return 1
        else: return 2
    if dim == 3:
        if oi1 == oi2: return 1
        elif oj1 == oj2: return 2
        else: return 3

def getCutDir(winz1, winz2):
    (oi1,oi2,oj1,oj2,ok1,ok2) = winz1
    (pi1,pi2,pj1,pj2,pk1,pk2) = winz2
    if oi2 == pi1: return 1
    elif oj2 == pj1: return 2
    elif ok2 == pk1: return 3
    elif pi2 == oi1: return 1
    elif pj2 == oj2: return 2
    else: return 3

# Retourne false si les fenetres sont identiques
# Retourne true si wins1 correspond a une sous fenetre de win
def isWinCut(wins1, win):
    if wins1[0] == win[0] and wins1[1] == win[1] and \
            wins1[2] == win[2] and wins1[3] == win[3] and \
            wins1[4] == win[4] and wins1[5] == win[5]:
        return False
    else: return True

# Reporte les BC matchs de z sur z1 et z2 (et modifie t)
def _adaptBCMatch(z, z1, z2, winz1, winz2, t=None):
    bcs = getBCMatchs__(z)
    for b in bcs:
        (oppBlock, winz, winDonor, trirac, periodic) = getBCMatchData__(b)
        triracopp = triracopp__(trirac)

        if oppBlock == z[0]: # self attached BCMatch
            #print("z,z1,z2",z[0],z1[0],z2[0])
            #print('windonor',winDonor, winz1, winz2)
            #print('windir',getWinDir(winz), getCutDir(winz1,winz2))

            #wins1 = intersectWins__(winz1, winDonor, ret=1)
            #wins2 = intersectWins__(winz2, winDonor, ret=1)
            #windonor1 = None; windonor2 = None
            #if wins1 is not None: windonor1 = wins1; oppBlock1 = z1[0]
            #if wins2 is not None: windonor2 = wins2; oppBlock2 = z2[0]

            # Reporte cette BC sur z1
            wini1 = intersectWins__(winz1, winz, ret=0)
            wini = intersectWins__(winz1, winz, ret=1)

            # Point de cut
            if wini is not None: # point de cut sur winz? Pt de cut sur winDonor

                #print('intersect z1 en ',wini1)
                # wini1 sur z, winopp sur z, winopp1 sur z1, winopp2 sur z2
                ind0 = donorIndex__(winz,winDonor,trirac,(wini1[0],wini1[2],wini1[4]))
                ind1 = donorIndex__(winz,winDonor,trirac,(wini1[1],wini1[3],wini1[5]))
                winopp = [min(ind0[0],ind1[0]),max(ind0[0],ind1[0]),
                          min(ind0[1],ind1[1]),max(ind0[1],ind1[1]),min(ind0[2],ind1[2]),max(ind0[2],ind1[2])]
                #print('correspond au donneur ',winopp)

                # DBX
                #winopp1 = intersectWins__(winz1, winopp, ret=0)
                #winopp2 = intersectWins__(winz2, winopp, ret=0)
                #print('donneur intersect z1=',winopp1,'et z2=',winopp2)
                # ENDDBX

                winopp1 = intersectWins__(winz1, winopp, ret=1)
                winopp2 = intersectWins__(winz2, winopp, ret=1)

                if winopp1 is not None and winopp2 is None:
                    if periodic is None: C._addBC2Zone(z1, 'match', 'BCMatch', wini, z1, winopp1, trirac)
                    else: C._addBC2Zone(z1, 'match', 'BCMatch', wini, z1, winopp1, trirac,
                                        rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])

                elif winopp2 is not None and winopp1 is None:
                    if periodic is None: C._addBC2Zone(z1, 'match', 'BCMatch', wini, z2, winopp2, trirac)
                    else: C._addBC2Zone(z1, 'match', 'BCMatch', wini, z2, winopp2, trirac,
                                        rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])

                else: # c'est le bazar, winopp1z sur z, winstart1 sur z
                    winopp1z = intersectWins__(winz1, winopp, ret=0)
                    triracoppn = numpy.array(triracopp)
                    # choper la fenetre inverse de winopp1z
                    ind0 = donorIndex__(winDonor,winz,triracoppn,(winopp1z[0],winopp1z[2],winopp1z[4]))
                    ind1 = donorIndex__(winDonor,winz,triracoppn,(winopp1z[1],winopp1z[3],winopp1z[5]))
                    winstart1 = [min(ind0[0],ind1[0]),max(ind0[0],ind1[0]),
                                 min(ind0[1],ind1[1]),max(ind0[1],ind1[1]),min(ind0[2],ind1[2]),max(ind0[2],ind1[2])]
                    winstarti1 = intersectWins__(winz1, winstart1, ret=1)
                    winstarti2 = intersectWins__(winz2, winstart1, ret=1)
                    if winstarti1 is not None:
                        if periodic is None: C._addBC2Zone(z1, 'match', 'BCMatch', winstarti1, z1, winopp1, trirac)
                        else: C._addBC2Zone(z1, 'match', 'BCMatch', winstarti1, z1, winopp1, trirac,
                                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
                    if winstarti2 is not None:
                        if periodic is None: C._addBC2Zone(z2, 'match', 'BCMatch', winstarti2, z1, winopp1, trirac)
                        else: C._addBC2Zone(z2, 'match', 'BCMatch', winstarti2, z1, winopp1, trirac,
                                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])

                    winopp2z = intersectWins__(winz2, winopp, ret=0)
                    # choper la fenetre inverse de winopp2z
                    ind0 = donorIndex__(winDonor,winz,triracoppn,(winopp2z[0],winopp2z[2],winopp2z[4]))
                    ind1 = donorIndex__(winDonor,winz,triracoppn,(winopp2z[1],winopp2z[3],winopp2z[5]))
                    winstart2 = [min(ind0[0],ind1[0]),max(ind0[0],ind1[0]),
                                 min(ind0[1],ind1[1]),max(ind0[1],ind1[1]),min(ind0[2],ind1[2]),max(ind0[2],ind1[2])]

                    winstarti1 = intersectWins__(winz1, winstart2, ret=1)
                    winstarti2 = intersectWins__(winz2, winstart2, ret=1)
                    if winstarti1 is not None:
                        if periodic is None: C._addBC2Zone(z1, 'match', 'BCMatch', winstarti1, z2, winopp2, trirac)
                        else: C._addBC2Zone(z1, 'match', 'BCMatch', winstarti1, z2, winopp2, trirac,
                                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
                    if winstarti2 is not None:
                        if periodic is None: C._addBC2Zone(z2, 'match', 'BCMatch', winstarti2, z2, winopp2, trirac)
                        else: C._addBC2Zone(z2, 'match', 'BCMatch', winstarti2, z2, winopp2, trirac,
                                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])


            #if wini is not None:
            #    if wins1 is None and wins2 is not None:
            #        oppBlock = z2[0]
            #        winopp = windonor2
            #        if periodic is None: C._addBC2Zone(z1, 'match', 'BCMatch', wini, oppBlock, winopp, trirac)
            #        else: C._addBC2Zone(z1, 'match', 'BCMatch', wini, oppBlock, winopp, trirac,
            #                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
            #    elif wins2 is None and wins1 is not None:
            #        oppBlock = z1[0]
            #        winopp = windonor1
            #        if periodic is None: C._addBC2Zone(z1, 'match', 'BCMatch', wini, oppBlock, winopp, trirac)
            #        else: C._addBC2Zone(z1, 'match', 'BCMatch', wini, oppBlock, winopp, trirac,
            #                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
            #    else: # il y en a deux a creer
            #        if periodic is None: C._addBC2Zone(z1, 'match', 'BCMatch', wini, z1[0], windonor1, trirac)
            #        else: C._addBC2Zone(z1, 'match', 'BCMatch', wini, z1[0], windonor1, trirac,
            #                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
            #        if periodic is None: C._addBC2Zone(z1, 'match', 'BCMatch', wini, z2[0], windonor2, trirac)
            #        else: C._addBC2Zone(z1, 'match', 'BCMatch', wini, z2[0], windonor2, trirac,
            #                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])


                # DBX
                #nodes = Internal.getNodesFromType(z1, 'GridConnectivity1to1_t')
                #for n in nodes:
                #    pr = Internal.getNodeFromName(n, 'PointRange')
                #    prd = Internal.getNodeFromName(n, 'PointRangeDonor')
                #    print('out1',z1[0],n[0],Internal.range2Window(pr[1]))
                #    print('out1',z1[0],n[0],Internal.range2Window(prd[1]))

            # Reporte cette BC sur z2
            wini1 = intersectWins__(winz2, winz, ret=0)
            wini = intersectWins__(winz2, winz, ret=1)

            if wini is not None:

                #print('intersect z2 en ',wini1)
                # wini1 sur z, winopp sur z, winopp1 sur z1, winopp2 sur z2
                ind0 = donorIndex__(winz,winDonor,trirac,(wini1[0],wini1[2],wini1[4]))
                ind1 = donorIndex__(winz,winDonor,trirac,(wini1[1],wini1[3],wini1[5]))
                winopp = [min(ind0[0],ind1[0]),max(ind0[0],ind1[0]),
                          min(ind0[1],ind1[1]),max(ind0[1],ind1[1]),min(ind0[2],ind1[2]),max(ind0[2],ind1[2])]
                #print('correspond au donneur ',winopp)

                # DBX
                #winopp1 = intersectWins__(winz1, winopp, ret=0)
                #winopp2 = intersectWins__(winz2, winopp, ret=0)
                #print('donneur intersect z1=',winopp1,'et z2=',winopp2)
                # ENDDBX

                winopp1 = intersectWins__(winz1, winopp, ret=1)
                winopp2 = intersectWins__(winz2, winopp, ret=1)

                if winopp1 is not None and winopp2 is None:
                    if periodic is None: C._addBC2Zone(z2, 'match', 'BCMatch', wini, z1, winopp1, trirac)
                    else: C._addBC2Zone(z2, 'match', 'BCMatch', wini, z1, winopp1, trirac,
                                        rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])

                elif winopp2 is not None and winopp1 is None:
                    if periodic is None: C._addBC2Zone(z2, 'match', 'BCMatch', wini, z2, winopp2, trirac)
                    else: C._addBC2Zone(z2, 'match', 'BCMatch', wini, z2, winopp2, trirac,
                                        rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])

                else: # c'est le bazar, winopp1z sur z, winstart1 sur z
                    triracoppn = numpy.array(triracopp)
                    winopp1z = intersectWins__(winz1, winopp, ret=0)
                    # choper la fenetre inverse de winopp1z
                    ind0 = donorIndex__(winDonor,winz,triracoppn,(winopp1z[0],winopp1z[2],winopp1z[4]))
                    ind1 = donorIndex__(winDonor,winz,triracoppn,(winopp1z[1],winopp1z[3],winopp1z[5]))
                    winstart1 = [min(ind0[0],ind1[0]),max(ind0[0],ind1[0]),
                                 min(ind0[1],ind1[1]),max(ind0[1],ind1[1]),min(ind0[2],ind1[2]),max(ind0[2],ind1[2])]

                    # DBX
                    #winstarti1 = intersectWins__(winz1, winstart1, ret=0)
                    #winstarti2 = intersectWins__(winz2, winstart1, ret=0)
                    #print("correspond au start1 ",winstarti1,winstarti2)
                    # ENDDBX

                    winstarti1 = intersectWins__(winz1, winstart1, ret=1)
                    winstarti2 = intersectWins__(winz2, winstart1, ret=1)
                    if winstarti1 is not None:
                        if periodic is None: C._addBC2Zone(z1, 'match', 'BCMatch', winstarti1, z1, winopp1, trirac)
                        else: C._addBC2Zone(z1, 'match', 'BCMatch', winstarti1, z1, winopp1, trirac,
                                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
                    if winstarti2 is not None:
                        if periodic is None: C._addBC2Zone(z2, 'match', 'BCMatch', winstarti2, z1, winopp1, trirac)
                        else: C._addBC2Zone(z2, 'match', 'BCMatch', winstarti2, z1, winopp1, trirac,
                                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])

                    winopp2z = intersectWins__(winz2, winopp, ret=0)
                    # choper la fenetre inverse de winopp2z
                    ind0 = donorIndex__(winDonor,winz,triracoppn,(winopp2z[0],winopp2z[2],winopp2z[4]))
                    ind1 = donorIndex__(winDonor,winz,triracoppn,(winopp2z[1],winopp2z[3],winopp2z[5]))
                    winstart2 = [min(ind0[0],ind1[0]),max(ind0[0],ind1[0]),
                                 min(ind0[1],ind1[1]),max(ind0[1],ind1[1]),min(ind0[2],ind1[2]),max(ind0[2],ind1[2])]

                    # DBX
                    #winstarti1 = intersectWins__(winz1, winstart2, ret=0)
                    #winstarti2 = intersectWins__(winz2, winstart2, ret=0)
                    #print("correspond au start2 ",winstarti1,winstarti2)
                    # ENDDBX

                    winstarti1 = intersectWins__(winz1, winstart2, ret=1)
                    winstarti2 = intersectWins__(winz2, winstart2, ret=1)

                    if winstarti1 is not None:
                        if periodic is None: C._addBC2Zone(z1, 'match', 'BCMatch', winstarti1, z2, winopp2, trirac)
                        else: C._addBC2Zone(z1, 'match', 'BCMatch', winstarti1, z2, winopp2, trirac,
                                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
                    if winstarti2 is not None:
                        if periodic is None: C._addBC2Zone(z2, 'match', 'BCMatch', winstarti2, z2, winopp2, trirac)
                        else: C._addBC2Zone(z2, 'match', 'BCMatch', winstarti2, z2, winopp2, trirac,
                                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])


            #    if wins1 is None and wins2 is not None:
            #        oppBlock = z2[0]
            #        winopp = windonor2
            #        if periodic is None: C._addBC2Zone(z2, 'match', 'BCMatch', wini, oppBlock, winopp, trirac)
            #        else: C._addBC2Zone(z2, 'match', 'BCMatch', wini, oppBlock, winopp, trirac,
            #                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
            #    elif wins2 is None and wins1 is not None:
            #        oppBlock = z1[0]
            #        winopp = windonor1
            #        if periodic is None: C._addBC2Zone(z2, 'match', 'BCMatch', wini, oppBlock, winopp, trirac)
            #        else: C._addBC2Zone(z2, 'match', 'BCMatch', wini, oppBlock, winopp, trirac,
            #                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
            #    else:
            #        if periodic is None: C._addBC2Zone(z2, 'match', 'BCMatch', wini, z1[0], windonor1, trirac)
            #        else: C._addBC2Zone(z2, 'match', 'BCMatch', wini, z1[0], windonor1, trirac,
            #                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
            #        if periodic is None: C._addBC2Zone(z2, 'match', 'BCMatch', wini, z2[0], windonor2, trirac)
            #        else: C._addBC2Zone(z2, 'match', 'BCMatch', wini, z2[0], windonor2, trirac,
            #                            rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])

                # DBX
                #nodes = Internal.getNodesFromType(z2, 'GridConnectivity1to1_t')
                #for n in nodes:
                #    pr = Internal.getNodeFromName(n, 'PointRange')
                #    prd = Internal.getNodeFromName(n, 'PointRangeDonor')
                #    print('out2',z2[0],n[0],Internal.range2Window(pr[1]))
                #    print('out2',z2[0],n[0],Internal.range2Window(prd[1]))

        else: # not self opposite

            # Reporte cette BC sur z1
            wini1 = intersectWins__(winz1, winz, ret=0)
            wini = intersectWins__(winz1, winz, ret=1)

            if wini is not None:
                ind0 = donorIndex__(winz,winDonor,trirac,(wini1[0],wini1[2],wini1[4]))
                ind1 = donorIndex__(winz,winDonor,trirac,(wini1[1],wini1[3],wini1[5]))
                winopp = [min(ind0[0],ind1[0]),max(ind0[0],ind1[0]),
                          min(ind0[1],ind1[1]),max(ind0[1],ind1[1]),min(ind0[2],ind1[2]),max(ind0[2],ind1[2])]

                if periodic is None: C._addBC2Zone(z1, 'match', 'BCMatch', wini, oppBlock, winopp, trirac)
                else: C._addBC2Zone(z1, 'match', 'BCMatch', wini, oppBlock, winopp, trirac,
                                    rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
                if t is not None:
                    zopp = Internal.getNodeFromName2(t, oppBlock)
                    if zopp is not None:
                        if periodic is None: C._addBC2Zone(zopp, 'match', 'BCMatch', winopp, z1[0], wini, triracopp)
                        else: C._addBC2Zone(zopp, 'match', 'BCMatch', winopp, z1[0], wini, triracopp,
                                            rotationCenter=periodic[1], rotationAngle=-1.*periodic[2], translation=-1.*periodic[0])
                        # DBX
                        #nodes = Internal.getNodesFromType(z2, 'GridConnectivity1to1_t')
                        #for n in nodes:
                        #    pr = Internal.getNodeFromName(n, 'PointRange')
                        #    prd = Internal.getNodeFromName(n, 'PointRangeDonor')
                        #    print('out3',z2[0],n[0],Internal.range2Window(pr[1]))
                        #    print('out3',z2[0],n[0],Internal.range2Window(prd[1]))

            # Reporte cette BC sur z2
            wini1 = intersectWins__(winz2, winz, ret=0)
            wini = intersectWins__(winz2, winz, ret=1)

            if wini is not None:
                ind0 = donorIndex__(winz,winDonor,trirac,(wini1[0],wini1[2],wini1[4]))
                ind1 = donorIndex__(winz,winDonor,trirac,(wini1[1],wini1[3],wini1[5]))
                winopp = [min(ind0[0],ind1[0]),max(ind0[0],ind1[0]),
                          min(ind0[1],ind1[1]),max(ind0[1],ind1[1]),min(ind0[2],ind1[2]),max(ind0[2],ind1[2])]
                if periodic is None: C._addBC2Zone(z2, 'match', 'BCMatch', wini, oppBlock, winopp, trirac)
                else: C._addBC2Zone(z2, 'match', 'BCMatch', wini, oppBlock, winopp, trirac,
                                    rotationCenter=periodic[1], rotationAngle=periodic[2], translation=periodic[0])
                if t is not None:
                    zopp = Internal.getNodeFromName2(t, oppBlock)
                    if zopp is not None:
                        if periodic is None: C._addBC2Zone(zopp, 'match', 'BCMatch', winopp, z2[0], wini, triracopp)
                        else: C._addBC2Zone(zopp, 'match', 'BCMatch', winopp, z2[0], wini, triracopp,
                                            rotationCenter=periodic[1], rotationAngle=-1.*periodic[2], translation=-1.*periodic[0])
                        # DBX
                        #nodes = Internal.getNodesFromType(z2, 'GridConnectivity1to1_t')
                        #for n in nodes:
                        #    pr = Internal.getNodeFromName(n, 'PointRange')
                        #    prd = Internal.getNodeFromName(n, 'PointRangeDonor')
                        #    print('out4',z2[0],n[0],Internal.range2Window(pr[1]))
                        #    print('out4',z2[0],n[0],Internal.range2Window(prd[1]))


    # Enleve les raccords de qui referent z[0] dans t
    if t is not None: _deleteBCMatchRef(t, z[0])
    return None

# split (z)
# Split une zone suivant une direction et un index
# z: zone ou zone skeleton
# Retourne deux zones avec les BCs + BCMatchs OK
# Si t est donne, t est modifie
def split(z, dir=1, index=1, t=None):
    """Split a structured zone given index and direction."""
    dimz = Internal.getZoneDim(z)
    ni = dimz[1]; nj = dimz[2]; nk = dimz[3]
    zoneName = z[0]

    # force loc2glob on z
    src, loc2glob = Internal.getLoc2Glob(z)
    if src is None:
        Internal._setLoc2Glob(z, z[0], win=[1,ni,1,nj,1,nk], sourceDim=[ni,nj,nk])

    if dir == 1: # direction i
        z1 = subzone(z, (1,1,1), (index,-1,-1))
        z1[0] = C.getZoneName(zoneName)
        z2 = subzone(z, (index,1,1), (-1,-1,-1))
        z2[0] = C.getZoneName(zoneName)
        z1[1] = Internal.array2PyTreeDim([None,None,index,nj,nk])
        z2[1] = Internal.array2PyTreeDim([None,None,ni-index+1,nj,nk])
        #z1[0] = z[0]+'A'; z2[0] = z[0]+'B'
        winz1 = [1,index,1,nj,1,nk]
        winz2 = [index,ni,1,nj,1,nk]
        w1 = winz1; dim = [ni,nj,nk]
        w2 = winz2; source = z[0]
        src, loc2glob = Internal.getLoc2Glob(z)
        if loc2glob is not None:
            ws = list(loc2glob[0:6]); dim = list(loc2glob[6:9])
            w1 = composeWins__(w1, ws)
            w2 = composeWins__(w2, ws)
            source = src
        Internal._setLoc2Glob(z1, source, win=w1, sourceDim=dim)
        Internal._setLoc2Glob(z2, source, win=w2, sourceDim=dim)

    elif dir == 2: # direction j
        z1 = subzone(z, (1,1,1), (-1,index,-1))
        z1[0] = C.getZoneName(zoneName)
        z2 = subzone(z, (1,index,1), (-1,-1,-1))
        z2[0] = C.getZoneName(zoneName)
        z1[1] = Internal.array2PyTreeDim([None,None,ni,index,nk])
        z2[1] = Internal.array2PyTreeDim([None,None,ni,nj-index+1,nk])
        #z1[0] = z[0]+'A'; z2[0] = z[0]+'B'
        winz1 = [1,ni,1,index,1,nk]
        winz2 = [1,ni,index,nj,1,nk]
        w1 = winz1; dim = [ni,nj,nk]
        w2 = winz2; source = z[0]
        src, loc2glob = Internal.getLoc2Glob(z)
        if loc2glob is not None:
            ws = list(loc2glob[0:6]); dim = list(loc2glob[6:9])
            w1 = composeWins__(w1, ws)
            w2 = composeWins__(w2, ws)
            source = src
        Internal._setLoc2Glob(z1, source, win=w1, sourceDim=dim)
        Internal._setLoc2Glob(z2, source, win=w2, sourceDim=dim)

    else: # direction k
        z1 = subzone(z, (1,1,1), (-1,-1,index))
        z1[0] = C.getZoneName(zoneName)
        z2 = subzone(z, (1,1,index), (-1,-1,-1))
        z2[0] = C.getZoneName(zoneName)
        z1[1] = Internal.array2PyTreeDim([None,None,ni,nj,index]) # force dim (skel)
        z2[1] = Internal.array2PyTreeDim([None,None,ni,nj,nk-index+1])
        #z1[0] = z[0]+'A'; z2[0] = z[0]+'B'
        winz1 = [1,ni,1,nj,1,index]
        winz2 = [1,ni,1,nj,index,nk]
        w1 = winz1; dim = [ni,nj,nk]
        w2 = winz2; source = z[0]
        src, loc2glob = Internal.getLoc2Glob(z)
        if loc2glob is not None:
            ws = list(loc2glob[0:6]); dim = list(loc2glob[6:9])
            w1 = composeWins__(w1, ws)
            w2 = composeWins__(w2, ws)
            source = src
        Internal._setLoc2Glob(z1, source, win=w1, sourceDim=dim)
        Internal._setLoc2Glob(z2, source, win=w2, sourceDim=dim)

    _createInternalBCMatch(z1, z2, dir)
    _adaptBCMatch(z, z1, z2, winz1, winz2, t)
    if t is not None: _replaceZoneWithSplit(t, z[0], z1, z2)
    return z1, z2

#=====================================================================
# IN: window: [imin,imax,jmin,jmax,kmin,kmax]
# Retourne 1, 2, 3
#=====================================================================
def getWinDir2__(win):
    if win[0] == win[1]: return 1
    if win[2] == win[3]: return 2
    if win[4] == win[5]: return 3
    return -1

#=====================================================================
# Casse les zones structurees pour que les raccords soient
# sur des faces pleines
#=====================================================================
def splitFullMatch(t):
    """Split all zones for matching on full faces."""
    tp = Internal.copyRef(t)
    _splitFullMatch(tp)
    return tp

def _splitFullMatch(t):
    """Split all zones for matching on full faces."""
    zones = Internal.getZones(t)
    reflen = len(zones)
    newlen = reflen
    while reflen < newlen:
        reflen = newlen
        stack = []
        for z in zones:
            dim = Internal.getZoneDim(z)
            if dim[0] == 'Structured': stack.append(z)

        while len(stack) > 0:
            z = stack.pop(0)
            splitFullMatch__(z, stack, t)

        zones = Internal.getZones(t)
        newlen = len(zones)

    return None

def splitFullMatch__(z, stack, t):
    dim = Internal.getZoneDim(z)
    ni = dim[1]; nj = dim[2]; nk = dim[3]
    bcs = getBCMatchs__(z)
    for bc in bcs:
        (oppBlock, rnge, donor, trf, periodic) = getBCMatchData__(bc)
        dir = getWinDir2__(rnge)

        # verifie si la fenetre est full
        if (dir == 2 or dir == 3) and rnge[0] > 1:
            z1,z2 = split(z, 1, rnge[0], t)
            stack.append(z1); stack.append(z2)
            return

        if (dir == 2 or dir == 3) and rnge[1] < ni:
            z1,z2 = split(z, 1, rnge[1], t)
            stack.append(z1); stack.append(z2)
            return

        if (dir == 1 or dir == 3) and rnge[2] > 1: # need split, stack
            z1,z2 = split(z, 2, rnge[2], t)
            stack.append(z1); stack.append(z2)
            return

        if (dir == 1 or dir == 3) and rnge[3] < nj: # need split, stack
            z1,z2 = split(z, 2, rnge[3], t)
            stack.append(z1); stack.append(z2)
            return

        if (dir == 1 or dir == 2) and rnge[4] > 1: # need split, stack
            z1,z2 = split(z, 3, rnge[4], t)
            stack.append(z1); stack.append(z2)
            return

        if (dir == 1 or dir == 2) and rnge[5] < nk: # need split, stack
            z1,z2 = split(z, 3, rnge[5], t)
            stack.append(z1); stack.append(z2)
            return

#======================================
# reorder the numerotation of the mesh
#======================================
# reorder BCs
def _reorderBC__(t, order):
    nodes = Internal.getZones(t)
    for z in nodes:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured' and len(order) == 3:
            oi = order[0]; oj = order[1]; ok = order[2]
            wins = Internal.getNodesFromType(z, 'BC_t')
            for w in wins:
                w0 = Internal.getNodeFromName1(w, 'PointRange')
                range0 = reorderIndices__(w0[1], dim, oi, oj, ok)
                r2 = Internal.window2Range(range0)
                w0[1] = r2
    return None

# Reorder du range d'une BCOverlap dans le cas ou
# la zone contenant ce BCOverlap est reordonnee
def _reorderBCOverlap__(a, order):
    nodes = Internal.getZones(a)
    for z in nodes:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured' and len(order) == 3:
            oi = order[0]; oj = order[1]; ok = order[2]
            connect = Internal.getNodesFromType(z, 'GridConnectivity_t')
            for i in connect:
                pr = Internal.getNodeFromName1(i, 'PointRange')
                r = Internal.getNodesFromType(i, 'GridConnectivityType_t')
                if r != []:
                    val = Internal.getValue(r[0])
                    if val == 'Overset':
                        range0 = reorderIndices__(pr[1], dim, oi, oj, ok)
                        r2 = Internal.window2Range(range0)
                        pr[1] = r2
    return None

def reorderIndices__(w0, dim, oi, oj, ok):
    r = Internal.range2Window(w0)
    i1 = r[0]; i2 = r[1]; j1 = r[2]; j2 = r[3]; k1 = r[4]; k2 = r[5]
    #i1 = w0[0,0]; j1 = w0[1,0]; k1 = w0[2,0]
    #i2 = w0[0,1]; j2 = w0[1,1]; k2 = w0[2,1]
    i1o = i1; j1o = j1; k1o = k1
    i2o = i2; j2o = j2; k2o = k2
    # oi
    if   oi == 2: j1o = i1; j2o = i2
    elif oi == 3: k1o = i1; k2o = i2
    elif oi == -1:
        if i1 == i2:
            if i1 == 1: i1o = dim[1]; i2o = dim[1]
            else: i1o = 1; i2o = 1
        else:
            i1o = dim[1]-i2+1
            i2o = dim[1]-i1+1

    elif oi == -2:
        j1o = i2; j2o = i1
        if i1 == i2:
            if i1 == 1: j1o = dim[1]; j2o = dim[1]
            else: j1o = 1; j2o = 1
        else:
            j1o = dim[1]-i2+1
            j2o = dim[1]-i1+1
    elif oi == -3:
        k1o = i2; k2o = i1
        if i1 == i2:
            if i1 == 1: k1o = dim[1]; k2o = dim[1]
            else: k1o = 1; k2o = 1
        else:
            k1o = dim[1]-i2+1
            k2o = dim[1]-i1+1
    # oj
    if   oj == 1: i1o = j1; i2o = j2
    elif oj == 3: k1o = j1; k2o = j2
    elif oj == -1:
        i1o = j2; i2o = j1
        if j1 == j2:
            if j1 == 1: i1o = dim[2]; i2o = dim[2]
            else: i1o = 1; i2o = 1
        else:
            i1o = dim[2]-j2+1
            i2o = dim[2]-j1+1
    elif oj == -2:
        if j1 == j2:
            if j1 == 1: j1o = dim[2]; j2o = dim[2]
            else: j1o = 1; j2o = 1
        else:
            j1o = dim[2]-j2+1
            j2o = dim[2]-j1+1
    elif oj == -3:
        k1o = j2; k2o = j1
        if j1 == j2:
            if j1 == 1: k1o = dim[2]; k2o = dim[2]
            else: k1o = 1; k2o = 1
        else:
            k1o = dim[2]-j2+1
            k2o = dim[2]-j1+1
    # ok
    if   ok== 1: i1o = k1; i2o = k2
    elif ok== 2: j1o = k1; j2o = k2
    elif ok==-1:
        i1o = k2; i2o = k1
        if k1 == k2:
            if k1 == 1: i1o = dim[3]; i2o = dim[3]
            else: i1o = 1; i2o = 1
        else:
            i1o = dim[3]-k2+1
            i2o = dim[3]-k1+1
    elif ok==-2:
        j1o = k2; j2o = k1
        if k1 == k2:
            if k1 == 1: j1o = dim[3]; j2o = dim[3]
            else: j1o = 1; j2o = 1
        else:
            j1o = dim[3]-k2+1
            j2o = dim[3]-k1+1
    elif ok==-3:
        if k1 == k2:
            if k1 == 1: k1o = dim[3]; k2o = dim[3]
            else: k1o = 1; k2o = 1
        else:
            k1o = dim[3]-k2+1
            k2o = dim[3]-k1+1

    return [min(i1o,i2o),max(i1o,i2o),min(j1o,j2o),max(j1o,j2o),min(k1o,k2o),max(k1o,k2o)]

#=============================================================================
# Update les indices du noeud Transform pour reordonner le trirac
#=============================================================================
def reorderTrirac__(transfo, order=[1,2,3]):
    nt = transfo.size
    TM = numpy.zeros((nt,nt), dtype=Internal.E_NpyInt) # matrice liee a la transfo
    OM = numpy.zeros((nt,nt), dtype=Internal.E_NpyInt) # matrice transposee de order:  son inverse
    BM = numpy.zeros((nt,nt), dtype=Internal.E_NpyInt) # T.R^-1 pour reordonner la transfo
    from numpy import linalg
    for i in range(nt):
        for j in range(nt):
            a = transfo[j]
            if abs(a) == i+1:
                if a > 0: TM[i,j] = 1
                else: TM[i,j] =-1
            b = order[j]
            OM[i,j] = 0
            if abs(b) == i+1:
                if b > 0: OM[i,j] = 1
                else: OM[i,j] =-1
    OM = linalg.inv(OM)
    BM = numpy.dot(TM,OM)
    trirac = [1,2,3]
    if nt == 2: trirac = [1,2]
    for i in range(nt):
        for j in range(nt):
            if BM[j,i] > 0: trirac[i] = j+1
            elif BM[j,i] < 0: trirac[i] = -(j+1)
    for i in range(nt): transfo[i] = trirac[i]
    return transfo

def reorderTriracOpp__(transfo, order=[1,2,3]):
    nt = transfo.size
    TM = numpy.zeros((nt,nt), dtype=Internal.E_NpyInt) # matrice liee a la transfo opp de z1 vers z1
    OM = numpy.zeros((nt,nt), dtype=Internal.E_NpyInt) # matrice transposee de order appliquee a z1:  son inverse
    BM = numpy.zeros((nt,nt), dtype=Internal.E_NpyInt) # R.T pour reordonner la transfo
    from numpy import linalg
    for i in range(nt):
        for j in range(nt):
            a = transfo[j]
            if abs(a) == i+1:
                if a > 0: TM[i,j] = 1
                else: TM[i,j] =-1
            b = order[j]
            OM[i,j] = 0
            if abs(b) == i+1:
                if b > 0: OM[i,j] = 1
                else: OM[i,j] =-1
    BM = numpy.dot(OM,TM)
    trirac = [1,2,3]
    if nt == 2: trirac = [1,2]
    for i in range(nt):
        for j in range(nt):
            if BM[j,i] > 0: trirac[i] = j+1
            elif BM[j,i] < 0: trirac[i] = -(j+1)
    for i in range(nt): transfo[i] = trirac[i]
    return transfo

# Reorder in a the PointRange for zones and the PointRangeDonor for donor zones of name in zoneNames
def _reorderBCNearMatch__(a, order, zoneNames):
    nodes = Internal.getZones(a)
    for z in nodes:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured' and len(order) == 3:
            oi = order[0]; oj = order[1]; ok = order[2]
            connect = Internal.getNodesFromType(z, 'GridConnectivity_t')
            for cn in connect:
                type = Internal.getNodeFromName1(cn, 'GridConnectivityType')
                if type is not None: val = Internal.getValue(type)
                else: val = 'Unknown'
                if val == 'Abutting':
                    # modif du PointRange si z dans la liste
                    for name in zoneNames:
                        if name == z[0]:# reorder le PointRange
                            wins = Internal.getNodesFromName1(cn, 'PointRange')
                            for w in wins:
                                (parent, d) = Internal.getParentOfNode(cn, w)
                                range0 = reorderIndices__(w[1],dim,oi,oj,ok)
                                r2 = Internal.window2Range(range0)
                                parent[2][d][1] = r2
                            transfos = Internal.getNodesFromName2(cn, 'Transform')
                            for transfo in transfos:
                                (parent, d) = Internal.getParentOfNode(cn, transfo)
                                trirac = reorderTrirac__(transfo[1],order)
                                parent[2][d][1] = trirac
                            nmratios = Internal.getNodesFromName2(cn, 'NMRatio')
                            for nmratio in nmratios:
                                (parent,d) = Internal.getParentOfNode(cn, nmratio)
                                triracnm = numpy.zeros((3), dtype=Internal.E_NpyInt)
                                for i in range(3): triracnm[i] = i+1
                                triracnm = reorderTrirac__(triracnm,[abs(oi),abs(oj),abs(ok)])
                                nmr = nmratio[1]
                                for i in range(3): nmr[i] = nmr[triracnm[i]-1]
                            break

                    # modif du PointRangeDonor si la donor est dans la liste
                    donorName = Internal.getValue(cn)
                    for name in zoneNames:
                        if name == donorName: # reorder le PointRangeDonor
                            zdnr = Internal.getNodesFromName(a, donorName)
                            if zdnr != []:
                                zdnr = zdnr[0]
                                dimDnr = Internal.getZoneDim(zdnr)
                                wins = Internal.getNodesFromName2(cn, 'PointRangeDonor')
                                for w in wins:
                                    (parent, d) = Internal.getParentOfNode(cn, w)
                                    range0 = reorderIndices__(w[1],dimDnr,oi,oj,ok)
                                    r2 = Internal.window2Range(range0)
                                    parent[2][d][1] = r2

                                transfos = Internal.getNodesFromName2(cn, 'Transform')
                                for transfo in transfos:
                                    (parent, d) = Internal.getParentOfNode(cn, transfo)
                                    trirac=reorderTriracOpp__(transfo[1],order)
                                    parent[2][d][1] = trirac
                                break

    return None
# Reorder in a the PointRange for zones and the PointRangeDonor for
# donor zones of name in zoneNames
def _reorderBCMatch__(a, order, zoneNames):
    nodes = Internal.getZones(a)
    for z in nodes:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured' and len(order) == 3:
            oi = order[0]; oj = order[1]; ok = order[2]
            connect = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')

            for cn in connect:
                # modif du PointRange si z dans la liste
                for name in zoneNames:
                    if name == z[0]: # reorder le PointRange
                        wins = Internal.getNodesFromName1(cn, 'PointRange')
                        for w in wins:
                            (parent, d) = Internal.getParentOfNode(cn, w)
                            range0 = reorderIndices__(w[1],dim,oi,oj,ok)
                            r2 = Internal.window2Range(range0)
                            parent[2][d][1] = r2
                        transfos = Internal.getNodesFromName1(cn, 'Transform')
                        for transfo in transfos:
                            (parent, d) = Internal.getParentOfNode(cn, transfo)
                            trirac = reorderTrirac__(transfo[1],order)
                            parent[2][d][1] = trirac
                        break

                # modif du PointRangeDonor si la donor est dans la liste
                donorName = Internal.getValue(cn)
                for name in zoneNames:
                    if name == donorName: # reorder le PointRangeDonor
                        zdnr = Internal.getNodeFromName(a, donorName)
                        if zdnr is not None:
                            dimDnr = Internal.getZoneDim(zdnr)
                            wins = Internal.getNodesFromName1(cn, 'PointRangeDonor')
                            for w in wins:
                                (parent, d) = Internal.getParentOfNode(cn, w)
                                range0 = reorderIndices__(w[1],dimDnr,oi,oj,ok)
                                r2 = Internal.window2Range(range0)
                                parent[2][d][1] = r2

                            transfos = Internal.getNodesFromName1(cn, 'Transform')
                            for transfo in transfos:
                                (parent, d) = Internal.getParentOfNode(cn, transfo)
                                trirac = reorderTriracOpp__(transfo[1],order)
                                parent[2][d][1] = trirac
                            break
    return None

# reorder GridConnectivity
def _reorderGC__(t, order):
    # GridConnectivity
    nodes = Internal.getZones(t)
    for z in nodes:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured' and len(order) == 3:
            oi = order[0]; oj = order[1]; ok = order[2]
            connect = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
            for cn in connect:
                wins = Internal.getNodesFromName2(cn, 'PointRange')
                for w in wins:
                    (parent, d) = Internal.getParentOfNode(cn, w)
                    range0 = reorderIndices__(w[1],dim,oi,oj,ok)
                    r2 = Internal.window2Range(range0)
                    parent[2][d][1] = r2
    return None

# Reorder the numbering of t. t and toptree are modified
def reorder(t, order, topTree=[]):
    """Reorder the numerotation of mesh.
    Usage: reorder(a, (2,1,-3))"""
    a = Internal.copyRef(t)
    _reorder(a, order, topTree)
    return a

# If topTree is given, must be used in place!
def _reorder(t, order, topTree=[]):
    if len(order) == 3: _reorderStruct__(t, order, topTree)
    else: _reorderUnstruct__(t, order)
    return None

def _reorderStruct__(t, order, topTree):
    loc = 'both'
    istoptree = Internal.isTopTree(t)
    if istoptree:
        _reorderBC__(t, order)
        _reorderBCOverlap__(t, order)
        zones = Internal.getZones(t)
        zoneNames=[] # zones dont le PointRange est a modifier ou en tant que PointRangeDonor si c'est la zoneDonor
        for z in zones: zoneNames.append(z[0])
        _reorderBCMatch__(t, order, zoneNames)
        _reorderBCNearMatch__(t, order, zoneNames)
        C._TZA1(t, loc, loc, True, Transform.reorder, order)
        return None
    else:
        istoptree = Internal.isTopTree(topTree)
        if not istoptree: # pas de modif des zoneDonors dans les BCMatch !
            _reorderBC__(t, order)
            _reorderGC__(t, order)
            C._TZA1(t, loc, loc, True, Transform.reorder, order)
            return None
        else: # toptree fourni
            _reorderBC__(t, order)
            _reorderBCOverlap__(t, order)
            C._TZA1(t, loc, loc, True, Transform.reorder, order)
            zones = Internal.getZones(t)
            zoneNames = []
            for z in zones: zoneNames.append(z[0])
            _reorderBCMatch__(topTree, order, zoneNames)
            _reorderBCNearMatch__(topTree, order, zoneNames)
            return None

def _reorderUnstruct__(t, order):
    C._TZA1(t, 'nodes', 'nodes', True, Transform.reorder, order)
    return None

#=============================================================
# reorder all zones to get the same orientation of the normals
#=============================================================
def reorderAll(t, dir=1):
    """Orientate normals of all surface blocks consistently in one direction (1) or the opposite (-1).
    For unstructured inputs, when dir is set to 1(-1), it means outward(inward).
    Usage: reorderAll(arrays, dir)"""
    tp = Internal.copyRef(t)
    _reorderAll(tp, dir)
    return tp

def _reorderAll(t, dir=1):
    """Orientate normals of all surface blocks consistently in one direction (1) or the opposite (-1).
    For unstructured inputs, when dir is set to 1(-1), it means outward(inward).
    Usage: reorderAll(arrays, dir)"""
    C._fillMissingVariables(t)
    zones = Internal.getZones(t)
    allBCInfos = C.extractBCInfo(zones)
    C._deleteZoneBC__(zones)
    C._deleteGridConnectivity__(zones)
    coords = []; fn = []; fc = []; indirn = []; indirc = []
    nofan=0; nofac=0
    for z in zones:
        # champs en noeuds
        coords.append(C.getFields(Internal.__GridCoordinates__, z)[0])
        fnz = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
        if fnz == []: indirn.append(0)
        else: indirn.append(1); nofan += 1
        fcz = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
        if fcz == []: indirc.append(0)
        else: indirc.append(1); nofac += 1
        fn.append(fnz); fc.append(fcz)

    if nofan == 0 and nofac == 0: # pas de champs en centres et noeuds
        an = Transform.reorderAll(coords, dir)
        C.setFields(an, zones, 'nodes')
    elif nofan == len(zones) and nofac == 0:
        an = Converter.addVars([coords,fn])
        an = Transform.reorderAll(an, dir)
        C.setFields(an,zones,'nodes')
    elif nofan == 0 and nofac == len(zones):
        fc = Converter.center2Node(fc)
        ac = Converter.addVars([coords,fc])
        ac = Transform.reorderAll(ac, dir)
        ac = Converter.rmVars(ac,['CoordinateX','CoordinateY','CoordinateZ'])
        ac = Converter.node2Center(ac)
        C.setFields(ac,zones,'centers')
    elif nofan == len(zones) and nofac == len(zones):
        an = Converter.addVars([coords,fn])
        an = Transform.reorderAll(an, dir)
        C.setFields(an, zones, 'nodes')

        fc = Converter.center2Node(fc)
        ac = Converter.addVars([coords,fc])
        ac = Transform.reorderAll(ac, dir)
        ac = Converter.rmVars(ac,['CoordinateX','CoordinateY','CoordinateZ'])
        ac = Converter.node2Center(ac)
        C.setFields(ac, zones, 'centers')

    # This function is in place
    zones = C.identifyBC(zones, allBCInfos)
    return None

#=============================================================================
# Align I,J,K directions of a Cartesian mesh with X,Y,Z axes
# Order the mesh such that it is direct
#=============================================================================
def makeCartesianXYZ(t, tol=1.e-10):
    t2 = Internal.copyRef(t)
    _makeCartesianXYZ(t2, tol)
    return t2

def _makeCartesianXYZ(t, tol=1.e-10):
    for z in Internal.getZones(t):
        dims = Internal.getZoneDim(z)
        if dims[0] == 'Structured':
            ni = dims[1]; nj = dims[2]; nk = dims[3]
            ind = 0; indi = ind+1; indj = ind+ni; indk = ind+ni*nj
            dx_i = C.getValue(z,'CoordinateX',indi)-C.getValue(z,'CoordinateX',ind)
            dy_i = C.getValue(z,'CoordinateX',indj)-C.getValue(z,'CoordinateX',ind)
            dz_i = C.getValue(z,'CoordinateX',indk)-C.getValue(z,'CoordinateX',ind)
            diri = 1; dirj = 2; dirk = 3
            diri = -1; dirj = -1; dirk = -1
            if abs(dx_i) > tol: diri = 1
            elif abs(dy_i) > tol: diri = 2
            elif abs(dz_i) > tol: diri = 3
            dx_j = C.getValue(z,'CoordinateY',indi)-C.getValue(z,'CoordinateY',ind)
            dy_j = C.getValue(z,'CoordinateY',indj)-C.getValue(z,'CoordinateY',ind)
            dz_j = C.getValue(z,'CoordinateY',indk)-C.getValue(z,'CoordinateY',ind)
            if abs(dx_j) > tol: dirj = 1
            elif abs(dy_j) > tol: dirj = 2
            elif abs(dz_j) > tol: dirj = 3
            dx_k = C.getValue(z,'CoordinateZ',indi)-C.getValue(z,'CoordinateZ',ind)
            dy_k = C.getValue(z,'CoordinateZ',indj)-C.getValue(z,'CoordinateZ',ind)
            dz_k = C.getValue(z,'CoordinateZ',indk)-C.getValue(z,'CoordinateZ',ind)
            if abs(dx_k) > tol: dirk = 1
            elif abs(dy_k) > tol: dirk = 2
            elif abs(dz_k) > tol: dirk = 3

            if diri==-1:
                sump = abs(dirj)+abs(dirk)
                if sump==3: diri = 3
                elif sump==4: diri = 2
                elif sump == 5: diri=1
            if dirj==-1:
                sump = abs(diri)+abs(dirk)
                if sump==3: dirj = 3
                elif sump==4: dirj = 2
                elif sump == 5: dirj=1
            if dirk==-1:
                sump = abs(dirj)+abs(diri)
                if sump==3: dirk = 3
                elif sump==4: dirk = 2
                elif sump == 5: dirk=1

            dirs = [0,0,0]
            if diri == 1: dirs[0] = 1
            elif diri==2: dirs[1] = 1
            else: dirs[2] = 1
            if dirj == 1: dirs[0] = 2
            elif dirj==2: dirs[1] = 2
            else: dirs[2] = 2
            if dirk == 1: dirs[0] = 3
            elif dirk==2: dirs[1] = 3
            else: dirs[2] = 3

            _reorder(z,(dirs[0], dirs[1], dirs[2]))
            dims = Internal.getZoneDim(z)
            ni = dims[1]; nj = dims[2]; nk = dims[3]
            ind = 0; indi = ind+1; indj = ind+ni; indk = ind+ni*nj
            diri = 1; dirj = 1; dirk = 1
            dx_i = C.getValue(z,'CoordinateX',indi)-C.getValue(z,'CoordinateX',ind)
            ok = 0
            diri = 1; dirj = 2; dirk = 3
            if dx_i < tol: diri =-1; ok = 1
            dy_j = C.getValue(z,'CoordinateY',indj)-C.getValue(z,'CoordinateY',ind)
            if dy_j < tol: dirj =-2; ok = 1
            dz_k = C.getValue(z,'CoordinateZ',indk)-C.getValue(z,'CoordinateZ',ind)
            if dz_k < tol: dirk =-3; ok = 1
            if ok == 1: _reorder(z,(diri,dirj,dirk))
    return None

def makeDirect(t):
    """Make a structured grid direct.
    Usage: makeDirect(t)"""
    tp = Internal.copyRef(t)
    _makeDirect(tp)
    return tp

def _makeDirect(t):
    """Make a structured grid direct."""
    import KCore.Vector as Vector
    zones = Internal.getZones(t)
    for z in zones:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured':
            i = dim[1]//2; j = dim[2]//2; k = dim[3]//2
            i = max(i,1); j = max(j,1); k = max(k,1)
            ip1 = i+1; jp1 = j+1; kp1 = k+1
            ip1 = min(ip1,dim[1]); jp1 = min(jp1,dim[2]); kp1 = min(kp1,dim[3])

            P0 = []; P1 = []; P2 = []; P3 = []

            x0 = C.getValue(z,'CoordinateX',(i,j,k)) ; P0.append(x0)
            y0 = C.getValue(z,'CoordinateY',(i,j,k)) ; P0.append(y0)
            z0 = C.getValue(z,'CoordinateZ',(i,j,k)) ; P0.append(z0)

            x1 = C.getValue(z,'CoordinateX',(ip1,j,k)) ; P1.append(x1)
            y1 = C.getValue(z,'CoordinateY',(ip1,j,k)) ; P1.append(y1)
            z1 = C.getValue(z,'CoordinateZ',(ip1,j,k)) ; P1.append(z1)

            x2 = C.getValue(z,'CoordinateX',(i,jp1,k)) ; P2.append(x2)
            y2 = C.getValue(z,'CoordinateY',(i,jp1,k)) ; P2.append(y2)
            z2 = C.getValue(z,'CoordinateZ',(i,jp1,k)) ; P2.append(z2)

            x3 = C.getValue(z,'CoordinateX',(i,j,kp1)) ; P3.append(x3)
            y3 = C.getValue(z,'CoordinateY',(i,j,kp1)) ; P3.append(y3)
            z3 = C.getValue(z,'CoordinateZ',(i,j,kp1)) ; P3.append(z3)

            l1 = Vector.sub(P1,P0); ln1 = Vector.norm2(l1)
            l2 = Vector.sub(P2,P0); ln2 = Vector.norm2(l2)
            l3 = Vector.sub(P3,P0); ln3 = Vector.norm2(l3)

            c = Vector.cross(l1, l2)
            if ln1 > 0 and ln2 > 0 and ln3 > 0: # 3D
                c = Vector.dot(c, l3)
                if c < 0: _reorder(z, (1,2,-3))
            elif ln1 > 0 and ln2 > 0: # 2D
                xmin = C.getMinValue(z, 'CoordinateX')
                xmax = C.getMaxValue(z, 'CoordinateX')
                ymin = C.getMinValue(z, 'CoordinateY')
                ymax = C.getMaxValue(z, 'CoordinateY')
                zmin = C.getMinValue(z, 'CoordinateZ')
                zmax = C.getMaxValue(z, 'CoordinateZ')
                if abs(zmax-zmin) < 1.e-6:
                    l3 = (0., 0., 1.)
                    c = Vector.dot(c, l3)
                    if c < 0: _reorder(z, (1,-2,3))
                elif abs(ymax-ymin < 1.e-6):
                    l3 = (0., 1., 0.)
                    c = Vector.dot(c, l3)
                    if c < 0: _reorder(z, (1,-2,3))
                elif abs(xmax-xmin < 1.e-6):
                    l3 = (1., 0., 0.)
                    c = Vector.dot(c, l3)
                    if c < 0: _reorder(z, (1,-2,3))
    return None

def addkplane(t, N=1):
    """Add N k-plane(s) to a mesh.
    Usage: addkplane(t, N)"""
    tp = Internal.copyRef(t)
    _addkplane(tp, N)
    return tp

def _addkplane(t, N=1):
    """Add N k-plane(s) to a mesh."""
    zones = Internal.getZones(t)
    for z in zones:
        nodes = C.getFields(Internal.__GridCoordinates__, z)[0]
        fn = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
        fc = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
        # Coordinates + fields located at nodes
        if fn != []:
            if nodes == []: nodes = fn
            else: nodes = Converter.addVars([nodes, fn])

        nodesN = Transform.addkplane(nodes, N=N)
        C.setFields([nodesN], z, 'nodes')
        # Fields located at centers
        if fc != []:
            centersN = Transform.addkplaneCenters(fc, nodesN, N=N)
            C.setFields([centersN], z, 'centers')
    # Modify BCs
    Internal._addOneLayer2BC(t, dir=3, N=N)
    return None

def projectAllDirs(t1, t2, vect=['nx','ny','nz'], oriented=0):
    """Project points defined in arrays to surfaces according to the direction provided by vect.
    Usage: projectAllDirs(arrays, surfaces, vect, oriented)
    """
    t = Internal.copyRef(t1)
    _projectAllDirs(t, t2, vect, oriented)
    return t

def _projectAllDirs(t1, t2, vect=['nx','ny','nz'], oriented=0):
    """Project points defined in arrays to surfaces according to the direction provided by vect."""
    zones = Internal.getZones(t1)
    a1 = C.getAllFields(zones, loc='nodes')
    a1 = Converter.extractVars(a1,['CoordinateX','CoordinateY','CoordinateZ']+vect)
    a2 = C.getFields(Internal.__GridCoordinates__, t2)
    res = Transform.projectAllDirs(a1, a2, vect, oriented)
    for noz in range(len(zones)):
        C.setFields([res[noz]], zones[noz], 'nodes')
    return None

def projectDir(t1, t2, dir, smooth=0, oriented=0):
    """Project a surface array onto surface arrays following dir.
    Usage: projectDir(t1, t2, dir)"""
    t = Internal.copyRef(t1)
    _projectDir(t, t2, dir, smooth, oriented)
    return t

def _projectDir(t1, t2, dir, smooth=0, oriented=0): # t1 is modified
    """Project a surface array onto surface arrays following dir."""
    zones = Internal.getZones(t1)
    a1 = C.getFields(Internal.__GridCoordinates__, zones)
    a2 = C.getFields(Internal.__GridCoordinates__, t2)
    res = Transform.projectDir(a1, a2, dir, smooth, oriented)
    for noz in range(len(zones)):
        C.setFields([res[noz]], zones[noz], 'nodes')
    return None

def projectOrtho(t1, t2):
    """Project a surface t1 onto surface t2 orthogonally.
    Usage: projectOrtho(t1, t2)"""
    t = Internal.copyRef(t1)
    _projectOrtho(t, t2)
    return t

def _projectOrtho(t1, t2): # t1 is modified
    """Project a surface t1 onto surface t2 orthogonally."""
    zones = Internal.getZones(t1)
    a1 = C.getFields(Internal.__GridCoordinates__, zones)
    a2 = C.getFields(Internal.__GridCoordinates__, t2)
    res = Transform.projectOrtho(a1, a2)
    for noz in range(len(zones)):
        C.setFields([res[noz]], zones[noz], 'nodes')
    return None

def projectOrthoSmooth(t1, t2, niter=1):
    """Project a surface array onto surface arrays following smoothed normals.
    Usage: projectOrthoSmooth(t1, t2)"""
    t = Internal.copyRef(t1)
    _projectOrthoSmooth(t, t2, niter)
    return t

def _projectOrthoSmooth(t1, t2, niter=1): # t1 is modified
    """Project a surface array onto surface arrays following smoothed normals."""
    zones = Internal.getZones(t1)
    a1 = C.getFields(Internal.__GridCoordinates__, zones)
    a2 = C.getFields(Internal.__GridCoordinates__, t2)
    res = Transform.projectOrthoSmooth(a1, a2, niter)
    for noz in range(len(zones)):
        C.setFields([res[noz]], zones[noz], 'nodes')
    return None

def projectRay(t1, t2, Pt):
    """Project a surface array onto surface arrays following rays.
    Usage: projectRay(t1, t2, Pt)"""
    t = Internal.copyRef(t1)
    _projectRay(t, t2, Pt)
    return t

def _projectRay(t1, t2, Pt): # t1 is modified
    """Project a surface array onto surface arrays following rays."""
    zones = Internal.getZones(t1)
    a1 = C.getFields(Internal.__GridCoordinates__, zones)
    a2 = C.getFields(Internal.__GridCoordinates__, t2)
    res = Transform.projectRay(a1, a2, Pt)
    for noz in range(len(zones)):
        C.setFields([res[noz]], zones[noz], 'nodes')
    return None

def alignVectorFieldWithRadialCylindricProjection(t, axisPassingPoint=(0,0,0), axisDirection=(1,0,0), vectorNames=[]):
    """Perform a cylindric radial projection of a vector field."""
    return C.TZA3(t, 'nodes', 'nodes', False, Transform.alignVectorFieldWithRadialCylindricProjection, axisPassingPoint, axisDirection, vectorNames)

def _alignVectorFieldWithRadialCylindricProjection(t, axisPassingPoint=(0,0,0), axisDirection=(1,0,0), vectorNames=[]):
    """Perform a cylindric radial projection of a vector field."""
    C.__TZA3(t, 'nodes', Transform._alignVectorFieldWithRadialCylindricProjection, axisPassingPoint, axisDirection, vectorNames)

# Split au milieu
def _splitSize__(z, N, multigrid, dirs, t, stack):
    dim = Internal.getZoneDim(z)
    if dim[0] == 'Unstructured':
        print('Warning: splitSize: unstructured zone not treated.')
        return
    ni = dim[1]; nj = dim[2]; nk = dim[3]
    if ni*nj*nk > N:
        dirl = Transform.getSplitDir__(ni, nj, nk, dirs)
        if dirl == 1:
            ns = Transform.findMGSplit__(ni, level=multigrid)
            if ns > 0:
                z1, z2 = split(z, 1, ns, t)
                stack.append(z1); stack.append(z2)
        elif dirl == 2:
            ns = Transform.findMGSplit__(nj, level=multigrid)
            if ns > 0:
                z1, z2 = split(z, 2, ns, t)
                stack.append(z1); stack.append(z2)
        elif dirl == 3:
            ns = Transform.findMGSplit__(nk, level=multigrid)
            if ns > 0:
                z1, z2 = split(z, 3, ns, t)
                stack.append(z1); stack.append(z2)
        else:
            ns = Transform.findMGSplit__(ni, level=multigrid)
            if ns > 0:
                z1, z2 = split(z, 1, ns, t)
                stack.append(z1); stack.append(z2)

# Split decentre
def _splitSizeUp__(z, N, multigrid, dirs, t, stack):
    dim = Internal.getZoneDim(z)
    if dim[0] == 'Unstructured':
        print('Warning: splitSize: unstructured zone not treated.')
        return

    ni = dim[1]; nj = dim[2]; nk = dim[3]
    nij = ni*nj; nik = ni*nk; njk = nj*nk
    if ni*nj*nk > N:
        dirl = Transform.getSplitDir__(ni, nj, nk, dirs)
        if dirl == 1:
            ns = Transform.findMGSplitUp__(ni, int(N/njk), level=multigrid)
            if ns > 0:
                z1, z2 = split(z, 1, ns, t)
                stack.append(z1); stack.append(z2)
        elif dirl == 2:
            ns = Transform.findMGSplitUp__(nj, int(N/nik), level=multigrid)
            if ns > 0:
                z1, z2 = split(z, 2, ns, t)
                stack.append(z1); stack.append(z2)
        elif dirl == 3:
            ns = Transform.findMGSplitUp__(nk, int(N/nij), level=multigrid)
            if ns > 0:
                z1, z2 = split(z, 3, ns, t)
                stack.append(z1); stack.append(z2)
        else:
            ns = Transform.findMGSplitUp__(ni, int(N/njk), level=multigrid)
            if ns > 0:
                z1, z2 = split(z, 1, ns, t)
                stack.append(z1); stack.append(z2)

# Return the number of cells in zone
def getNCells(z):
    """Return the number of cells in zone."""
    dim = Internal.getZoneDim(z)
    if dim[0] == 'Unstructured': return dim[2]
    else:
        ni = dim[1]; nj = dim[2]; nk = dim[3]
        ni1 = max(1, ni-1); nj1 = max(1, nj-1); nk1 = max(1, nk-1)
        return ni1*nj1*nk1

# Split size decentre avec ressources
def _splitSizeUpR__(t, N, R, multigrid, dirs, minPtsPerDir, topTree):
    bases = Internal.getBases(t)
    SP = []; Nl = 0
    for b in bases:
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for z in zones:
            dim = Internal.getZoneDim(z)
            if dim[0] == 'Unstructured':
                print('Warning: splitSize: unstructured zone skipped.')
            if dim[0] == 'Structured':
                ni = dim[1]; nj = dim[2]; nk = dim[3]
                ni1 = max(1, ni-1); nj1 = max(1, nj-1); nk1 = max(1, nk-1)
                SP.append((ni1*nj1*nk1,z,b)); Nl += ni1*nj1*nk1
    if N == 0: N = Nl*1. / R
    #print('average cells ', N)
    from operator import itemgetter

    # Init le vecteur des ressources
    Rs = [0 for i in range(R)]
    procs = [ [0,[]] for i in range(R)]
    mins = minPtsPerDir-1 # nbre de cellules mini des blocs par direction

    while len(SP) > 0:
        SP = sorted(SP, key=itemgetter(0), reverse=True)
        Rs = sorted(Rs)
        procs = sorted(procs, key=itemgetter(0))
        #print('ress', Rs[0], C.getNCells(SP[0][1]))
        a = SP[0][1] # le plus gros
        base = SP[0][2]
        dim = Internal.getZoneDim(a)
        ni = dim[1]; nj = dim[2]; nk = dim[3]
        ni1 = max(1, ni-1); nj1 = max(1, nj-1); nk1 = max(1, nk-1)
        nik = ni1*nk1; njk = nj1*nk1; nij = ni1*nj1
        Nr = min(N, N-Rs[0])
        ncells = ni1*nj1*nk1
        if ncells > Nr:
            # Calcul le meilleur split
            nc = int(round(Nr*1./njk,0))+1
            ns = Transform.findMGSplitUp__(ni, nc, level=multigrid)
            if ns-1 < mins: ns = 5
            delta1 = ns-1
            delta2 = ni-ns
            if delta2 < mins: delta2 -= 1.e6
            delta3 = abs((ns-1)*njk - Nr)/njk
            deltai = delta3-delta1-delta2
            nc = int(round(Nr*1./nik,0))+1
            ns = Transform.findMGSplitUp__(nj, nc, level=multigrid)
            if ns-1 < mins: ns = 5
            delta1 = ns-1
            delta2 = nj-ns
            if delta2 < mins: delta2 -= 1.e6
            delta3 = abs(ni1*(ns-1)*nk1 - Nr)/nik
            deltaj = delta3-delta1-delta2
            nc = int(round(Nr*1./nij,0))+1
            ns = Transform.findMGSplitUp__(nk, nc, level=multigrid)
            if ns-1 < mins: ns = 5
            delta1 = ns-1
            delta2 = nk-ns
            if delta2 < mins: delta2 -= 1.e6
            delta3 = abs(ni1*nj1*(ns-1) - Nr)/nij
            deltak = delta3-delta1-delta2
            dirl = 1
            if deltai <= deltaj  and deltai <= deltak:
                if 1 in dirs: dirl = 1
                elif deltaj <= deltak and 2 in dirs: dirl = 2
                elif 3 in dirs: dirl = 3
            elif deltaj <= deltai and deltaj <= deltak:
                if 2 in dirs: dirl = 2
                elif deltai <= deltak and 1 in dirs: dirl = 1
                elif 3 in dirs: dirl = 3
            elif deltak <= deltai and deltak <= deltaj:
                if 3 in dirs: dirl = 3
                elif deltai <= deltaj and 1 in dirs: dirl = 1
                elif 2 in dirs: dirl = 2

            trynext = 1
            if dirl == 1:
                nc = int(round(Nr*1./njk,0))+1
                ns = Transform.findMGSplitUp__(ni, nc, level=multigrid)
                if ns-1 >= mins and ni-ns >= mins:
                    a1, a2 = split(a, 1, ns, topTree)
                    SP[0] = (getNCells(a2), a2, base)
                    Rs[0] += getNCells(a1)
                    procs[0][1].append(a1[0])
                    procs[0][0] = Rs[0]
                    trynext = 0
            elif dirl == 2:
                nc = int(round(Nr*1./nik,0))+1
                ns = Transform.findMGSplitUp__(nj, nc, level=multigrid)
                if ns-1 >= mins and nj-ns >= mins:
                    a1, a2 = split(a, 2, ns, topTree)
                    SP[0] = (getNCells(a2), a2, base)
                    Rs[0] += getNCells(a1)
                    procs[0][1].append(a1[0])
                    procs[0][0] = Rs[0]
                    trynext = 0
            elif dirl == 3:
                nc = int(round(Nr*1./nij,0))+1
                ns = Transform.findMGSplitUp__(nk, nc, level=multigrid)
                if ns-1 >= mins and nk-ns >= mins:
                    a1, a2 = split(a, 3, ns, topTree)
                    SP[0] = (getNCells(a2), a2, base)
                    Rs[0] += getNCells(a1)
                    procs[0][1].append(a1[0])
                    procs[0][0] = Rs[0]
                    trynext = 0
            if trynext == 1:
                Rs[0] += getNCells(a); del SP[0]
                procs[0][1].append(a[0])
                procs[0][0] = Rs[0]
        else:
            Rs[0] += getNCells(a); del SP[0]
            procs[0][1].append(a[0])
            procs[0][0] = Rs[0]

    # Affectation des procs
    try:
        import Distributor2.PyTree as D2
        for np, p in enumerate(procs):
            for zname in p[1]:
                z = Internal.getNodeFromName(t, zname)
                D2._addProcNode(z, np)
    except: pass
    return None

# Split size decentre avec ressources
def splitSizeUpR_OMP__(t, N, R, multigrid, dirs, minPtsPerDir):
    ific = 0
    bases = Internal.getBases(t)
    SP = []; Nl = 0
    zsplit = []
    for b in bases:
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for z in zones:
            dim = Internal.getZoneDim(z)
            if dim[0] == 'Unstructured':
                print('Warning: splitSize: unstructured zone not treated.')
            if dim[0] == 'Structured':
                ni = dim[1]-2*ific; nj = dim[2]-2*ific; nk = dim[3]-2*ific
                ni1 = max(1, ni-1); nj1 = max(1, nj-1); nk1 = max(1, nk-1)
                SP.append((ni1*nj1*nk1,[z[0],[ni,nj,nk]])); Nl += ni1*nj1*nk1
                zsplit.append(([1,ni1,1,nj1,1,nk1],(ni,nj,nk),z[0],"master",z[0]))

    if N == 0: N = Nl*1. / R
    #print 'average cells ', N
    from operator import itemgetter

    # Init le vecteur des ressources
    Rs = [0]*R
    Thread_z = [([0],1,[])]
    for ith in range(2,R+1): Thread_z.append(([0],ith,[])) # Thread_z = [Nbpoints,ithread,corresponding work]

    mins = minPtsPerDir-1 # nbre de cellules mini des blocs par direction

    out = []
    nbl = 0
    while len(SP) > 0:
        SP = sorted(SP, key=itemgetter(0), reverse=True)
        Rs = sorted(Rs)
        Thread_z = sorted(Thread_z, key=itemgetter(0))
        #print 'ress', Rs[0], C.getNCells(SP[0][1])
        a = SP[0][1] # le plus gros
        dim = a[1]

        ni = dim[0]; nj = dim[1]; nk = dim[2]
        ni1 = max(1, ni-1); nj1 = max(1, nj-1); nk1 = max(1, nk-1)
        nik = ni1*nk1; njk = nj1*nk1; nij = ni1*nj1
        Nr = min(N, N-Rs[0])
        ncells = ni1*nj1*nk1
        if ncells > Nr:
            # Calcul le meilleur split
            nc = int(round(Nr*1./njk,0))+1
            ns = Transform.findMGSplitUp__(ni, nc, level=multigrid)
            if ns-1 < mins: ns = 5
            delta1 = ns-1
            delta2 = ni-ns
            if delta2 < mins: delta2 -= 1.e6
            delta3 = abs((ns-1)*njk - Nr)/njk
            deltai = delta3-delta1-delta2
            nc = int(round(Nr*1./nik,0))+1
            ns = Transform.findMGSplitUp__(nj, nc, level=multigrid)
            if ns-1 < mins: ns = 5
            delta1 = ns-1
            delta2 = nj-ns
            if delta2 < mins: delta2 -= 1.e6
            delta3 = abs(ni1*(ns-1)*nk1 - Nr)/nik
            deltaj = delta3-delta1-delta2
            nc = int(round(Nr*1./nij,0))+1
            ns = Transform.findMGSplitUp__(nk, nc, level=multigrid)
            if ns-1 < mins: ns = 5
            delta1 = ns-1
            delta2 = nk-ns
            if delta2 < mins: delta2 -= 1.e6
            delta3 = abs(ni1*nj1*(ns-1) - Nr)/nij
            deltak = delta3-delta1-delta2
            dirl = 1
#            if (deltai <= deltaj  and deltai <= deltak):
#                if (1 in dirs): dirl = 1
#                elif (deltaj <= deltak and 2 in dirs): dirl = 2
#                elif (3 in dirs): dirl = 3
#            elif (deltaj <= deltai and deltaj <= deltak):
#                if (2 in dirs): dirl = 2
#                elif (deltai <= deltak and 1 in dirs): dirl = 1
#                elif (3 in dirs): dirl = 3
#            elif (deltak <= deltai and deltak <= deltaj):
#                if (3 in dirs): dirl = 3
#                elif (deltai <= deltaj and 1 in dirs): dirl = 1
#                elif (2 in dirs): dirl = 2
#
            if deltak <= deltai and deltak <= deltaj:            # Favor k split over j and i
                if 3 in dirs: dirl = 3
                elif deltai <= deltaj and 1 in dirs: dirl = 1
                elif 2 in dirs: dirl = 2
            elif deltaj <= deltai and deltaj <= deltak:          # then try j split
                if 2 in dirs: dirl = 2
                elif deltai <= deltak and 1 in dirs: dirl = 1
                elif 3 in dirs: dirl = 3
            elif deltai <= deltaj  and deltai <= deltak:         # and i split if needed
                if 1 in dirs: dirl = 1
                elif deltaj <= deltak and 2 in dirs: dirl = 2
                elif 3 in dirs: dirl = 3

            trynext = 1

            if dirl == 3:
                nc = int(round(Nr*1./nij,0))+1
                ns = Transform.findMGSplitUp__(nk, nc, level=multigrid)
                if ns-1 >= mins and nk-ns >= mins:
                    #a1 = subzone(a, (1,1,1), (ni,nj,ns))
                    #a2 = subzone(a, (1,1,ns), (ni,nj,nk))
                    a1 = ["leafl"+str(('%05d' % nbl)),[ni,nj,ns]]
                    a2 = ["leafr"+str(('%05d' % nbl)),[ni,nj,nk-ns+1]]
                    SP[0] = ((ni-1)*(nj-1)*(nk-ns),a2)
                    nbl=nbl+1
                    Rs[0] += (ni-1)*(nj-1)*(ns-1)
                    Thread_z[0][0][0] += (ni-1)*(nj-1)*(ns-1)
                    Thread_z[0][2].append(a1[0])
                    trynext = 0
            elif dirl == 2:
                nc = int(round(Nr*1./nik,0))+1
                ns = Transform.findMGSplitUp__(nj, nc, level=multigrid)
                if ns-1 >= mins and nj-ns >= mins:
                    #a1 = subzone(a, (1,1,1), (ni,ns,nk))
                    #a2 = subzone(a, (1,ns,1), (ni,nj,nk))
                    #SP[0] = (getNCells(a2), a2, base)
                    a1 = ["leafl"+str(('%05d' % nbl)),[ni,ns,nk]]
                    a2 = ["leafr"+str(('%05d' % nbl)),[ni,nj-ns+1,nk]]
                    SP[0] = ((ni-1)*(nj-ns)*(nk-1),a2)
                    nbl += 1
                    Rs[0] += (ni-1)*(ns-1)*(nk-1)
                    Thread_z[0][0][0] += (ni-1)*(ns-1)*(nk-1)
                    Thread_z[0][2].append(a1[0])
                    trynext = 0
            elif dirl == 1:
                nc = int(round(Nr*1./njk,0))+1
                ns = Transform.findMGSplitUp__(ni, nc, level=multigrid)
                if ns-1 >= mins and ni-ns >= mins:
                    #a1 = subzone(a, (1,1,1), (ns,nj,nk))
                    #a2 = subzone(a, (ns,1,1), (ni,nj,nk))
                    #SP[0] = (getNCells(a2), a2, base)
                    a1 = ["leafl"+str(('%05d' % nbl)),[ns,nj,nk]]
                    a2 = ["leafr"+str(('%05d' % nbl)),[ni-ns+1,nj,nk]]
                    SP[0] = ((ni-ns)*(nj-1)*(nk-1),a2)
                    nbl += 1
                    Rs[0] += (ns-1)*(nj-1)*(nk-1)
                    Thread_z[0][0][0] += (ns-1)*(nj-1)*(nk-1)
                    Thread_z[0][2].append(a1[0])
                    trynext = 0

            if (dirl>=1)&(trynext==0):
                test=[s for s in zsplit if a[0] in s]

                dima1 = a1[1]
                dima2 = a2[1]

                if test != []:

                    indexleft  = numpy.zeros(6, dtype=Internal.E_NpyInt)
                    indexright = numpy.zeros(6, dtype=Internal.E_NpyInt)
                    indexfather= test[0][0]
                    dimfather  = test[0][1]

                    indexleft[0]=indexfather[0];indexleft[1]=indexfather[1]
                    indexleft[2]=indexfather[2];indexleft[3]=indexfather[3]
                    indexleft[4]=indexfather[4];indexleft[5]=indexfather[5]

                    indexright[0]=indexfather[0];indexright[1]=indexfather[1]
                    indexright[2]=indexfather[2];indexright[3]=indexfather[3]
                    indexright[4]=indexfather[4];indexright[5]=indexfather[5]

                    if dima1[0] != dimfather[0]:
                        indexleft[1]  = indexleft[0] + dima1[0] -2-2*ific
                        indexright[0] = indexleft[1] + 1
                    if dima1[1] != dimfather[1]:
                        indexleft[3]  = indexleft[2] + dima1[1] -2-2*ific
                        indexright[2] = indexleft[3] + 1
                    if dima1[2] != dimfather[2]:
                        indexleft[5]  = indexleft[4] + dima1[2] -2-2*ific
                        indexright[4] = indexleft[5] + 1

                    zsplit.append(([ indexleft[0] ,indexleft[1] ,indexleft[2] ,indexleft[3], indexleft[4], indexleft[5]],(dima1[0]-2*ific,dima1[1]-2*ific,dima1[2]-2*ific),a1[0],a[0],test[0][4]))
                    zsplit.append(([indexright[0],indexright[1],indexright[2],indexright[3],indexright[4],indexright[5]],(dima2[0]-2*ific,dima2[1]-2*ific,dima2[2]-2*ific),a2[0],a[0],test[0][4]))

            if trynext == 1:
                Thread_z[0][0][0] += (a[1][0]-1)*(a[1][1]-1)*(a[1][2]-1)
                Thread_z[0][2].append(a[0])
                Rs[0] += (a[1][0]-1)*(a[1][1]-1)*(a[1][2]-1); del SP[0]

        else:
            Rs[0] += (a[1][0]-1)*(a[1][1]-1)*(a[1][2]-1); del SP[0]
            Thread_z[0][0][0] += (a[1][0]-1)*(a[1][1]-1)*(a[1][2]-1)
            Thread_z[0][2].append(a[0])

    zomp_threads={}

    for z in zones:
        th={}
        for r in range(1,R+1):
            th[r]=[]
        zomp_threads[z[0]]=th

    for ith, listth in enumerate(Thread_z): # get the subzone list for each threads
        for zleaf in enumerate(Thread_z[ith][2]):
            zind = [s for s in zsplit if zleaf[1] in s[2][:] ] # get the indexes
            zomp_threads[zind[0][4]][ith+1].append((zleaf[1],zind[0][0]))

    t = C.addVars(t,'centers:thread_number')
    t = C.addVars(t,'centers:thread_subzone')
    bases = Internal.getBases(t)
    for b in bases:
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for z in zones:
            solverParam=Internal.createChild(z,'.Solver#Param','UserDefinedData_t',value=None,children=[],pos=-1)
            ompthread  =Internal.createChild(solverParam,  'omp_threads'  ,'UserDefinedData_t',value=None,children=[],pos=-1)
            sol  = Internal.getNodeFromName1(z, 'FlowSolution#Centers')
            solth=Internal.getNodeFromName1(sol, 'thread_number')
            solsz=Internal.getNodeFromName1(sol, 'thread_subzone')
            for i in range(1,R+1):
                thnode=[]
                thnode=Internal.createChild(ompthread,str(i),'UserDefinedData_t',value=None,children=[],pos=-1)
                isb=0
                for ilisth in range(0,len(zomp_threads[z[0]][i])):
                    isb=isb+1
                    if zomp_threads[z[0]][i][ilisth][1] != []:
                        Internal.createChild(thnode,'subzone'+str(isb),'DataArray_t',value=zomp_threads[z[0]][i][ilisth][1],children=[],pos=-1)
                        ind=zomp_threads[z[0]][i][ilisth][1]
                        l=0
                        for iths in range(ind[0]-1,ind[1]):
                            for jths in range(ind[2]-1,ind[3]):
                                for kths in range(ind[4]-1,ind[5]):
                                    solth[1][iths][jths][kths]=i
                                    solsz[1][iths][jths][kths]=isb

    print('ress:', Rs)
    Tot = 0
    for i in Rs: Tot += i
    print('Tot', Tot)
    print("Imbalance",min(Rs)/(max(Rs)*1.0))
    return t

def _splitNParts(t, N, multigrid=0, dirs=[1,2,3], recoverBC=True, topTree=None):
    if topTree is None: topTree = t
    zones = Internal.getZones(t)
    # Fait des paquets de zones structurees et non structurees
    zonesS = []; zonesN = []
    NpS = []; NpN = [] # nbre de points
    NeS = []; NeN = [] # nbre de cellules
    outO = []
    for z in zones:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured':
            zonesS.append(z)
            NpS.append(dim[1]*dim[2]*dim[3])
            NeS.append(max(dim[1]-1,1)*max(dim[2]-1,1)*max(dim[3]-1,1))
        else:
            zonesN.append(z)
            NpN.append(dim[1])
            NeN.append(dim[2])

    SumS = 0.; SumN = 0.
    for i in NeS: SumS += i
    for i in NeN: SumN += i
    if SumS+SumN < 0.01: return outO
    alpha = N*1./(SumS+SumN)
    NbN = len(NeN) # nbre de grilles non structurees
    NPart = [0]*(NbN+1); Nt = 0
    for i in range(NbN): NPart[i] = max(int(alpha*NeN[i]),1); Nt += NPart[i]
    if SumS != 0: NPart[NbN] = max(N-Nt,1)
    else: NPart[NbN-1] = max(N-Nt+NPart[NbN-1], 1)

    # Blocs non structures
    for i in range(len(zonesN)):
        z = zonesN[i]
        dim = Internal.getZoneDim(z)
        a = C.getFields(Internal.__GridCoordinates__, z)[0]
        zlist = []
        if NPart[i] > 1:
            if recoverBC: bcs = C.getBCs(z)
            if dim[3] == 'NGON':
                elts = transform.splitNGon(a, NPart[i])
            else: elts = transform.splitElement(a, NPart[i])
            for e in elts:
                zL = subzone(z, e, type='elements')
                if recoverBC: C._recoverBCs(zL, bcs)
                zlist.append(zL)
            _replaceZoneWithSplit(t, z[0], zlist[0], zlist[1:])

    # Blocs structures
    l = len(zonesS)
    if l == 0: return None
    NPa = NPart[NbN]
    Ns = Transform.findNsi__(l, NPa, NpS)
    #for z in zonesS: print(z[0])
    for i in range(l):
        a = zonesS[i]
        dimL = Internal.getZoneDim(a)
        ni = dimL[1]; nj = dimL[2]; nk = dimL[3]
        splits = Transform.findSplits__(ni, nj, nk, Ns[i], dirs, multigrid)
        # Find split direction in splits
        size = len(splits)
        nks = []
        k = 0
        while k < size-1:
            (i1,i2,j1,j2,k1,k2) = splits[k]
            (ip1,ip2,jp1,jp2,kp1,kp2) = splits[k+1]
            if i1 == ip1 and i2 == ip2 and j1 == jp1 and j2 == jp2: dir = 3
            elif i1 == ip1 and i2 == ip2 and k1 == kp1 and k2 == kp2: dir = 2
            elif j1 == jp1 and k1 == kp1 and j2 == jp2 and k2 == kp2: dir = 1
            else: dir = 4
            if dir == 3: nks.append(k2)
            else: break
            k += 1
        nk = len(nks); nk1 = nk+1
        njs = []
        k = 0
        while k < size//nk1-1:
            (i1,i2,j1,j2,k1,k2) = splits[k*nk1]
            (ip1,ip2,jp1,jp2,kp1,kp2) = splits[k*nk1+nk1]
            if i1 == ip1 and i2 == ip2 and j1 == jp1 and j2 == jp2: dir = 3
            elif i1 == ip1 and i2 == ip2 and k1 == kp1 and k2 == kp2: dir = 2
            elif j1 == jp1 and k1 == kp1 and j2 == jp2 and k2 == kp2: dir = 1
            else: dir = 4
            if dir == 2: njs.append(j2)
            else: break
            k += 1
        nj = len(njs); nj1 = nj+1
        nis = []
        k = 0
        while k < size//(nj1*nk1)-1:
            (i1,i2,j1,j2,k1,k2) = splits[k*nj1*nk1]
            (ip1,ip2,jp1,jp2,kp1,kp2) = splits[k*nj1*nk1+nj1*nk1]
            if i1 == ip1 and i2 == ip2 and j1 == jp1 and j2 == jp2: dir = 3
            elif i1 == ip1 and i2 == ip2 and k1 == kp1 and k2 == kp2: dir = 2
            elif j1 == jp1 and k1 == kp1 and j2 == jp2 and k2 == kp2: dir = 1
            else: dir = 4
            if dir == 1: nis.append(i2)
            else: break
            k += 1
        ni = len(nis)
        # decalage des nis,njs,nks
        prev = 0
        for e, i in enumerate(nis):
            nis[e] -= prev
            prev += nis[e]-1
        prev = 0
        for e, i in enumerate(njs):
            njs[e] -= prev
            prev += njs[e]-1
        prev = 0
        for e, i in enumerate(nks):
            nks[e] -= prev
            prev += nks[e]-1
        store = []
        ap = a

        for e, i in enumerate(nis):
            a1,a2 = split(ap, 1, i, topTree)
            store.append(a1)
            if e == len(nis)-1: store.append(a2)
            ap = a2

        if len(nis) == 0: store = [a]
        store2 = []
        for ap in store:
            for e, j in enumerate(njs):
                a1,a2 = split(ap, 2, j, topTree)
                store2.append(a1)
                if e == len(njs)-1: store2.append(a2)
                ap = a2
        if len(njs) == 0: store2 = store
        for ap in store2:
            for e, k in enumerate(nks):
                a1,a2 = split(ap, 3, k, topTree)
                ap = a2
    return None

def splitNParts(t, N, multigrid=0, dirs=[1,2,3], recoverBC=True, topTree=None):
    """Split zones in t in N parts. 
    Usage: splitNParts(t, N, multigrid, dirs, recoverBC)"""
    tp = Internal.copyRef(t)
    tpp, typen = Internal.node2PyTree(tp)
    _splitNParts(tpp, N, multigrid, dirs, recoverBC, topTree)
    tp = Internal.pyTree2Node(tpp, typen)
    return tp

def splitSize(t, N=0, multigrid=0, dirs=[1,2,3], type=0, R=None,
              minPtsPerDir=5, topTree=None):
    """Split zones following size."""
    tp = Internal.copyRef(t)
    tpp, typen = Internal.node2PyTree(tp)
    _splitSize(tpp, N, multigrid, dirs, type, R, minPtsPerDir, topTree)
    tp = Internal.pyTree2Node(tpp, typen)
    return tp

def _splitSize(t, N=0, multigrid=0, dirs=[1,2,3], type=0, R=None,
               minPtsPerDir=5, topTree=None):
    minPtsPerDir = max(minPtsPerDir, 2**(multigrid+1)+1)
    if R is not None: type = 2
    if topTree is None: topTree = t
    if type == 0:
        stack = []
        for z in Internal.getZones(t): stack.append(z)
        while len(stack)>0:
            z = stack.pop()
            _splitSize__(z, N, multigrid, dirs, topTree, stack)
    elif type == 1:
        stack = []
        for z in Internal.getZones(t): stack.append(z)
        while len(stack)>0:
            z = stack.pop()
            _splitSizeUp__(z, N, multigrid, dirs, topTree, stack)
    elif type == 2:
        _splitSizeUpR__(t, N, R, multigrid, dirs, minPtsPerDir, topTree)
    return None

def splitCurvatureAngle(t, sensibility):
    """Split a curve following curvature angle.
    Usage: splitCurvatureAngle(t, sensibility)"""
    a = C.getAllFields(t, 'nodes')[0]
    arrays = Transform.splitCurvatureAngle(a, sensibility)
    zones = []
    for i in arrays:
        zone = C.convertArrays2ZoneNode('split', [i])
        zones.append(zone)
    return zones

def splitConnexity(t):
    """Split zone into connex zones.
    Usage: splitConnexity(t)"""
    a = C.getAllFields(t, 'nodes')[0]
    A = Transform.splitConnexity(a)
    zones = []
    for i in A:
        zone = C.convertArrays2ZoneNode('split', [i])
        zones.append(zone)
    return zones

def breakElements(t):
    """Break a NGON array in a set of arrays of BAR, TRI, ... elements.
    Usage: breakElements(t)"""
    a = C.getAllFields(t, 'nodes')
    A = Transform.breakElements(a)
    zones = []
    for i in A:
        if len(i) == 5: name = 'Struct'
        else: name = i[3]
        zone = C.convertArrays2ZoneNode(name, [i])
        zones.append(zone)
    return zones

def dual(t, extraPoints=1):
    """Return the dual mesh of a conformal mesh.
    Usage: dual(t, extraPoints)"""
    t = C.deleteFlowSolutions__(t, 'centers')
    return C.TZA1(t, 'nodes', 'nodes', True, Transform.dual, extraPoints)

def _dual(t, extraPoints=1):
    C._deleteFlowSolutions__(t, 'centers')
    return C._TZA1(t, 'nodes', 'nodes', True, Transform.dual, extraPoints)

def splitSharpEdges(t, alphaRef=30.):
    """Split zone into smooth zones.
    Usage: splitSharpEdges(t, alphaRef)"""
    a = C.getAllFields(t, 'nodes')[0]
    A = Transform.splitSharpEdges(a, alphaRef)
    zones = []
    for i in A:
        zone = C.convertArrays2ZoneNode('split',[i])
        zones.append(zone)
    return zones

def splitCurvatureRadius(t, Rs=100.):
    """Return the indices of the array where the curvature radius is low.
    Usage: splitCurvatureRadius(t, Rs)"""
    a = C.getAllFields(t, 'nodes')[0]
    arrays = Transform.splitCurvatureRadius(a, Rs)
    zones = []
    for i in arrays:
        zone = C.convertArrays2ZoneNode('split',[i])
        zones.append(zone)
    return zones

# -- SplitMultiplePts --
# detruit les BCMatch qui sont en raccord avec zoneName
def _deleteInconsistentBCMatch__(t, zoneName):
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
        for i in nodes:
            donorName = Internal.getValue(i)
            if donorName == zoneName:
                (parent, d) = Internal.getParentOfNode(z, i)
                del parent[2][d]
    return None

def splitMultiplePts3D__(t):
    restart = 0
    bases = Internal.getBases(t)
    for b in bases:
        (parent0,nob) = Internal.getParentOfNode(t, b)
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for z in zones:
            taga = C.getField('definedBC', z)[0]
            parent,d2 = Internal.getParentOfNode(b, z)
            dims = Internal.getZoneDim(z); ni = dims[1]; nj = dims[2]; nk = dims[3]; ninj = ni*nj
            isplit = -1; jsplit = -1; ksplit = -1

            # detecte si un pt interieur est de tag > 1
            ni1 = max(2,ni-1); nj1 = max(2,nj-1); nk1 = max(2,nk-1)
            # fenetre i = 1 et ni
            for i in [1,ni]:
                for k in range(1,nk1):
                    for j in range(1,nj1):
                        ind = i-1+j*ni+k*ninj
                        if taga[1][0,ind] > 1.:
                            jsplit = j+1; ksplit = k+1
                            # split en k ?
                            if nk > 2:
                                indm = ind-ni; indp = ind+ni
                                if taga[1][0,indm] > 1. and taga[1][0,indp] > 1.:
                                    z1 = subzone(z,(1,1,1),(ni,nj,ksplit))
                                    z2 = subzone(z,(1,1,ksplit),(ni,nj,nk))
                                    del parent[2][d2]; parent0[2][nob] = parent
                                    _deleteInconsistentBCMatch__(parent0[2][nob][2],z[0])
                                    parent0[2][nob][2]+=[z1,z2]
                                    restart = 1
                                    return t, restart
                            # split en j ?
                            indm = ind-ninj; indp = ind+ninj
                            if indp > ninj*nk-1: indp = ind
                            if taga[1][0,indm] > 1. and taga[1][0,indp] > 1.:
                                z1 = subzone(z,(1,1,1),(ni,jsplit,nk))
                                z2 = subzone(z,(1,jsplit,1),(ni,nj,nk))
                                del parent[2][d2]; parent0[2][nob] = parent
                                _deleteInconsistentBCMatch__(parent0[2][nob][2],z[0])
                                parent0[2][nob][2]+=[z1,z2]
                                restart = 1
                                return t, restart

            # fenetre j = 1 et nj
            for j in [1,nj]:
                for k in range(1,nk1):
                    for i in range(1,ni1):
                        ind = i+(j-1)*ni+k*ninj
                        if taga[1][0,ind] > 1.:
                            isplit = i+1; ksplit = k+1
                            # split en k ?
                            if nk > 2:
                                indm = ind-1; indp = ind+1
                                if taga[1][0,indm] > 1. and taga[1][0,indp] > 1.:
                                    z1 = subzone(z,(1,1,1),(ni,nj,ksplit))
                                    z2 = subzone(z,(1,1,ksplit),(ni,nj,nk))
                                    del parent[2][d2]; parent0[2][nob] = parent
                                    _deleteInconsistentBCMatch__(parent0[2][nob][2],z[0])
                                    parent0[2][nob][2]+=[z1,z2]
                                    restart = 1
                                    return t, restart
                            # split en i ?
                            indm = ind-ninj; indp = ind+ninj
                            if indp > ninj*nk-1: indp = ind
                            if taga[1][0,indm] > 1. and taga[1][0,indp] > 1.:
                                z1 = subzone(z,(1,1,1),(isplit,nj,nk))
                                z2 = subzone(z,(isplit,1,1),(ni,nj,nk))
                                del parent[2][d2]; parent0[2][nob] = parent
                                _deleteInconsistentBCMatch__(parent0[2][nob][2],z[0])
                                parent0[2][nob][2]+=[z1,z2]
                                restart = 1
                                return t, restart
            # fenetre k = 1
            if nk > 2:
                for k in [1,nk]:
                    for j in range(1,nj1):
                        for i in range(1,ni1):
                            ind = i+j*ni+(k-1)*ninj
                            if taga[1][0,ind] > 1.:
                                isplit = i+1; jsplit = j+1
                                # split en i ?
                                indm = ind-ni; indp = ind+ni
                                if taga[1][0,indm] > 1. and taga[1][0,indp] > 1.:
                                    z1 = subzone(z,(1,1,1),(isplit,nj,nk))
                                    z2 = subzone(z,(isplit,1,1),(ni,nj,nk))
                                    del parent[2][d2]; parent0[2][nob] = parent
                                    _deleteInconsistentBCMatch__(parent0[2][nob][2],z[0])
                                    parent0[2][nob][2]+=[z1,z2]
                                    restart = 1
                                    return t, restart
                                # split en j ?
                                indm = ind-1; indp = ind+1
                                if taga[1][0,indm] > 1. and taga[1][0,indp] > 1.:
                                    z1 = subzone(z,(1,1,1),(ni,jsplit,nk))
                                    z2 = subzone(z,(1,jsplit,1),(ni,nj,nk))
                                    del parent[2][d2]; parent0[2][nob] = parent
                                    _deleteInconsistentBCMatch__(parent0[2][nob][2],z[0])
                                    parent0[2][nob][2]+=[z1,z2]
                                    restart = 1
                                    return t, restart
    return t, restart

def splitMultiplePts2D__(t):
    restart = 0
    bases = Internal.getBases(t)
    for b in bases:
        (parent0,nob) = Internal.getParentOfNode(t, b)
        zones = Internal.getNodesFromType1(b,'Zone_t')
        for z in zones:
            taga = C.getField('definedBC', z)[0]
            parent,d2 = Internal.getParentOfNode(b, z)
            dims = Internal.getZoneDim(z); ni = dims[1]; nj = dims[2]; nk = dims[3]; ninj = ni*nj
            isplit = -1; jsplit = -1
            # detecte si un pt interieur est de tag > 1
            # fenetre i = 1
            for j in range(1,nj-1):
                if taga[1][0,j*ni] > 1.:
                    isplit = 1; jsplit = j+1
                    z1 = subzone(z,(1,1,1),(ni,jsplit,nk))
                    z2 = subzone(z,(1,jsplit,1),(ni,nj,nk))
                    del parent[2][d2]; parent0[2][nob] = parent
                    _deleteInconsistentBCMatch__(parent0[2][nob][2],z[0])
                    parent0[2][nob][2]+=[z1,z2]
                    restart = 1
                    return t, restart
            # fenetre i = ni
            for j in range(1,nj-1):
                if taga[1][0,ni-1+j*ni] > 1.:
                    isplit = ni; jsplit = j+1
                    z1 = subzone(z,(1,1,1),(ni,jsplit,nk))
                    z2 = subzone(z,(1,jsplit,1),(ni,nj,nk))
                    del parent[2][d2]; parent0[2][nob] = parent
                    _deleteInconsistentBCMatch__(parent0[2][nob][2],z[0])
                    parent0[2][nob][2]+=[z1,z2]
                    restart = 1
                    return t, restart
            # fenetre j = 1
            for i in range(1,ni-1):
                if taga[1][0,i] > 1.:
                    isplit = i+1; jsplit = 1
                    z1 = subzone(z,(1,1,1),(isplit,nj,nk))
                    z2 = subzone(z,(isplit,1,1),(ni,nj,nk))
                    del parent[2][d2]; parent0[2][nob] = parent
                    _deleteInconsistentBCMatch__(parent0[2][nob][2],z[0])
                    parent0[2][nob][2]+=[z1,z2]
                    restart = 1
                    return t, restart
            # fenetre j = nj
            for i in range(1,ni-1):
                if taga[1][0,i+(nj-1)*ni] > 1.:
                    isplit = i+1; jsplit = nj
                    z1 = subzone(z,(1,1,1),(isplit,nj,nk))
                    z2 = subzone(z,(isplit,1,1),(ni,nj,nk))
                    del parent[2][d2]; parent0[2][nob] = parent
                    _deleteInconsistentBCMatch__(parent0[2][nob][2],z[0])
                    parent0[2][nob][2]+=[z1,z2]
                    restart = 1
                    return t, restart
    return t, restart

def splitMultiplePts__(tp, dim):
    try: import Connector.PyTree as X
    except:
        raise ImportError("splitMultiplePts requires Connector.PyTree module.")
    tp = C.rmVars(tp, 'definedBC')
    tp = X.connectMatch(tp, dim=dim, tol=1.e-6)

    # tag des zones de tp: tag est incremente pour des qu une fenetre est en raccord avec une autre fenetre
    zones = Internal.getZones(tp)
    for z in zones:
        dims = Internal.getZoneDim(z); ni = dims[1]; nj = dims[2]; nk = dims[3]; ninj = ni*nj
        taga = Converter.array('definedBC', ni, nj, nk)
        # BCNearMatch
        bnds = Internal.getNodesFromType2(z, 'GridConnectivity_t')
        for bc in bnds:
            type = Internal.getNodeFromName1(bc, 'GridConnectivityType')
            if type is not None:
                val = Internal.getValue(type)
                if val == 'Abutting':
                    range0 = Internal.getNodeFromName1(bc, 'PointRange')
                    r = range0[1]; [i1,i2,j1,j2,k1,k2] = Internal.range2Window(r)
                    for k in range(k1-1,k2):
                        for j in range(j1-1,j2):
                            for i in range(i1-1,i2):
                                ind = i + j *ni + k *ninj
                                taga[1][0,ind] += 1
        # BC match
        bnds = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
        for bc in bnds:
            range0 = Internal.getNodeFromName1(bc, 'PointRange')
            r = range0[1]
            [i1,i2,j1,j2,k1,k2] = Internal.range2Window(r)
            for k in range(k1-1,k2):
                for j in range(j1-1,j2):
                    for i in range(i1-1,i2):
                        ind = i+j*ni+k*ninj
                        taga[1][0,ind] += 1

        z = C.setFields([taga],z,'nodes')

    # split des zones si tagDefinedBC est superieur a 1 en un pt (2D) ou sur une ligne i,j ou k (3D)
    split = 1; count = 0
    while split == 1:
        count += 1
        if dim == 2: tp, split = splitMultiplePts2D__(tp)
        else: tp, split = splitMultiplePts3D__(tp)
    tp = C.rmVars(tp,'definedBC')
    return tp, count

def splitMultiplePts(t, dim=3):
    """Split blocks in tree if connected to several blocks at a given border.
    Usage: splitMultiplePts(t, dim)"""
    type = 0 # O: toptree, 1: liste de bases, 2: base, 3: liste de zones, 4: zone
    toptree = Internal.isTopTree(t)
    if toptree:
        type = 0; tp = Internal.copyRef(t)
    else:
        # t base ou non
        bases = Internal.getBases(t)
        if bases != []:
            stdNode = Internal.isStdNode(t)
            if stdNode == 0: type = 1 # liste de bases
            else: type = 2 # une base
            tp = C.newPyTree(); tp[2][1:] = bases

        else: # liste zones ou zone ?
            stdNode = Internal.isStdNode(t)
            if stdNode == 0: type = 3 # liste de zones
            else: type = 4 # une zone
            zones = Internal.getNodesFromType(t, 'Zone_t')
            tp = C.newPyTree(['Base']); tp[2][1][2] = zones

    count = 2
    while count > 1:
        tp, count = splitMultiplePts__(tp, dim)
    if type == 0: return tp
    elif type == 1: return Internal.getBases(tp)
    elif type == 2: return Internal.getBases(tp)[0]
    elif type == 3: return Internal.getZones(tp)
    else: return Internal.getZones(tp)

def splitBAR(t, N, N2=-1):
    """Split a BAR at index N (start 0).
    Usage: splitBAR(t, N)"""
    a = C.getAllFields(t, 'nodes')[0]
    A = Transform.splitBAR(a, N, N2)
    zones = []
    for i in A:
        zone = C.convertArrays2ZoneNode('split',[i])
        zones.append(zone)
    return zones

def splitTBranches(t, tol=1.e-10):
    """Split a BAR at vertices where T-branches exist.
    Usage: splitTBranches(t, tol)"""
    a = C.getAllFields(t, 'nodes')
    A = Transform.splitTBranches(a, tol)
    zones = []
    for i in A:
        zone = C.convertArrays2ZoneNode('split',[i])
        zones.append(zone)
    return zones

def splitTRI(t, idxList):
    """Split a TRI into several TRIs delimited by the input poly line.
    Usage: splitTRI(t, idxList)"""
    a = C.getAllFields(t, 'nodes')[0]
    A = Transform.splitTRI(a, idxList)
    zones = []
    for i in A:
        zone = C.convertArrays2ZoneNode('split',[i])
        zones.append(zone)
    return zones

def splitManifold(t):
    """Split an unstructured mesh (only TRI or BAR currently) into several manifold pieces.
    Usage: splitManifold(array)"""
    a = C.getAllFields(t, 'nodes')
    A = Transform.splitManifold(a)
    zones = []
    for i in A:
        zone = C.convertArrays2ZoneNode('manifold',[i])
        zones.append(zone)
    return zones

def _splitNGon(t, N, N2=-1, shift=1000):
    """Split a NGon putting result in a "part" field"""
    zones = Internal.getZones(t)
    for z in zones:
        a1 = C.getAllFields(z, 'nodes')[0]
        a2 = C.getAllFields(z, 'centers')[0] # must contain "part" field
        Transform.transform.splitNGon2(a1, a2, N, N2, shift)
        C.setFields([a2], z, 'centers', writeDim=False)
    return None

def stick(t, psurf, stickBCName='FamilySpecified:stick', nitSmooth=0):
    """Stick a mesh on a surface."""
    tp = Internal.copyRef(t)
    _stick(tp, psurf, stickBCName, nitSmooth)
    return tp

def _stick(t, psurf, stickBCName='FamilySpecified:stick', nitSmooth=0):
    """Stick a mesh on a surface."""
    from . import Stick
    Stick._stick(t, psurf, stickBCName, nitSmooth)
    return None
