"""Grid generation module.
"""
__version__ = '4.0'
__author__ = "Stephanie Peron, Sam Landier, Christophe Benoit, Gaelle Jeanfaivre, Pascal Raud, Luis Bernardos"

from . import generator
import numpy
import Converter as C

from .TFIs import TFITri, TFIO, TFIHalfO, TFIMono, TFIStar, TFIStar2, mergeEdges

__all__ = ['cart', 'cartr1', 'cartr2', 'cartHexa', 'cartTetra', 'cartPenta',
    'cartPyra', 'cartNGon', 'cylinder', 'cylinder2', 'cylinder3', 'delaunay',
    'checkDelaunay', 'constrainedDelaunay', 'check', 'bbox', 'BB',
    'barycenter', 'CEBBIntersection', 'bboxIntersection', 'checkPointInCEBB',
    'enforceX', 'enforceY', 'enforceZ', 'enforcePlusX', 'enforcePlusY',
    'enforcePlusZ', 'enforceMoinsX', 'enforceMoinsY', 'enforceMoinsZ',
    'enforceLine', 'enforcePoint', 'enforceCurvature', 'enforceCurvature2',
    'addPointInDistribution', 'map', 'map1d', 'map1dpl', 'map2d',
    'mapCurvature', 'refine', 'defineSizeMapForMMGs', 'mmgs', 'densify',
    'hyper2D', 'hyper2D2', 'hyper2D3', 'hyper2D4', 'close', 'closeLegacy', 'zip',
    'pointedHat', 'stitchedHat', 'plaster', 'selectInsideElts', 'grow', 'stack', 
    'TFI', 
    'TFITri', 'TFIO', 'TFIHalfO', 'TFIMono', 'TFIStar', 'TFIStar2', 'mergeEdges',
    'TTM', 'bboxOfCells', 'getCellPlanarity', 'getVolumeMap',
    'getCellCenters', 'getFaceCentersAndAreas', 'getNormalMap',
    'getSmoothNormalMap', 'getEdgeRatio', 'getMaxLength', 'collarMesh',
    'surfaceWalk', 'buildExtension', 'getCircumCircleMap', 'getInCircleMap',
    'addNormalLayers', 'gencartmb', 'mapSplit', 'T3mesher2D', 'tetraMesher',
    'fittingPlaster', 'gapfixer', 'gapsmanager', 'front2Hexa', 'front2Struct',
    'snapFront', 'snapSharpEdges', 'fillWithStruct', 'octree2Struct',
    'cutOctant', 'octree', 'conformOctree3', 'adaptOctree', 'expandLayer',
    'forceMatch', '_forceMatch', 'getOrthogonalityMap', 'getRegularityMap',
    'getAngleRegularityMap', 'getTriQualityMap', 'getTriQualityStat',
    'quad2Pyra', 'extendCartGrids', 'checkMesh']

def cart(Xo, H, N, api=1):
    """Create a cartesian mesh defined by a structured array.
    Usage: cart((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    return generator.cart(Xo, H, N, api)

def cartr1(Xo, H, R, N, doubleLeft=(0,0,0), doubleRight=(0,0,0), api=1):
    """Create a structured cartesian mesh with geometric distribution.
    Usage: cartt1((xo,yo,zo), (hi,hj,hk), (ri,rj,rk), (ni,nj,nk))"""
    return generator.cartr1(Xo, H, R, N, doubleLeft, doubleRight, api)

def cartr2(Xo, H, R, Xf, doubleLeft=(0,0,0), doubleRight=(0,0,0), api=1, skeleton=False):
    """Create a structured cartesian mesh with geometric distribution fixing last point.
    Usage: cartr2((xo,yo,zo), (hi,hj,hk), (ri,rj,rk), (xf,yf,zf))"""
    return generator.cartr2(Xo, H, R, Xf, doubleLeft, doubleRight, api, skeleton)

def cartHexa(Xo, H, N, api=1):
    """Create a cartesian mesh defined by an hexaedrical array.
    Usage: cartHexa((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    return generator.cartHexa(Xo, H, N, api)

def cartTetra(Xo, H, N, api=1):
    """Create a cartesian mesh defined by a tetraedrical array.
    Usage: cartTetra((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    return generator.cartTetra(Xo, H, N, api)

def cartPenta(Xo, H, N, api=1):
    """Create a cartesian mesh defined by a prismatic array.
    Usage: cartPenta((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    return generator.cartPenta(Xo, H, N, api)

def cartPyra(Xo, H, N, api=1):
    """Create a cartesian mesh defined by a pyramidal array.
    Usage: cartPyra((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    return generator.cartPyra(Xo, H, N, api)

def cartNGon(Xo, H, N, api=1):
    """Create a cartesian mesh defined by a NGON array.
    Usage: cartNGon((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    return generator.cartNGon(Xo, H, N, api)

def cylinder(Xo, R1, R2, tetas, tetae, H, N):
    """Create a portion of regular cylindrical grid.
    Usage: cylinder((xo,yo,zo), R1, R2, tetas, tetae, H, (ni,nj,nk))"""
    return generator.cylinder(Xo, R1, R2, tetas, tetae, H, N)

def cylinder2(Xo, R1, R2, tetas, tetae, H, arrayR, arrayT, arrayZ):
    """Create a portion of cylindrical grid.
    Usage: cylinder2((xo,yo,zo), R1, R2, tetas, tetae, H, arrayR. arrayT, arrayZ)"""
    return generator.cylinder2(Xo, R1, R2, tetas, tetae,
                               H, arrayR, arrayT, arrayZ)

def cylinder3(arrayxz, tetas, tetae, arrayteta):
    """Create a portion of cylindrical grid.
    Usage: cylinder3(arrayxz, tetas, tetae, arrayteta)"""
    return generator.cylinder3(arrayxz, tetas, tetae, arrayteta)

def delaunay(array, tol=1.e-10, keepBB=0):
    """Create a delaunay mesh given a set of points defined by array.
    Usage: delaunay(array, tol, keepBB)"""
    array = close(array, tol)
    return generator.delaunay(array, keepBB)

def checkDelaunay(contour, tri):
    """Check if the Delaunay triangulation defined by tri is inside the contour.
    Usage: checkDelaunay(contour, tri)"""
    return generator.checkDelaunay(contour, tri)

def constrainedDelaunay(cont0, tol=1.e-10, keepBB=0):
    """Create a constrained-Delaunay mesh starting from a BAR-array defining
    the contour.
    Usage: constrainedDelaunay(array, tol, keepBB)"""
    tri = delaunay(cont0, tol, keepBB)
    cont = generator.checkDelaunay(cont0, tri)
    cntmax = 100
    cnt = 0
    while (cont0[2][0].shape[0] != cont[2][0].shape[0] and cnt < cntmax):
        cont0 = cont
        tri = delaunay(cont0, tol, keepBB)
        cont = generator.checkDelaunay(cont0, tri)
        cnt += 1
    return tri
    
def check(array):
    """Check a mesh for regularity, orthogonality...
    Usage: check(array)"""
    generator.check(array)

def bbox(arrays):
    """Returns the bounding box of a list of arrays.
    Usage: bbox(arrays)"""
    import KCore
    if len(arrays) == 0: return [1.e256, 1.e256, 1.e256, -1.e256, -1.e256, -1.e256]
    if not isinstance(arrays[0], list): ars = [arrays]
    else: ars = arrays
    xmin = 1.e256; ymin = 1.e256; zmin = 1.e256
    xmax =-1.e256; ymax =-1.e256; zmax =-1.e256
    for a in ars:
        varx = KCore.isNamePresent(a, 'CoordinateX')
        if varx == -1:
            varx = KCore.isNamePresent(a, 'x')
            if varx == -1:
                raise ValueError("bbox: x-coordinate not found.")
            else: varx = 'x'
        else: varx = 'CoordinateX'

        vary = KCore.isNamePresent(a, 'CoordinateY')
        if vary == -1:
            vary = KCore.isNamePresent(a, 'y')
            if vary == -1:
                raise ValueError("bbox: y-coordinate not found.")
            else: vary = 'y'
        else: vary = 'CoordinateY'

        varz = KCore.isNamePresent(a, 'CoordinateZ')
        if varz == -1:
            varz = KCore.isNamePresent(a, 'z')
            if varz == -1 :
                raise ValueError("bbox: z-coordinate not found.")
            else: varz = 'z'
        else: varz = 'CoordinateZ'  
            
        xmin0 = C.getMinValue(a, varx) 
        ymin0 = C.getMinValue(a, vary)
        zmin0 = C.getMinValue(a, varz)
        xmax0 = C.getMaxValue(a, varx)
        ymax0 = C.getMaxValue(a, vary)
        zmax0 = C.getMaxValue(a, varz)

        if xmin0 < xmin: xmin = xmin0
        if ymin0 < ymin: ymin = ymin0
        if zmin0 < zmin: zmin = zmin0
        if xmax0 > xmax: xmax = xmax0
        if ymax0 > ymax: ymax = ymax0
        if zmax0 > zmax: zmax = zmax0

    return [xmin, ymin, zmin, xmax, ymax, zmax]

def BB(array, method='AABB', weighting=0, tol=0.):
    """Return the axis-aligned or oriented bounding box of an array as an array.
    Usage: b = BB(a, method='AABB', weighting=0, tol=0.)"""
    if isinstance(array[0], list):
        out = []
        for a in array:
            if method == 'AABB':  # Computes AABB
                sbb = bbox(a)
                ar = C.array('x,y,z', 2, 2, 2)
                C.setValue(ar, (1,1,1), [sbb[0]-tol, sbb[1]-tol, sbb[2]-tol])
                C.setValue(ar, (2,1,1), [sbb[3]+tol, sbb[1]-tol, sbb[2]-tol])
                C.setValue(ar, (1,2,1), [sbb[0]-tol, sbb[4]+tol, sbb[2]-tol])
                C.setValue(ar, (1,1,2), [sbb[0]-tol, sbb[1]-tol, sbb[5]+tol])
                C.setValue(ar, (2,2,1), [sbb[3]+tol, sbb[4]+tol, sbb[2]-tol])
                C.setValue(ar, (1,2,2), [sbb[0]-tol, sbb[4]+tol, sbb[5]+tol])
                C.setValue(ar, (2,1,2), [sbb[3]+tol, sbb[1]-tol, sbb[5]+tol])
                C.setValue(ar, (2,2,2), [sbb[3]+tol, sbb[4]+tol, sbb[5]+tol])
                out.append(ar)
            elif method == 'OBB':  # Computes OBB
                out.append(generator.obbox(a, weighting))
            else:
                print('BB: Warning, method=%s not implemented, making an OBB.'%method)
                out.append(generator.obbox(a, weighting))
        return out
    else:
        if method == 'AABB':  # Computes AABB
            sbb = bbox(array)
            ar = C.array('x,y,z', 2, 2, 2)
            C.setValue(ar, (1,1,1), [sbb[0]-tol, sbb[1]-tol, sbb[2]-tol])
            C.setValue(ar, (2,1,1), [sbb[3]+tol, sbb[1]-tol, sbb[2]-tol])
            C.setValue(ar, (1,2,1), [sbb[0]-tol, sbb[4]+tol, sbb[2]-tol])
            C.setValue(ar, (1,1,2), [sbb[0]-tol, sbb[1]-tol, sbb[5]+tol])
            C.setValue(ar, (2,2,1), [sbb[3]+tol, sbb[4]+tol, sbb[2]-tol])
            C.setValue(ar, (1,2,2), [sbb[0]-tol, sbb[4]+tol, sbb[5]+tol])
            C.setValue(ar, (2,1,2), [sbb[3]+tol, sbb[1]-tol, sbb[5]+tol])
            C.setValue(ar, (2,2,2), [sbb[3]+tol, sbb[4]+tol, sbb[5]+tol])
        elif method == 'OBB':  # Computes OBB
            ar = generator.obbox(array, weighting)
        else:
            print('BB: Warning, method=%s not implemented, making an OBB.'%method)
            ar = generator.obbox(array, weighting)
        return ar       
     
def barycenter(array, weight=None):
    """Get the barycenter of an array.
    Usage: barycenter(a, w)"""
    if isinstance(array[0], list):
        N = 0; xb = 0; yb = 0; zb = 0
        for i, a in enumerate(array):
            if weight is not None:
                X = generator.barycenter(a, weight[i])
            else: X = generator.barycenter(a, None)
            if isinstance(a[1], list): # array2/3
               n = a[1][0].size
            else: n = a[1].shape[1]
            xb += X[0]*n
            yb += X[1]*n
            zb += X[2]*n
            N += n
        xb = xb / N; yb = yb / N; zb = zb / N
        return [xb, yb, zb]
    else:
        return generator.barycenter(array, weight)
    
def CEBBIntersection(array1, array2, tol=1.e-10):
    """Get the Cartesian Elements bounding box intersection of 2 arrays."""
    return generator.CEBBIntersection(array1, array2, tol)

def bboxIntersection(array1, array2, tol=1.e-6, isBB=False, method='AABB'):
    """Return the intersection of bounding boxes of 2 arrays."""
    if not isBB:
        array1 = BB(array1, method='OBB')
        array2 = BB(array2, method='OBB')
    if method == 'AABB':  # Computes the intersection between 2 AABB
        return generator.bboxIntersection(array1, array2, tol)
    elif method == 'OBB':  # Computes the intersection between 2 OBB
        return generator.obboxIntersection(array1, array2)
    elif method == 'AABBOBB':  # Computes the intersection between an AABB and an OBB
        return generator.crossIntersection(array1, array2)        
    else:
        print('Warning: bboxIntersection: method %s not implemented, switching to AABB.'%method)
        return generator.bboxIntersection(array1, array2)

def checkPointInCEBB(array, P):
    """Check if point P is in the Cartesian Elements Bounding Box of array."""
    return generator.checkPointInCEBB(array, P)

def enforceX(array, x0, enforcedh, N, add=0, verbose=True):
    """Enforce a x0-centered line in a distribution defined by an array.
    Usage: enforceX(array, x0, enforcedh, supp, add) -or-
    Usage: enforceX(array, x0, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforceX(array, x0, enforcedh, N, verbose)
    else:
        return generator.enforce(array, "enforceX", x0, enforcedh, N, add)
    
def enforceY(array, y0, enforcedh, N, add=0, verbose=True):
    """Enforce a j line in a distribution defined by an array.
    Usage: enforceY(array, y0, enforcedh, supp, add) -or-
    Usage: enforceY(array, y0, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforceY(array, y0, enforcedh, N, verbose)
    else:
        return generator.enforce(array, "enforceY", y0, enforcedh, N, add)
    
def enforceZ(array, z0, enforcedh, N, add=0):
    """Enforce a k line in a distribution defined by an array.
    Usage: enforceZ(array, z0, enforcedh, supp, add)"""
    return generator.enforce(array, "enforceZ", z0, enforcedh, N, add)
    
def enforcePlusX(array, enforcedh, N, add=0, verbose=True):
    """Enforce the first X-line in a distribution defined by an array.
    (one sided distribution, right).
    Usage: enforcePlusX(array, enforcedh, supp, add) -or-
    Usage: enforcePlusX(array, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforcePlusX(array, enforcedh, N, verbose)
    else:
        return generator.enforce(array, "enforcePlusX", 0., enforcedh, N, add)
    
def enforcePlusY(array, enforcedh, N, add=0, verbose=True):
    """Enforce a j line in a distribution  defined by an array.
    (one sided distribution, top).
    Usage: enforcePlusY(array, enforcedh, supp, add) -or-
    Usage: enforcePlusY(array, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforcePlusY(array, enforcedh, N, verbose)
    else:
        return generator.enforce(array, "enforcePlusY", 0., enforcedh, N, add)

def enforcePlusZ(array, enforcedh, N, add=0):
    """Enforce a k line in a distribution  defined by an array.
    (one sided distribution, top).
    Usage: enforcePlusZ(array, enforcedh, supp, add)"""
    return generator.enforce(array, "enforcePlusZ", 0., enforcedh, N, add)
    
def enforceMoinsX(array, enforcedh, N, add=0, verbose=True):
    """Enforce the last X-line in a distribution (one sided, left).
    Usage: enforceMoinsX(array, enforcedh, supp, add) -or-
    Usage: enforceMoinsX(array, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforceMoinsX(array, enforcedh, N, verbose)
    else:
        return generator.enforce(array, "enforceMoinsX", 0., enforcedh, N, add)
    
def enforceMoinsY(array, enforcedh, N, add=0, verbose=True):
    """Enforce a j line in a distribution (one sided distribution, bottom).
    Usage: enforceMoinsY(array, enforcedh, supp, add) -or-
    Usage: enforceMoinsY(array, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforceMoinsY(array, enforcedh, N, verbose)
    else:
        return generator.enforce(array, "enforceMoinsY", 0., enforcedh, N, add)

def enforceMoinsZ(array, enforcedh, N, add=0):
    """Enforce a k line in a distribution (one sided distribution, bottom).
    Usage: enforceMoinsY(array, enforcedh, supp, add)"""
    return generator.enforce(array, "enforceMoinsZ", 0., enforcedh, N, add)

def enforceLine(array, arrayline, enforcedh, N):
    """Enforce a line in a distribution.
    Usage: enforceLine(array, arrayline, enforcedh, (supp,add))"""
    return generator.enforceLine(array, arrayline, enforcedh, N)
    
def getTanhDist__(Nx,dy1,dy2,ityp):
    """Private function. Wrapper around k6stretch function of StretchF.for.
    Returns an array of a 1D line between x=0.0 and x=1.0 of Nx points, with a
    first cell length of dy1 if ityp=1 or ityp=2 and a last cell height of dy2
    if ityp=2. The distribution law is an hyperbolic tangent. """
    return generator.tanhdist__( Nx, dy1, dy2, ityp)

def enforcePoint(array, x0):
    """Enforce a point in a distribution.
    Usage: enforcePoint(array, x0)"""
    return generator.enforcePoint(array, x0)

def enforceCurvature(arrayD, arrayC, power=0.5):
    """Enforce curvature of a curve in a distribution.
    The distribution is defined by arrayD, and the curvature by arrayC
    Usage: enforceCurvature(arrayD, arrayC, power)"""
    return generator.enforceCurvature(arrayD, arrayC, power)

#--------------------------------------------------------------
# enforce a 1D distribution wrt the curvature radius
# alpha is the factor of stretch for point of maximum curvature
#--------------------------------------------------------------
def enforceCurvature2(arrayD, arrayC, alpha=1.e-2):
    try: import Geom as D; import math; import KCore
    except: raise ImportError("enforceCurvature2: requires Converter and Geom modules.")
    
    tol = 1.e-12 # tolerance sur les pts confondues(=close)
    loop = 0 # le contour est il une boucle

    posx = KCore.isNamePresent(arrayC, 'x')
    if posx != -1: xt = C.extractVars(arrayC, ['x'])[1]
    else:
        posx = KCore.isNamePresent(arrayC, 'CoordinateX')
        if posx != -1: xt = C.extractVars(arrayC, ['CoordinateX'])[1]
        else: raise ValueError("enforceCurvature2: coordinates must be present in array.")
    posy = KCore.isNamePresent(arrayC, 'y')
    if posy != -1: yt = C.extractVars(arrayC, ['y'])[1]
    else:
        posy = KCore.isNamePresent(arrayC, 'CoordinateY')
        if posy != -1: yt = C.extractVars(arrayC, ['CoordinateY'])[1]
        else: raise ValueError("enforceCurvature2: coordinates must be present in array.")
    posz = KCore.isNamePresent(arrayC, 'z')
    if posz != -1: zt = C.extractVars(arrayC, ['z'])[1]
    else:
        posz = KCore.isNamePresent(arrayC, 'CoordinateZ')
        if posz != -1: zt = C.extractVars(arrayC, ['CoordinateZ'])[1]
        else: raise ValueError("enforceCurvature2: coordinates must be present in array.")
    
    nmax = xt.shape[1]-1
    if abs(xt[0,0]-xt[0,nmax])<tol and abs(yt[0,0]-yt[0,nmax])<tol and abs(zt[0,0]-zt[0,nmax])<tol: loop = 1
    
    rc = D.getCurvatureRadius(arrayC)
    ni = rc[1].shape[1]
    rcmin = C.getMinValue(rc,'radius'); rcmax = 1.
    for i in range(ni): rc[1][0,i] = min(1.,rc[1][0,i])
    rcmean = C.getMeanValue(rc, 'radius')
    dh = C.initVars(rc, 'dhloc', 1.); dh = C.extractVars(dh, ['dhloc'])
    coefa = math.log(alpha)/(rcmin-1.)
    for i in range(ni):
        rad = rc[1][0,i]
        if rad < 0.2*rcmean: dh[1][0,i] = math.exp(coefa*rc[1][0,i]-coefa)
    if loop == 1: rc[1][0,ni-1] = rc[1][0,0]
    minima = []; dht = dh[1]
    dhmax = C.getMaxValue(dh, 'dhloc'); dhmin = C.getMinValue(dh, 'dhloc')
    if abs(dhmax-dhmin) < 1.e-10: return arrayD
    
    # premier point
    if loop == 0:
        if dht[0,0] < dht[0,1]: minima.append(0)
    else: # loop = 1
        if dht[0,0]<dht[0,1] and dht[0,0]<dht[0,nmax-1]: minima.append(0)
    for i in range(1,ni-1):
        dhloc = dht[0,i]
        im1=i-1; ip1=i+1; dhm1 = dht[0,im1]; dhp1 = dht[0,ip1]
        if dhloc < dhm1 and dhloc<dhp1: minima.append(i)

    # dernier point
    if loop == 0 and dht[0,nmax]<dht[0,nmax-1]: minima.append(nmax)
    if minima == [] : return arrayD
    if len(minima) > 1: # verifier que 2 pts separes d'un pt ne sont pas minima
        minima0 = []
        m = minima[0]; mp1 = minima[1]; mm1 = m
        if mp1-m > 2: minima0.append(m)
        for nom in range(1,len(minima)-1):
            m = minima[nom]; mp1 = minima[nom+1]
            if m-mm1 > 2 and mp1-m > 2: minima0.append(m)
            elif m-mm1 <= 2 and mp1-m <= 2: minima0.append(m)
            elif m-mm1 == 2: minima0.append(m-1)
            mm1 = m
        # dernier pt
        m = mp1; mm1 = minima0[len(minima0)-1]
        if m-mm1 > 2: minima0.append(m) # le precedent si adjacent n a pas ete pris en compte
        minima = minima0
    if minima == []: return arrayD
    #---------------------
    # abscisse curviligne
    #---------------------
    posx = KCore.isNamePresent(arrayD, 'x')
    if posx!=-1: varx='x'
    else:
        posx = KCore.isNamePresent(arrayD, 'CoordinateX')
        if posx == -1: print('Warning: enforceCurvature2: x variable not found.'); return arrayD
        varx = 'CoordinateX'
    xs = C.getMinValue(arrayD, varx); xe = C.getMaxValue(arrayD, varx)
    s = D.getCurvilinearAbscissa(arrayC)[1]
    stot = s[0, ni-1]
    for i in range(ni): s[0,i] = (s[0,i]/stot)*(xe-xs)+xs

    # recherche du point a forcer dans distrib
    distrib = C.copy(arrayD)
    distrib[0] = 'x,y,z'
    for indc in minima:
        xt = C.extractVars(distrib, ['x'])[1]; nid = distrib[2]
        sc = s[0,indc]; ind = D.getDistantIndex(distrib, 1,sc-xs)
        x0 = xt[0,ind]; dhloc = dht[0,indc]
        if ind == 0:
            dhloc = dhloc*(xt[0,1]-x0); distrib = enforcePlusX(distrib, dhloc, (5,15))
            if loop == 1: distrib = enforceMoinsX(distrib, dhloc, (5,15))
        elif ind == xt.shape[1]-1: dhloc = dhloc*(x0-xt[0,nid-2]); distrib = enforceMoinsX(distrib, dhloc, (5,15))
        else: dhloc = dhloc*(xt[0,ind+1]-xt[0,ind]); distrib = enforceX(distrib,sc, dhloc, 5,15)
    return distrib

def addPointInDistribution(array, ind):
    """Add a point in a distribution defined by array.
    Usage: addPointInDistribution(array, ind)"""
    return generator.addPointInDistribution(array, ind)

def map(array, d, dir=0, h1=None, h2=None, isAvg=False, nAvg=2):
    """Map a distribution on a curve or a surface.
    Usage: map(array, d, dir, h1, h2, isAvg, nAvg)"""
    if len(d) == 5 and d[3] != 1 and d[4] == 1 and dir == 0:
        return map2d(array, d)
    elif len(d) == 5 and dir != 0:
        return map1dpl(array, d, dir, h1, h2, isAvg, nAvg)
    else: return map1d(array, d)
    
# map sur une courbe
def map1d(array, d):
    """Map on a curve."""
    return generator.map(array, d)

# map par lignes dans la direction dir
def map1dpl(array, d, dir, h1, h2, isAvg, pnts):
    try: import Transform as T
    except:
        raise ImportError("map: requires Transform and Converter modules.")

    islocationdependent = False
    if h1 is not None and h2 is not None:
        islocationdependent = True
        import numpy
        import Geom as D
        import Geom.MapEdge as MapE
        N = len(d[1][0])
    
    if dir == 2: m = T.reorder(array, (2,1,3))
    elif dir == 3: m = T.reorder(array, (3,2,1))
    elif dir == 1: m = array
    ni = m[2]; nj = m[3]; nk = m[4]; ndi = d[2]; ndi2 = ndi*nj
    a = C.array('x,y,z', ndi, nj, nk)

    if islocationdependent:
        for k in range(nk):
            for j in range(nj):                
                l = T.subzone(m, (1,j+1,k+1), (ni,j+1,k+1))
                h1_local = h1
                h2_local = h2                
                if h1_local<0:h1_local=numpy.sqrt((l[1][0][0]   -l[1][0][1]   )**2+
                                                  (l[1][1][0]   -l[1][1][1]   )**2+
                                                  (l[1][2][0]   -l[1][2][1]   )**2)
                if h2_local<0:h2_local=numpy.sqrt((l[1][0][ni-1]-l[1][0][ni-2])**2+
                                                  (l[1][1][ni-1]-l[1][1][ni-2])**2+
                                                  (l[1][2][ni-1]-l[1][2][ni-2])**2)
                length_local = D.getLength(l)                
                d = MapE.buildDistrib(h1_local/length_local,h2_local/length_local,N+1)

                if isAvg:
                    d_sum=d[1][0]
                    if j>pnts and j<nj-(pnts+1):
                        for i in range(1,pnts+1):
                            d_local_val = d_local(m,j-i,k,ni,h1,h2,N)
                            d_sum      += d_local_val[1][0]
                            d_local_val = d_local(m,j+i,k,ni,h1,h2,N)
                            d_sum      += d_local_val[1][0]
                        d[1][0] = d_sum/(2*pnts+1)

                ind = j*ndi+k*ndi2
                am = map1d(l, d)
                a[1][:,ind:ndi+ind] = am[1][:,0:ndi]
    else:
        for k in range(nk):
            for j in range(nj):
                l = T.subzone(m, (1,j+1,k+1), (ni,j+1,k+1))
                am = map1d(l, d)
                ind = j*ndi+k*ndi2
                a[1][:,ind:ndi+ind] = am[1][:,0:ndi]

    if dir == 2: a = T.reorder(a, (2,1,3))
    elif dir == 3: a = T.reorder(a, (3,2,1))
    return a

def d_local(m,j,k,ni,h1,h2,N):
    import Transform as T
    import numpy
    import Geom as D
    import Geom.MapEdge as MapE
    l = T.subzone(m, (1,j+1,k+1), (ni,j+1,k+1))
    h1_local = h1
    h2_local = h2                
    if h1_local<0:h1_local=numpy.sqrt((l[1][0][0]   -l[1][0][1]   )**2+
                                      (l[1][1][0]   -l[1][1][1]   )**2+
                                      (l[1][2][0]   -l[1][2][1]   )**2)
    if h2_local<0:h2_local=numpy.sqrt((l[1][0][ni-1]-l[1][0][ni-2])**2+
                                      (l[1][1][ni-1]-l[1][1][ni-2])**2+
                                      (l[1][2][ni-1]-l[1][2][ni-2])**2)
    length_local = D.getLength(l)                
    
    d = MapE.buildDistrib(h1_local/length_local,h2_local/length_local,N+1)
    return d

# map sur une surface
def map2d(array, d):
    try: import Transform as T
    except:
        raise ImportError("map: requires Transform and Converter modules.")
    di = T.subzone(d, (1,1,1), (d[2],1,1)); ndi = di[2]
    b = T.rotate(d, (0,0,0), (0,0,1), -90.)
    b = T.reorder(b, (2,1,3))
    dj = T.subzone(b, (1,1,1), (b[2],1,1)); ndj = dj[2]
    # Passage en i
    ni = array[2]; nj = array[3]
    m = C.array('x,y,z', ndi, nj, 1)
    for j in range(nj):
        a = T.subzone(array, (1,j+1,1), (ni,j+1,1))
        am = map1d(a, di)
        for i in range(ndi):
            v = C.getValue(am, (i+1,1,1))
            C.setValue(m, (i+1,j+1,1), v)
    m = T.reorder(m, (2,1,3))
    p = C.array('x,y,z', ndi, ndj, 1)
    for i in range(ndi):
        a = T.subzone(m, (1,i+1,1), (nj,i+1,1))
        am = map1d(a, dj)
        for j in range(ndj):
            v = C.getValue(am, (j+1,1,1))
            C.setValue(p, (i+1,j+1,1), v)
    return p

# Map un maillage structure suivant sa courbure
def mapCurvature(array, N, power, dir):
    """Remesh with a step proportional to curvature."""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(mapCurvature___(i, N, power, dir))
        return b
    else:
        return mapCurvature___(array, N, power, dir)

def mapCurvature___(array, N, power, dir):
    try: import Transform as T
    except: 
        raise ImportError("mapCurvature: requires Transform, Converter modules.")
    if dir == 2: m = T.reorder(array, (2,1,3))
    elif dir == 3: m = T.reorder(array, (3,2,1))
    elif dir == 1: m = array
    ni = m[2]; nj = m[3]; nk = m[4]
    N2 = N*nj
    a = C.array('x,y,z', N, nj, nk)
    d = cart((0,0,0), (1./(N-1),1,1), (N,1,1))
    for k in range(nk):
        for j in range(nj):
            l = T.subzone(m, (1,j+1,k+1), (ni,j+1,k+1))
            e = enforceCurvature(a, l, power)
            am = map1d(l, e)
            ind = j*N+k*N2
            a[1][:,ind:N+ind] = am[1][:,0:N]
    if dir == 2: a = T.reorder(a, (2,1,3))
    elif dir == 3: a = T.reorder(a, (3,2,1))
    return a

#==============================================================================
# Raffine un maillage structure (dir=0,1,2,3)
#==============================================================================
def refine(array, power, dir=0):
    """Refine a mesh of power power along all the directions or on a specified one. Usage: refine(a,power,dir)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            if dir == 0:
                for dir0 in [1,2,3]: i = refinePerDir__(i, power, dir0)
                b.append(i)
            else: b.append(refinePerDir__(i, power, dir))
        return b
    else: 
        if dir == 0:
            a2 = refinePerDir__(array, power, 1)
            a2 = refinePerDir__(a2, power, 2)
            a2 = refinePerDir__(a2, power, 3)
            return a2
        else: return refinePerDir__(array, power, dir)

def refinePerDir__(a, power, dir):
    if a[2] == 1 and dir == 1: return a
    if a[3] == 1 and dir == 2: return a
    if a[4] == 1 and dir == 3: return a
    if not isinstance(power,float): factor = float(power)
    else: factor = power
    try:
        import Transform as T
        import Geom as D; import Post as P
    except: raise ImportError("refine: requires Transform, Converter, Geom, Post modules.")
    if dir != 1 and dir != 2 and dir != 3: raise ValueError("refine: dir must be 1, 2 or 3.")

    array = C.copy(a)
    if dir == 2: array = T.reorder(array,(2,1,3))
    elif dir == 3: array = T.reorder(array,(3,2,1))
    
    ni = array[2]; nj = array[3]; nk = array[4]
    NI = int(round(factor*(ni-1)+1))

    NJ = nj; NK = nk

    distribI = C.array('x,y,z', NI, NJ, NK)  
    vvI = 1./max(NI-1, 1.e-12)
    vvJ = 1./max(NJ-1, 1.e-12)
    vvK = 1./max(NK-1, 1.e-12)
    NINJ = NI*NJ; ninj = ni*nj
    # Distribution fine a interpoler: cartesienne a pas constant par direction
    for k in range(NK):
        for j in range(NJ):
            for i in range(NI): 
                ind = i + j*NI + k*NINJ
                distribI[1][0,ind] = i*vvI
                distribI[1][1,ind] = j*vvJ
                distribI[1][2,ind] = k*vvK                

    # Calcul de la distribution en i
    dhi =  1./max(ni-1, 1.e-12)
    dhj =  1./max(nj-1, 1.e-12)
    dhk =  1./max(nk-1, 1.e-12)
    coord = cart((0,0,0),(dhi,dhj,dhk),(ni,nj,nk))
    coord = C.addVars(coord,['s'])
    for k in range(nk):
        for j in range(nj):
            l = T.subzone(array, (1,j+1,k+1), (ni,j+1,k+1))
            cu = D.getCurvilinearAbscissa(l)
            l0 = T.subzone(coord, (1,j+1,k+1), (ni,j+1,k+1))
            cu = C.addVars([l0,cu])
            ind = j*ni+k*ninj
            coord[1][:,ind:ni+ind] = cu[1]

    distribI = P.extractMesh([coord], distribI, order=2)
    distribI = C.extractVars(distribI, ['s','y','z']); distribI[0] = 'x,y,z'
    
    aout = C.array('x,y,z', NI, nj, nk)
    for k in range(nk):
        for j in range(nj):
            l = T.subzone(array, (1,j+1,k+1), (ni,j+1,k+1))
            dis = T.subzone(distribI,(1,j+1,k+1), (NI,j+1,k+1))
            am = map1d(l, dis)
            ind = j*NI+k*NINJ
            a1 = aout[1]; am1 = am[1]
            a1[:,ind:NI+ind] = am1[:,0:NI]
    
    if dir == 2: return T.reorder(aout, (2,1,3))
    elif dir == 3: return T.reorder(aout, (3,2,1))
    else: return aout

def defineSizeMapForMMGs(array, hmax, sizeConstraints):
    import KCore; import Generator; import Transform
    if hmax > 0: array = C.initVars(array, 'sizemap=%f'%hmax)
    else: 
        vol = Generator.getVolumeMap(array)
        vol = C.initVars(vol, '{vol}=(1.15*{vol})**0.5')
        vol = C.center2Node(vol)
        vol[0] = 'sizemap'
        array = C.addVars([array, vol])
    pos = KCore.isNamePresent(array, 'sizemap')

    szcs = C.convertArray2Hexa(sizeConstraints)
    c = Transform.join(szcs)
    v = Generator.getVolumeMap(c)
    v = C.center2Node(v) # devrait etre max
    c = C.addVars([c,v])
    hook = C.createHook(c, function='nodes')
    ret = C.nearestNodes(hook, array)
    n = ret[0]; d = ret[1] 
    pt = array[1]
    alpha = numpy.empty(d.size, dtype=numpy.float64)
    if isinstance(pt, list): alpha[:] = 0.3*d[:] / pt[pos,:]
    else: alpha[:] = 0.3*d[:] / pt[pos][:]
    alpha[:] = (alpha[:] > 1.)*1.+(alpha[:] <= 1.)*alpha[:]
    if isinstance(pt, list): pt[pos,:] = alpha[:]*pt[pos,:]+(1.-alpha[:])*v[1][0,n[:]-1]
    else: pt[pos][:] = alpha[:]*pt[pos][:]+(1.-alpha[:])*v[1][0,n[:]-1]
    #C.convertArrays2File(array, 'array.plt')
    return array

# Remaille une surface avec mmgs
def mmgs(array, ridgeAngle=45., hmin=0., hmax=0., hausd=0.01, grow=1.1, 
         anisotropy=0, optim=0, fixedConstraints=[], sizeConstraints=[]):
    """Surface remeshing using MMGS."""
    if isinstance(array[0], list):
        l = []
        for i in array:
            if fixedConstraints != []:
                fixedNodes = []; fixedEdges = []
                hook = C.createHook(i, function='nodes')
                for c in fixedConstraints:
                    loc = C.nearestNodes(hook, c)[0]
                    fixedNodes.append(loc)
                    if c[3] == 'BAR':
                        pt = c[2].copy()
                        pt[:,:] = loc[pt[:,:]-1]
                        fixedEdges.append(pt)
            else: fixedNodes = None; fixedEdges = None
            if sizeConstraints != []:
                i = defineSizeMapForMMGs(i, hmax, sizeConstraints)
                hmaxl = 1000.
            else: hmaxl = hmax
            l.append(generator.mmgs(i, ridgeAngle, hmin, hmaxl, hausd,
                                    grow, anisotropy, optim, fixedNodes, fixedEdges))
        return l
    else:
        if fixedConstraints != []:
            fixedNodes = []; fixedEdges = []
            hook = C.createHook(array, function='nodes')
            for c in fixedConstraints:
                loc = C.nearestNodes(hook, c)[0]
                fixedNodes.append(loc)
                if c[3] == 'BAR':
                    pt = c[2].copy()
                    pt[:,:] = loc[pt[:,:]-1]
                    fixedEdges.append(pt)
        else: fixedNodes = None; fixedEdges = None
        if sizeConstraints != []:
            array = defineSizeMapForMMGs(array, hmax, sizeConstraints)
            hmaxl = 1000.
        else: hmaxl = hmax
        return generator.mmgs(array, ridgeAngle, hmin, hmaxl, hausd, 
                              grow, anisotropy, optim, fixedNodes, fixedEdges)

#==============================================================================
# Densifie le maillage d'un i-array
# IN: array: i-array
# IN: h: pas de discretisation
# OUT: i-array avec la nouvelle discretisation
#==============================================================================
def densify(array, h):
    """Densify a mesh."""
    if isinstance(array[0], list):
        l = []
        for i in array:
            l.append(generator.densify(i, h))
        return l
    else:
        return generator.densify(array, h)

def hyper2D(array, arrayd, type, 
            eta_start=10, eta_end=-1, beta=0.0):
    """Generate an hyperbolic mesh. 
    Usage: hyper2D(array, arrayd, type)"""
    return generator.hyper2D(array, arrayd, type, eta_start, eta_end, beta)

def hyper2D2(array, arrayd, type, alpha):
    """Generate an hyperbolic mesh with a constant alpha angle.
    Usage: hyper2D2(array, arrayd, type, alpha)"""
    return generator.hyper2D2(array, arrayd, type, alpha)

def hyper2D3(array, arrayd, type, alpha1, alpha2):
    """Generate an hyperbolic mesh with boundary alpha angles.
    Usage: hyper2D3(array, arrayd, type, alpha1, alpha2)"""
    return generator.hyper2D3( array, arrayd, type, alpha1, alpha2)
    
def hyper2D4(array, arrayd, type):
    """Generate an hyperbolic mesh.
    Usage: hyper2D4(array, arrayd, type)"""
    return generator.hyper2D4(array, arrayd, type)

def closeLegacy(array, tol=1.e-12, suppressDegeneratedNGons=False):
    """Close an unstructured mesh defined by an array gathering points closer than tol.
    Usage: close(array, tol)"""
    if isinstance(array[0], list):
        out = []
        for a in array:

            if len(a)==5: # merge intra-borders (C-type meshes)
                outl = generator.closeBorders([a], [], tol)[0]
            else:
                outl = generator.closeMeshLegacy(a, tol, suppressDegeneratedNGons)
            out.append(outl)
        return out
    else:
        return generator.closeMeshLegacy(array, tol, suppressDegeneratedNGons)
        
def close(array, tol=1.e-12, rmOverlappingPts=True, rmOrphanPts=True,
          rmDuplicatedFaces=True, rmDuplicatedElts=True,
          rmDegeneratedFaces=True, rmDegeneratedElts=True,
          indices=None):
    """Close an unstructured mesh defined by an array gathering points closer than tol.
    Usage: close(array, tol)"""
    exportIndirPts = False
    if isinstance(indices, list) and not indices: exportIndirPts = True
    if isinstance(array[0], list):
        out = []
        for a in array:
            indirl = None
            if len(a) == 5: # merge intra-borders (C-type meshes)
                outl = generator.closeBorders([a], [], tol)[0]
            else:
                outl = generator.closeMesh(a, tol, rmOverlappingPts,
                                           rmOrphanPts, rmDuplicatedFaces,
                                           rmDuplicatedElts, rmDegeneratedFaces,
                                           rmDegeneratedElts, exportIndirPts)
                if exportIndirPts: outl, indirl = outl
            out.append(outl)
            if exportIndirPts: indices.append(indirl)
        return out
    else:
        out = generator.closeMesh(array, tol, rmOverlappingPts,
                                  rmOrphanPts, rmDuplicatedFaces,
                                  rmDuplicatedElts, rmDegeneratedFaces,
                                  rmDegeneratedElts, exportIndirPts)
        if exportIndirPts:
            out, indirl = out
            indices.append(indirl)
        return out

def zip(array, tol=1e-12):
    """Zip a set of meshes defined by gathering exterior points closer than tol.
    Usage: zip(array, tol)"""
    if isinstance(array[0], list):
        extFaces = []
        try: 
            import Post as P
            for a in array: 
                if len(a) == 4: extFaces.append(P.exteriorFaces(a))
        except: pass
        return generator.closeBorders(array, extFaces, tol)
    else:
        return generator.closeBorders([array], [], tol)[0]

def pointedHat(array, coord):
    """Create a structured surface defined by a contour and a point (x,y,z).
    Usage: pointedHat(array, (x,y,z))"""
    return generator.pointedHat(array, coord)

def stitchedHat(array, offset, tol=1.e-6, tol2=1.e-5):
    """Create a structured surface defined by a contour and an offset (dx,dy,dz).
    Usage: stitchedHat(array, (dx,dy,dz))"""
    try: import Transform as T
    except: return generator.stitchedHat(array, tol)
    c = generator.stitchedHat(array, tol)
    t = T.subzone(c, (1,2,1),(c[2],2,1))
    t = close(t, tol2)
    t = T.translate(t, offset)
    c = T.patch(c, t, (1,2,1))
    return c

def plaster(contours, surfaces, side=0):
    """Create a sticky plaster around contours surrounded by surfaces.
    Usage: plaster(contours, surfaces)"""
    ni = 100; nj = 100
    try: import Transform as T
    except: raise ImportError("plaster: requires Converter, Transform modules.")

    c = C.convertArray2Tetra(contours); c = T.join(c)
    s = C.convertArray2Tetra(surfaces); s = T.join(s)
   
    bb = bbox(contours)
    bb[0] = bb[0] - 1.e-2; bb[1] = bb[1] - 1.e-2; bb[2] = bb[2] - 1.e-2
    bb[3] = bb[3] + 1.e-2; bb[4] = bb[4] + 1.e-2; bb[5] = bb[5] + 1.e-2
    lx = bb[3]-bb[0]; ly = bb[4]-bb[1]; lz = bb[5]-bb[2]
    s1 = lx*ly; s2 = lx*lz; s3 = ly*lz
    if s1 >= s2 and s1 >= s3:
        if side == 0:
            p = cart( (bb[0],bb[1],bb[2]), (lx/(ni-1),ly/(nj-1),1), (ni,nj,1) )
        else:
            p = cart( (bb[0],bb[1],bb[5]), (lx/(ni-1),ly/(nj-1),1), (ni,nj,1) )
        dir = (0,0,1)
    elif s2 >= s1 and s2 >= s3:
        if side == 0:
            p = cart( (bb[0],bb[1],bb[2]), (lx/(ni-1),1,lz/(nj-1)), (ni,1,nj) )
        else:
            p = cart( (bb[0],bb[4],bb[2]), (lx/(ni-1),1,lz/(nj-1)), (ni,1,nj) )
        dir = (0,1,0)
    else:
        if side == 0:
            p = cart( (bb[0],bb[1],bb[2]), (1,ly/(ni-1),lz/(nj-1)), (1,ni,nj) )
        else:
            p = cart( (bb[3],bb[1],bb[2]), (1,ly/(ni-1),lz/(nj-1)), (1,ni,nj) )
        dir = (1,0,0)
    p = T.projectDir(p, surfaces, dir, smooth=1)
    return p
    
def selectInsideElts(array, curvesList):
    """Select elements whose center is in the surface delimited by curves.
    Usage: selectInsideElts(array, curvesList)"""
    return generator.selectInsideElts(array, curvesList)

def grow(array, vector):
    """Grow a surface array of one layer by deplacing points of vector.
    Usage: grow(array, vector)"""
    return generator.grow(array, vector)

def stack(array1, array2=None):
    """Stack two meshes (with same nixnj) into a single mesh.
    Usage: stack(array1, array2)"""
    if array2 is not None: return generator.stack([array1, array2])
    else: return generator.stack(array1)

def TFI(arrays):
    """Generate a transfinite interpolation mesh from boundaries.
    Usage: TFI(arrays)"""
    return generator.TFI(arrays)

def TTM(array, niter=100):
    """Smooth a mesh with Thompson-Mastin elliptic generator.
    Usage: TTM(array, niter)"""
    return generator.TTM(array, niter)

def bboxOfCells(array):
    """Return the bounding box of all cells of an array.
    Usage: getBBoxOfCells(array)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            b.append(generator.bboxOfCells(i))
        return b
    else:
        return generator.bboxOfCells(array)

def getCellPlanarity(array):
    """Return the cell planarity of a surface mesh in an array.
    Usage: getCellPlanarity(array)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            b.append(generator.getCellPlanarity(i))
        return b
    else:
        return generator.getCellPlanarity(array)

def getVolumeMap(array, method=0):
    """Return the volume map in an array.
    Usage: getVolumeMap(array, method)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            b.append(generator.getVolumeMap(i, method))
        return b
    else:
        return generator.getVolumeMap(array, method)

def getFaceCentersAndAreas(array):
    """Return the face centers and areas in an NGon array.
    Usage: getFaceCentersAndAreas(array)"""
    return generator.getFaceCentersAndAreas(array)

def getCellCenters(array, fc, fa, own, nei):
    """Return the cell centers in an NGon array.
    Usage: getCellCenters(array, fc, fa, own, nei)"""
    return generator.getCellCenters(array, fc, fa, own, nei)

def getNormalMap(array):
    """Return the map of surface normals in an array.
    Usage: getNormalMap(array)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            b.append(generator.getNormalMap(i))
        return b
    else:
        return generator.getNormalMap(array)
  
# IN: niter: nbre d'iterations de lissage
# IN: eps: possible float ou numpy de taille des vertex: eps lissage
# IN: cellN: si present, extrapole en cas de cellN=0
# IN: algo: si 0: lisse n, si 1: lisse dn
def getSmoothNormalMap(array, niter=2, eps=0.4, cellN=None, algo=0):
    """Return the map of smoothed surface normals in an array.
    Usage: getSmoothNormalMap(array, niter, eps)"""
    it = 1
    n = getNormalMap(array)
    n = C.normalize(n, ['sx','sy','sz'])

    if cellN is not None:
        fake = ['cellN',cellN[1],n[2],n[3]]
        n = C.addVars([n, fake])
        generator.extrapWithCellN(array, n)
        n = C.extractVars(n, ['sx','sy','sz'])

    n = C.center2Node(n)
    n = C.normalize(n, ['sx','sy','sz'])
    
    while it < niter:
        np = C.node2Center(n)
        np = C.normalize(np, ['sx','sy','sz'])
        
        if cellN is not None: 
            fake = ['cellN',cellN[1],n[2],n[3]]
            np = C.addVars([np, fake])
            generator.extrapWithCellN(array, np)
            np = C.extractVars(np, ['sx','sy','sz'])

        np = C.center2Node(np)
        np = C.normalize(np, ['sx','sy','sz'])
        it += 1
        if algo == 0:
            if isinstance(eps, float): n[1][:] += eps*np[1][:]
            else: n[1][:] += eps[:]*np[1][:]
        else:
            if isinstance(eps, float): n[1][:] += eps*(np[1][:]-n[1][:])
            else: n[1][:] += eps[:]*(np[1][:]-n[1][:])
        n = C.normalize(n, ['sx','sy','sz'])
    return n

# identique a getSmoothNormalMap mais utilisant smoothField
def getSmoothNormalMap2(array, niter=2, eps=0.4, algo=0):
    try: import Transform as T
    except: raise ImportError("getSmoothNormalMap: requires Converter, Transform module.")
    n = getNormalMap(array)
    n = C.center2Node(n)
    T._smoothField(n, eps=eps, niter=niter, type=0, varNames=['sx','sy','sz'])
    n = C.normalize(n, vars=['sx','sy','sz'])
    return n

#=============================================================================
# Computes the ratio between the max and min length of all the edges of cells
#=============================================================================
def getEdgeRatio(array, dim=3):
    """Computes the ratio between the max and min lengths of all the edges of
    cells in an array.
    Usage: getEdgeRatio(a)"""
    if isinstance(array[0], list): 
        b = []
        for i in array: b.append(generator.getEdgeRatio(i, dim))
        return b
    else: return generator.getEdgeRatio(array, dim)

#=============================================================================
# Computes the max length of all the edges of cells
#=============================================================================
def getMaxLength(array, dim=3):
    """Computes the max length of all the edges of cells in an array.
    Usage: getMaxLength(a)"""
    if isinstance(array[0], list): 
        b = []
        for i in array: b.append(generator.getMaxLength(i, dim))
        return b
    else: return generator.getMaxLength(array, dim)

#=============================================================================
# Generate a list of collar grids depending on the assembly type
#=============================================================================
def collarMesh(s1, s2, distribj,distribk, niterj=100, niterk=100, ext=10,
               alphaRef=30., type='union',
               contour=[], constraints1=[], constraints2=[], toldist=1.e-6):
    """Generates a collar mesh starting from s1 and s2 surfaces, distributions along the surfaces
    and along the normal direction, with respect to the assembly type between grids.
    Usage: collarMesh(s1,s2,distribj,distribk,niterj,niterk,ext, alphaRef,type,contour,constraints1,constraints2,toldist)"""
    try: from . import Collar
    except: raise ImportError("collarMesh: requires Collar module.")
    if isinstance(s1[0], list): surfaces1 = s1
    else: surfaces1 = [s1]
    if isinstance(s2[0], list): surfaces2 = s2
    else: surfaces2 = [s2]

    infos = Collar.createCollarMesh__(surfaces1, surfaces2, distribj, distribk,\
                                      niterj, niterk, ext, alphaRef, type,\
                                      contour, constraints1, constraints2, toldist)
    A = []
    for info in infos: A+=[info[0]]
    return A

#=============================================================================
# generates a structured collar surface starting from a list of surfaces
# and the curve c shared by them
# les i-array contraintes doivent etre orientes dans le sens de la marche
# alphaRef is the deviation angle wrt 180 deg above which the walk is stopped
#=============================================================================
def surfaceWalk(surfaces, c, distrib, constraints=[], niter=0,
                alphaRef=180., check=0, toldist=1.e-6):
    """Generate a surface mesh by a walk on a list of surfaces, starting from
    a contour c and following constraints. niter is the number of smoothings.
    if check=1, walk stops before negative cells appear.
    Usage: surfaceWalk(surfaces, c, distrib, constraints, niter, alphaRef,
    check, toldist)"""
    from . import SurfaceWalk
    res = SurfaceWalk.surfaceWalk__(surfaces, c, distrib, constraints, niter,
                                    alphaRef, check, toldist)
    return res
    
#==============================================================================
# builds an extension starting from a contour c using the normals to surfaces
# and a point distribution dh
#==============================================================================
def buildExtension(c, surfaces, dh, niter=0):
    """Build an extension zone starting from contour c with respect to normals (smoothed
    niter times) to surfaces"""
    from . import SurfaceWalk as SW
    return SW.buildExtension__(c, surfaces, dh, niter)

#==============================================================================
# Determines the factor h0 to multiply the hloc in the vicinity of sharp
# angles. 
# IN : sn: normales a la surface s non structuree
# OUT: ht the amplification factor of h, the step in the marching direction  
#==============================================================================
def getLocalStepFactor__(s, sn, smoothType, nitLocal, kappaType, kappaS, algo):
    #import Transform
    sc = C.addVars([s,sn])
    sc = C.node2Center(sc)
    if algo == 0:
        ht = generator.getLocalStepFactor(s, sc)
        # Lissage hauteur du pas
        niter = 10; it = 0
        while it < niter:
            ht = C.node2Center(ht)
            ht = C.center2Node(ht)
            it += 1
        #ht = Transform.smoothField(ht, eps=0.25, niter=niter, varNames=['ht'])
            
    else:
        # Nouvelle version
        ht = generator.getLocalStepFactor2(s, sc, kappaType, kappaS[0], kappaS[1])
        # Lissage epsilon local
        it = 0
        hl = C.extractVars(ht, ['hl'])
        if smoothType == 2:
            hl2 = C.node2Center(hl)
            hl2 = C.center2Node(hl2)
            dif = abs(hl2[1]-hl[1])<1.

        while it < nitLocal:
            hl = C.node2Center(hl)
            hl = C.center2Node(hl)
            it += 1
        #hl = Transform.smoothField(hl, eps=0.25, niter=nitLocal, varNames=['hl'])
        
        if smoothType == 2: hl[1] = hl[1]*dif

        ht[1][1,:] = hl[1][0,:]
        if smoothType == 0: ht[1][1,:] = 1. # force constant epsilon
    return ht

#===============================================================================
# Regle la hauteur des normales, retourne aussi le champ pour le lisseur (algo=1)
def modifyNormalWithMetric(array, narray, algo=0, smoothType=0, eps=0.4, nitLocal=3, kappaType=0, kappaS=[0.2,1.6]):
    """Correct the normals located at nodes with respect to metric."""
    a = C.copy(array); n = C.copy(narray)
    if len(a) == 5: a = C.convertArray2Hexa(a); n = C.convertArray2Hexa(n)
    ht = getLocalStepFactor__(a, n, smoothType, nitLocal, kappaType, kappaS, algo)
    nx = C.extractVars(n,['sx']); ny = C.extractVars(n,['sy']); nz = C.extractVars(n,['sz'])
    nx[1][0,:] = ht[1][0,:]* nx[1][0,:]
    ny[1][0,:] = ht[1][0,:]* ny[1][0,:]
    nz[1][0,:] = ht[1][0,:]* nz[1][0,:]
    if algo == 1: epsl = C.extractVars(ht,['hl']); epsl[1][:] = eps*epsl[1][:]
    else: epsl = None
    n = C.addVars([n, nx, ny, nz])
    return n, epsl

def getCircumCircleMap(array):
    """Return the map of circum circle radius of any cell in a TRI array.
    Usage: getCircumCircleMap(array)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            b.append(generator.getCircumCircleMap(i))
        return b
    else:
        return generator.getCircumCircleMap(array)

def getInCircleMap(array):
    """Return the map of inscribed circle radius of any cell in a TRI array.
    Usage: getInCircleMap(array)"""
    if isinstance(array[0], list): 
        b = []
        for i in array: b.append(generator.getInCircleMap(i))
        return b
    else: return generator.getInCircleMap(array)

# -- Add normal layers to a surface --
def addNormalLayers(surface, distrib, check=0, 
                    niterType=0, niter=0, niterK=[], 
                    smoothType=0, eps=0.4, nitLocal=3, 
                    kappaType=0, kappaS=[0.2,1.6], blanking=False, cellNs=[],
                    algo=0):
    """Generate N layers to a surface following normals. distrib is the 
    height of the layers.
    If niter=0, the normal are not smoothed; else niter is the number of
    smoothing iterations applied to normals.
    Usage: addNormalLayers(surface, distrib, check, niter)"""
    if isinstance(surface[0], list): # liste d'arrays
        if len(surface[0]) == 5: type = 0 # structured
        else: type = 1 # unstructured
        for i in range(1,len(surface)):
            if (len(surface[i]) == 5 and type != 0) or (len(surface[i]) == 4 and type != 1):
                raise ValueError("addNormalLayers: all the surfaces must be structured or unstructured.")
        if type == 0: return addNormalLayersStruct__(surface, distrib, check, niterType, niter, niterK, smoothType, eps, nitLocal, kappaType, kappaS, blanking, cellNs, algo)
        else: # NS
            out = []
            for s in surface: out.append(addNormalLayersUnstr__(s, distrib, check, niterType, niter, niterK, smoothType, eps, nitLocal, kappaType, kappaS, blanking, cellNs, algo))
            return out
    else: # 1 array
        if len(surface) == 5: return addNormalLayersStruct__([surface], distrib, check, niterType, niter, niterK, smoothType, eps, nitLocal, kappaType, kappaS, blanking, cellNs, algo)[0]
        else: return addNormalLayersUnstr__(surface, distrib, check, niterType, niter, niterK, smoothType, eps, nitLocal, kappaType, kappaS, blanking, cellNs, algo) # NS

#-----------------------------------------------------------------------------
# Generation de grilles cartesiennes multibloc a partir de:
# grilles de corps: bodies
# h: pas d'espace sur la grille la plus fine
# Dfar: distance d eloignement au corps
# nlvl: nb de points par niveau, sauf sur le niveau le plus fin
# nlvl[0]: niveau le plus grossier
#-----------------------------------------------------------------------------
def gencartmb(bodies, h, Dfar, nlvl):
    """Generate a muliblock Cartesian mesh."""
    import KCore
    try: import Transform as T
    except:
        raise ImportError("gencartmb: requires Transform and Converter.")
    for b in bodies:
        ni = b[2]; nj = b[3]; nk = b[4]
        if ni < 2 or nj < 2 or nk < 2:
            raise ValueError("gencartmb: arrays must be 3D.")
        
    # Cree un bloc
    # IN: pmin: indices min sur la grille composite
    # IN: pmax: indices max sur la grille composite
    # IN: level: niveau de raffinement (1 = grossier)
    # IN: ref: grille composite
    def createBlock( pmin, pmax, level, ref):
        Href = ref[1][0,1] - ref[1][0,0]
        xmin = ref[1][0,0] + (pmin[0]-1)*Href
        ymin = ref[1][1,0] + (pmin[1]-1)*Href
        zmin = ref[1][2,0] + (pmin[2]-1)*Href
        #xmax = ref[1][0,0] + (pmax[0]-1)*Href
        #ymax = ref[1][1,0] + (pmax[1]-1)*Href
        #zmax = ref[1][2,0] + (pmax[2]-1)*Href
        hloc = 2**(level-1)
        Ni = (pmax[0]-pmin[0])*hloc+1
        Nj = (pmax[1]-pmin[1])*hloc+1
        Nk = (pmax[2]-pmin[2])*hloc+1
        return cart((xmin,ymin,zmin), (Href/hloc, Href/hloc, Href/hloc),
                    (Ni,Nj,Nk))
    # Cree un niveau
    # IN: nb: epaisseur en nombre de cellules de la grille composite
    # IN: level: level du niveau
    # IN/OUT: ref: grille composite
    def createLevel(nb, level, ref):
        out = []
        pmin = (1,1,1)
        pmax = (ref[2], ref[3], nb)
        if ref[4]-2*nb < 3  or  ref[3]-2*nb < 3 or ref[2]-2*nb < 3:
            print('Warning: number of points for level %d is too big: %d'%(level, nb))
            print('composite grid: %d %d %d.'%(ref[2],ref[3],ref[4])) 
            return out
        out.append(createBlock(pmin, pmax, level, ref))

        ref = T.subzone(ref, (1,1,nb), (ref[2],ref[3],ref[4]))
        
        pmin = (1,1,ref[4]-nb+1)
        pmax = (ref[2], ref[3], ref[4])
        out.append(createBlock(pmin, pmax, level, ref))
        ref = T.subzone(ref, (1,1,1), (ref[2], ref[3], ref[4]-nb+1))
        
        pmin = (1,1,1)
        pmax = (nb, ref[3], ref[4])
        out.append(createBlock(pmin,pmax,level,ref))
        ref = T.subzone(ref, (nb,1,1), (ref[2], ref[3], ref[4]))
        
        pmin = (ref[2]-nb+1,1,1)
        pmax = (ref[2], ref[3], ref[4])
        out.append(createBlock(pmin,pmax,level,ref))
        ref = T.subzone(ref, (1,1,1), (ref[2]-nb+1, ref[3], ref[4]))
        
        pmin = (1,1,1)
        pmax = (ref[2], nb, ref[4])
        out.append(createBlock(pmin,pmax,level,ref))
        ref = T.subzone(ref, (1,nb,1), (ref[2], ref[3], ref[4]))
        
        pmin = (1,ref[3]-nb+1,1)
        pmax = (ref[2], ref[3], ref[4])
        out.append(createBlock(pmin,pmax,level,ref))
        ref = T.subzone(ref, (1,1,1), (ref[2], ref[3]-nb+1, ref[4]))
        return (out, ref)

    # Bounding box
    bb = bbox(bodies)
    xmin = bb[0] - Dfar
    ymin = bb[1] - Dfar
    zmin = bb[2] - Dfar
    xmax = bb[3] + Dfar
    ymax = bb[4] + Dfar
    zmax = bb[5] + Dfar
    
    # Grille composite
    nc = 2**len(nlvl)
    Href = h * nc
    Ni = (xmax - xmin)/Href; Ni = nc*(int(Ni/nc)+1)+1
    Nj = (ymax - ymin)/Href; Nj = nc*(int(Nj/nc)+1)+1
    Nk = (zmax - zmin)/Href; Nk = nc*(int(Nk/nc)+1)+1
    ref = cart((xmin,ymin,zmin), (Href,Href,Href), (Ni,Nj,Nk))
    
    # Niveaux
    out = []
    c = 1
    for i in nlvl:
        lev,ref = createLevel(i, c, ref)
        c += 1
        out = out + lev

    # Derniere grille
    pmin = (1,1,1)
    pmax = (ref[2],ref[3],ref[4])
    level = 4
    out.append(createBlock(pmin, pmax, level, ref))

    tx = xmax-xmin-Href*(Ni-1)
    ty = ymax-ymin-Href*(Nj-1)
    tz = zmax-zmin-Href*(Nk-1)
    out = T.translate(out, (0.5*tx,0.5*ty,0.5*tz))
    return out 

#------------------------------------------------------------------------------
# Split a i-array and map a distribution on the splitted i-array
# splitCrit is the sensibility
#------------------------------------------------------------------------------
def mapSplit(array, dist, splitCrit=100., densMax=1000):
    """Split a curve and map a distribution on the set of split curves.
    Usage: mapSplit(a, dist, splitCrit, densMax)"""
    if len(array) == 5: # structure
        if array[3] != 1 or array[4] != 1:
            raise TypeError("mapSplit: requires a i-array.")
        else:
            return mapSplitStruct__(array, dist, splitCrit, densMax)
    elif len(array) == 4: # non structure
        raise TypeError("mapSplit: requires a i-array.") 
    else:
        raise TypeError("mapSplit: requires a i-array.")        

# mapSplit pour un i-array structure
def mapSplitStruct__(array, dist, splitCrit, densMax):
    import KCore
    try: 
        import math
        import Transform as T
        import Geom as D
    except:
        raise ImportError("mapSplit: requires Converter, Geom, Transform modules.")

    a = C.copy(array); d = C.copy(dist)

    # Compute total length of dist
    ldt = D.getLength(d); ldti = 1./ldt

    posx = KCore.isNamePresent(d, "x")
    posy = KCore.isNamePresent(d, "y")
    posz = KCore.isNamePresent(d, "z")
    if posx == -1: posx = KCore.isNamePresent(a, "CoordinateX")
    if posy == -1: posy = KCore.isNamePresent(a, "CoordinateY")
    if posz == -1: posz = KCore.isNamePresent(a, "CoordinateZ")
    if posx == -1 or posy == -1 or posz == -1:
        raise TypeError("mapSplit: coordinates not found in an distribution.")
    x = d[1][posx]; y = d[1][posy]; z = d[1][posz]
    nbpoints = len(x)
    ld = numpy.zeros(nbpoints, dtype='float64'); ld[0] = 0.

    for i in range(1, nbpoints):
        dx = x[i]-x[i-1]; dy = y[i]-y[i-1]; dz = z[i]-z[i-1]
        ld[i] = ld[i-1] + math.sqrt(dx*dx+dy*dy+dz*dz)*ldti

    # Compute total length of array
    lt = D.getLength(a); lti = 1./lt
    # Compute minimum discretization step so as to densify array
    stepmin = lt
    posx = KCore.isNamePresent(d, "x")
    posy = KCore.isNamePresent(d, "y")
    posz = KCore.isNamePresent(d, "z")
    if posx == -1: posx = KCore.isNamePresent(a, "CoordinateX")
    if posy == -1: posy = KCore.isNamePresent(a, "CoordinateY")
    if posz == -1: posz = KCore.isNamePresent(a, "CoordinateZ")
    if posx == -1 or posy == -1 or posz == -1:
        raise TypeError("mapSplit: coordinates not found in array.")

    xa = a[1][posx]; ya = a[1][posy]; za = a[1][posz]

    for i in range(1,len(xa)):
        dxa = xa[i]-xa[i-1]; dya = ya[i]-ya[i-1]; dza = za[i]-za[i-1]
        la = math.sqrt(dxa*dxa+dya*dya+dza*dza) *lti 
        if la < stepmin and la > 0.: stepmin = la
    if lt/stepmin > densMax: stepmin = lt/densMax
    # Densify array
    a = densify(a, stepmin)
    # Split array
    a = T.splitCurvatureRadius(a, splitCrit)

    # Build array with distance of each "split point"
    L = numpy.zeros((len(a)), dtype='float64')
    ltinv = 1./lt
    L[0] = D.getLength(a[0]) * ltinv
   
    for i in range(1,len(a)):
        L[i] = L[i-1] + D.getLength(a[i]) * ltinv
     
    # Find indices in dist which correspond to "split points"
    indsplit_previous = 0
    for i in range(len(L)):
        ind2 = D.getDistantIndex(dist, 1, L[i]*ldt)
        ind1 = ind2-1

        d1 = math.fabs(L[i]-ld[ind1])
        d2 = math.fabs(L[i]-ld[ind2])
        if d1 < d2: indsplit = ind1
        else: indsplit = ind2
        # Split distribution with "split points"
        if indsplit_previous < indsplit:
            ds = T.subzone(d, (indsplit_previous+1,1,1), (indsplit+1,1,1))
            hom = 1./((ld[indsplit]-ld[indsplit_previous])*ldt)
            xs = ds[1][posx]; ys = ds[1][posy]; zs = ds[1][posz]
            newdist = T.homothety(ds, (xs[0],ys[0],zs[0]),hom)
            newdist = T.translate(newdist, (-xs[0],-ys[0],-zs[0]))
            a[i] = map(a[i], newdist)
            indsplit_previous = indsplit
    return a

def T3mesher2D(a, grading=1.2, triangulateOnly=0, metricInterpType=0):
    """Create a delaunay mesh given a set of points defined by a.
    Usage: T3mesher2D(a, grading, triangulateOnly, metricInterpType)"""
    try:
        b = C.convertArray2Tetra(a); b = close(b)
        return generator.T3mesher2D(b, grading, triangulateOnly, metricInterpType)
    except:
        return generator.T3mesher2D(a, grading, triangulateOnly, metricInterpType)

def tetraMesher(a, maxh=-1., grading=0.4, triangulateOnly=0, 
                remeshBoundaries=0, algo=1, optionString=""):
    """Create a TRI/TETRA mesh given a set of BAR or surfaces in a.
    Usage: tetraMesher(a, maxh, grading)"""
    try:
        import Transform as T
        a = C.convertArray2Tetra(a)
        a = T.join(a); a = close(a)
    except: pass
    import math
    if a[3] == 'BAR':
        p = fittingPlaster(a)
        b = gapfixer(a, p)
        return b
        #return generator.T3mesher2D(a, triangulateOnly)
    elif a[3] == 'TRI':
        if maxh == -1.: # auto maxh
            vol = getVolumeMap(a)
            maxh = C.getMeanValue(vol, 'vol')
            maxh = math.sqrt(2*maxh)
        if algo == 0: # netgen
            if remeshBoundaries == 0:
                if maxh < 0: maxh = 1.e6
                return generator.netgen1(a, maxh, grading)
            else:
                C.convertArrays2File([a], 'netgen.0120.stl', 'fmt_stl')
                if maxh < 0: maxh = 1.e6
                return generator.netgen2(a, maxh, grading)
        else: # tetgen
            # find holes coords
            holes = []
            try: 
                import Transform as T
                sp = T.splitConnexity(a)
                for s in sp:
                    px = C.isNamePresent(s, 'x')
                    if px != -1: ext = C.extractVars(s, ['x','y','z'])
                    else: ext = C.extractVars(s, ['CoordinateX','CoordinateY','CoordinateZ'])
                    n = getNormalMap(s)
                    n = C.center2Node(n)
                    n = C.normalize(n, ['sx','sy','sz'])
                    pt = C.getValue(ext, 0)
                    n = C.getValue(n, 0)
                    eps = 1.e-10
                    holes.append([pt[0]+eps*n[0], pt[1]+eps*n[1], pt[2]+eps*n[2]])
            except: pass
            return generator.tetgen(a, maxh, grading, remeshBoundaries, holes, optionString)
    else:
        raise TypeError("tetraMesher: requires BAR or TRI mesh.")

def fittingPlaster(contour, bumpFactor=0.):
    """Generate a structured points cloud over a BAR.
    Usage: fittingPlaster(contour, bumpFactor)"""
    try: contour = C.convertArray2Tetra(contour)
    except: pass
    contour = close(contour)
    return generator.fittingPlaster(contour, bumpFactor)

def gapfixer(contour, cloud, hardPoints=None, refine=1):
    """Fix a gap defined by a contour bar and a point cloud representing the gap surface.
    Some hard points can be specified to force the constructed surface to pass by.
    If the optional refine argument is set to 0, the resulting surface will be a constrained triangulation of the contour [and the additional hard nodes].
    Usage: gapFixer(contour, cloud, hardPoints, refine)"""
    try: contour = C.convertArray2Tetra(contour)
    except: pass
    contour = close(contour)
    return generator.gapfixer(contour, cloud, hardPoints, refine)

def gapsmanager(components, mode=0, refine=0, coplanar=0):
    """Fix a gap between several component surfaces (list of arrays).
    The mode sets the case (0 for POST with a center mesh, 1 for POST with a nodal mesh, 2 otherwise).
    The planar argument tells whether the components are coplanar or not (1 for coplanar case).
    Usage: gapsmanager(components, mode, refine, coplanar)"""
    try:
        components = C.convertArray2Tetra(components)
        components = close(components)
    except: pass
    return generator.gapsmanager(components, mode, refine, coplanar)

def front2Hexa(a, surf, h, hf, hext, density=50):
    """Generate an hexa grid starting from a front a, a surface surf,
    and h, hf, hext the height of the mesh, of the first and last cells."""
    try: import Transform as T; import math
    except: raise ImportError("front2Hexa: requires Transform module.")

    # projection de a sur la surface ortho
    b = T.projectOrtho(a, [surf]) # projection sur la surface: quad surf
    
    # Distribution dans la direction k
    npts = a[1][0].shape[0]; h0 = 0.
    a1 = a[1]; b1 = b[1]
    for ind in range(npts):
        dx=a1[0,ind]-b1[0,ind]; dy=a1[1,ind]-b1[1,ind]; dz=a1[2,ind]-b1[2,ind]
        h0 = max(h0, dx*dx+dy*dy+dz*dz)
    h0 = math.sqrt(h0)
    h0 = max(h, h0) # hauteur max du maillage (sans extension)
    nk = int(h0*density)+1; nk = max(nk,5)
    distrib = cart((0.,0.,0.),(1./nk,1.,1.),(nk+1,1,1))
    supp = nk-1; add = 10
    distrib = enforceMoinsX(distrib, hf/h, supp, add)
    distrib = enforcePlusX(distrib, hext/h, supp, add)
    # for debug
    distrib = cart((0.,0.,0.),(1./1,1.,1.),(2,1,1))
    return generator.front2Hexa(b, a, distrib, h0)

def front2Struct(front, surf, distrib, Vmin, dist):
    """Generate struct grids starting from a front, a surface surf,
    and a point distribution."""  
    return generator.front2Struct(front, surf, distrib, Vmin, dist)

def snapFront(meshes, surfaces, optimized=1):
    """Adapt meshes to a given surface (cellN defined). 
    Usage: snapFront(meshes, surfaces)"""
    try: import Transform as T
    except:
        raise ImportError("snapFront: requires Converter, Transform module.")

    try:
        import Transform as T
        b = C.convertArray2Tetra(surfaces); b = T.join(b)
    except: b = surfaces[0]

    if isinstance(meshes[0], list):
        return generator.snapFront(meshes, b, optimized)
    else:
        return generator.snapFront([meshes], b, optimized)[0]

def snapSharpEdges(meshes, surfaces, step=None, angle=30.):
    """Adapt meshes to a given surface sharp edges. 
    Usage: snapSharpEdges(meshes, surfaces)"""
    out = refinedSharpEdges__(surfaces, step, angle)
    contours = out[1]
    if contours == []: contours = None
    corners = out[2]
    if corners == []: corners = None

    #if contours is not None:
    #    C.convertArrays2File(contours, 'contours.plt')
    #if corners is not None:
    #    C.convertArrays2File(corners, 'corners.plt')

    if isinstance(meshes[0], list):
        if contours is not None:
            meshes = generator.snapSharpEdges(meshes, contours)
        if corners  is not None:
            meshes = generator.snapSharpEdges(meshes, corners)
        return meshes
    else:
        if contours is not None:
            meshes = generator.snapSharpEdges([meshes], contours)[0]
        if corners is not None:
            meshes = generator.snapSharpEdges([meshes], corners)[0]
        return meshes

#------------------------------------------------------------------------------
# Get refined sharp edges from a given surface
# IN: surfaces ou contours
# IN: step: step for refinement
# IN: angle: angle for surfaces splitting
# OUT: outList: list which contains:
#               - the surfaces joined in a unique array
#               - the contours (if not nul)
#               - the corners (if not nul)
#------------------------------------------------------------------------------
def refinedSharpEdges__(surfaces, step, angle):
    """Get refined sharp edges from a given surface. 
    Usage: snapSharpEdges(meshes, surfaces)"""
    try:
        import Post as P; import Geom as D 
        import Transform as T; from . import Generator as G
    except:
        raise ImportError("snapSharpEdges: requires Post, Geom, Converter, Transform module.")
    b = C.convertArray2Tetra(surfaces); b = T.join(b); b = close(b)

    # dimension de surfaces: 1D ou 2D
    dim = 2
    if b[3] == 'BAR': dim = 1
    
    if dim == 2:
        # get contours and corners from 2D-surfaces
        try: contours = P.sharpEdges(b, angle)
        except: contours = []; corners = []
           
        if contours != []:
            contours = T.splitConnexity(contours)
            # split les contours par rapport aux angles
            try: contours = T.splitSharpEdges(contours, angle)
            except: pass
            # split les contours non structures par rapport aux branches
            try: contours = T.splitTBranches(contours)
            except: pass
            # get corners from contours
            try: corners = P.sharpEdges(contours, angle)
            except: corners = []
    else: # 1D
        # get contours and corners from 1D-surfaces
        try: contours = T.splitConnexity(b)
        except: contours = b 
        # split les contours par rapport aux angles
        try: contours = T.splitSharpEdges(contours, angle)
        except: pass
        # split les contours non structures par rapport aux branches
        try: contours = T.splitTBranches(contours)
        except: pass
        try: corners = P.sharpEdges(b, angle)
        except: corners = []
        #contours = [] # force pour debug

    # remaillage des contraintes (contours) en fonction
    # d'un pas defini par l'utilisateur
    ncontours = []
    if step is not None:
        for c in contours:
            # conversion en structure pour map
            c = C.convertBAR2Struct(c)
            l = D.getLength(c)
            N = int(l/step)+1
            N = max(N, 3)
            distrib = cart((0,0,0), (1./(N-1),1,1), (N,1,1))
            c = map(c, distrib)
            ncontours.append(c)
    else: ncontours = contours

    # les contraintes contours et corners sont passes en un array unique
    # non-structure
    if ncontours != []:
        contours =  C.convertArray2Tetra(ncontours)
        contours = T.join(contours)
        contours = G.close(contours)
    if corners != []: corners = T.join(corners)
    return [b, contours, corners]
    
def fillWithStruct(a, Vmin):
    """Generates struct grids in quad mesh."""  
    return generator.fillWithStruct(a, Vmin)

def octree2Struct(a, vmin=15, ext=0, optimized=1, merged=1, AMR=0,
                  sizeMax=1000000000):
    """Generates a structured set of regular Cartesian grid starting from
    an octree HEXA or QUAD mesh. vmin is the number of minimum points per grid,
    and can be specified for each level;
    ext is the extension of grids in all the directions.
    If optimized=1, then the extension can be reduced for minimum overlapping.
    merged=1 means that Cartesian grids are merged when possible.
    If AMR=1, a list of AMR grids is generated.
    Usage: octree2Struct(a, vmin, ext, optimized, merged, AMR)"""
    if not isinstance(vmin, list): vmin = [vmin]
    for nov in range(len(vmin)):
        if vmin[nov] < 2:
            print('Warning: octree2Struct, vmin is set to 2.'); vmin[nov] = 2
        if ext == 0 and vmin[nov]%2 == 0:
            vmin[nov] += 1
            print('Warning: octree2Struct: vmin must be odd, vmin set to %d.'%vmin[nov])
    if AMR == 0: cartzones = generator.octree2Struct(a, vmin)
    else: cartzones = generator.octree2AMR(a, vmin[0])
    if merged == 1:
        try:
            import Transform
            #cartzones = Transform.mergeCartByRefinementLevel(cartzones,sizeMax)
            cartzones = Transform.mergeCart(cartzones, sizeMax)
        except: pass
    if optimized != 1 and optimized != 0:
        print('Warning: octree2Struct: optimized must be 0 or 1. Set to 1.')
        optimized = 1

    if ext == 0: return cartzones
    elif ext > 0: return extendOctreeGrids__(cartzones,
                                             ext=ext, optimized=optimized, extBnd=0)
    else:
        print('Warning: octree2Struct: ext must be equal or greater than 0. Set to 0.')
    return cartzones

# Decoupe un octant (octant)
# N multiple de 4 (2D) ou de 8 (3D) 
def cutOctant(octant, N, ind, dim=3):
    if dim == 2:
        M = N**0.5
        dx = octant[2]-octant[0]
        dy = octant[3]-octant[1]
        dx = dx/M; dy = dy/M
        j = int(ind/M)
        i = ind - j*M
        x0 = octant[0]+dx*i
        y0 = octant[1]+dy*j
        return [x0,y0,x0+dx,y0+dy]
    if dim == 3:
        M = N**(1./3.)
        dx = octant[2]-octant[0]
        dy = octant[3]-octant[1]
        dz = octant[4]-octant[2]
        dx = dx/M; dy = dy/M; dw = dz/M
        k = int(ind/(M*M))
        j = int((ind-k*M*M)/M)
        i = ind - j*M - k*M*M
        x0 = octant[0]+dx*i
        y0 = octant[1]+dy*j
        z0 = octant[2]+dz*k
        return [x0,y0,z0,x0+dx,y0+dy,z0+dz]

def octree(stlArrays, snearList=[], dfarList=[], dfar=-1., balancing=0, levelMax=1000, ratio=2, octant=None, dfarDir=0, mode=0):
    """Generate an octree (or a quadtree) mesh starting from a list of TRI (or BAR) arrays defining bodies,
    a list of corresponding snears, and the extension dfar of the mesh."""
    try: s = C.convertArray2Tetra(stlArrays)
    except: s = stlArrays
    if ratio == 2:
        o = generator.octree(s, snearList, dfarList, dfar, levelMax, octant, dfarDir, mode)
        if balancing == 0: return o
        elif balancing == 1: return balanceOctree__(o, 2, corners=0)
        elif balancing == 2: return balanceOctree__(o, 2, corners=1)
    else: #27-tree
        o = generator.octree3(s, snearList, dfar, levelMax, octant)
        if balancing == 0: return o
        else: return balanceOctree__(o, 3, corners=0)

def conformOctree3(octree):
    """Conformize an octree3.
    Usage: conformOctree3(octree)"""
    return generator.conformOctree3(octree)
    
def balanceOctree__(octree, ratio=2, corners=0):
    return generator.balanceOctree(octree, ratio, corners)

def extendCartGrids(A, ext=0, optimized=0, extBnd=0):
    A, rinds = generator.extendCartGrids(A, ext, optimized, extBnd)
    return A, rinds

def extendOctreeGrids__(A, ext, optimized, extBnd=0):
    """Extend grids with ext cells. If optimized is ext, the minimum overlapping is ensured.
    Usage: extendOctreeGrids__(cartGrids, ext, optimized, extBnd)"""
    A, rinds = extendCartGrids(A, ext, optimized, extBnd)
    return A

def adaptOctree(octreeHexa, indicField, balancing=1, ratio=2):
    """Adapt an unstructured octree w.r.t. the indicatorField -n,0,n defined for each element.
    Usage: adaptOctree(o, indicField, balancing)"""
    ok = 1
    hexa = octreeHexa; indic = indicField
    while ok != 0:
        if ratio == 2: res = generator.adaptOctree(hexa, indic)
        else: res = generator.adaptOctree3(hexa, indic)
        hexa = res[0]; indic = res[1]
        ok = 0
        indic2 = indic[1][0]
        ret = numpy.count_nonzero(indic2)
        if ret > 0: ok = 1
    if balancing == 0: return hexa
    elif balancing==1: return balanceOctree__(hexa, ratio, corners=0)
    elif balancing==2: return balanceOctree__(hexa, ratio, corners=1)
    else: raise ValueError("adaptOctree: bad value for balancing argument.")

def expandLayer(octreeHexa, level=0, corners=0, balancing=0):
    """Expand the layer of octree elements of level l of one additional layer.
    The initial octree must be balanced.
    Usage: expandLayer(octree, level, corners, balancing)"""
    typeExpand=2 # check neigbours only
    indic = C.node2Center(octreeHexa)
    indic = C.initVars(indic, 'indicator', 0.)
    indic = generator.modifyIndicToExpandLayer(octreeHexa, indic,
                                               level, corners, typeExpand)
    return adaptOctree(octreeHexa, indic, balancing)

def forceMatch(a1, a2, tol=1.):
    b1 = C.copy(a1)
    _forceMatch(b1, a2, tol)
    return b1
    
# find the best of two numpys
def findBest(diff, bary):
    ind1 = numpy.argmin(diff)
    ind2 = numpy.argmin(bary)
    if ind1 == ind2: return ind1
    diff2 = numpy.copy(diff)
    diff[ind1] = 1.e6
    ind1 = numpy.argmin(diff)
    if ind1 == ind2: return ind1
    diff = diff2
    bary[ind2] = 1.e6
    ind2 = numpy.argmin(bary)
    if ind1 == ind2: return ind1
    return ind1
    
# force near boundary (<tol) of a1 to match with a2
# in place on a1 (TRI surfaces)
# if P0 and P1: point of ext of a1, used to find the
# piece of ext that must match

# Match l'exterieur de a1 sur l'exterieur de a2 si la distance est
# inferieure a tol
def _forceMatch1(a1, a2, tol):
    import Post; import KCore; import Geom; import Transform; import Generator
        
    # exterior of a1
    ext1 = Post.exteriorFaces(a1)
    
    # exterior of a2
    ext2 = Post.exteriorFaces(a2)

    # Get pos
    posx1 = KCore.isCoordinateXPresent(a1)
    posy1 = KCore.isCoordinateYPresent(a1)
    posz1 = KCore.isCoordinateZPresent(a1)
    posx2 = KCore.isCoordinateXPresent(a2)
    posy2 = KCore.isCoordinateYPresent(a2)
    posz2 = KCore.isCoordinateZPresent(a2)
    
    vol1 = getVolumeMap(ext1)
    vol1 = C.center2Node(vol1)[1]
    vol2 = getVolumeMap(ext2)
    vol2 = C.center2Node(vol2)[1]
    
    # identifie ext1 sur a1
    hook = C.createHook(a1, function='nodes')
    indices1 = C.identifyNodes(hook, ext1)
    
    # identifie ext2 sur a2
    hook = C.createHook(a2, function='nodes')        
    indices2 = C.identifyNodes(hook, ext2)
    
    # match ext1 sur ext2
    hook = C.createHook(ext2, function='nodes')
    nodes,dist = C.nearestNodes(hook, ext1)
    npts = C.getNPts(ext1)
    for i in range(npts):
        if dist[i] < tol*0.55*vol1[0,i]:
            ind1 = indices1[i]-1
            ind2 = nodes[i]-1
            ext1[1][posx1,i] = ext2[1][posx2,ind2]
            ext1[1][posy1,i] = ext2[1][posy2,ind2]
            ext1[1][posz1,i] = ext2[1][posz2,ind2]
            a1[1][posx1,ind1] = ext2[1][posx2,ind2]
            a1[1][posy1,ind1] = ext2[1][posy2,ind2]
            a1[1][posz1,ind1] = ext2[1][posz2,ind2]        
        
    # match ext2 sur ext1
    hook = C.createHook(ext1, function='nodes')
    nodes,dist = C.nearestNodes(hook, ext2)
    npts = C.getNPts(ext2)
    for i in range(npts):
        if dist[i] < tol*0.55*vol2[0,i]:
            ind1 = nodes[i]-1
            ind2 = indices2[i]-1
            a2[1][posx2,ind2] = ext1[1][posx1,ind1]
            a2[1][posy2,ind2] = ext1[1][posy1,ind1]
            a2[1][posz2,ind2] = ext1[1][posz1,ind1]        
    return None

# Force match sur la bande delimitee par P1-P2
def _forceMatch2(a1, a2, P1, P2):
    import Post; import KCore; import Geom; import Transform; import Generator
    
    # exterior of a1
    ext1 = Post.exteriorFaces(a1)
    
    # exterior of a2
    ext2 = Post.exteriorFaces(a2)
    
    # Find split index of P1 and P2
    hook = C.createHook(ext1, function='nodes')
    nodes,dist = C.nearestNodes(hook, Geom.point(P1))
    ind1s1 = nodes[0]-1
    nodes,dist = C.nearestNodes(hook, Geom.point(P2))    
    ind2s1 = nodes[0]-1
    hook = C.createHook(ext2, function='nodes')
    nodes,dist = C.nearestNodes(hook, Geom.point(P1))    
    ind1s2 = nodes[0]-1
    nodes,dist = C.nearestNodes(hook, Geom.point(P2))    
    ind2s2 = nodes[0]-1
    ext1 = Transform.splitBAR(ext1, ind1s1, ind2s1)
    ext2 = Transform.splitBAR(ext2, ind1s2, ind2s2)
    C.convertArrays2File(ext1+ext2, 'exts.plt')
    # on garde ceux qui ont a peu pres la meme longeur
    # et un barycentre semblable
    n1 = len(ext1); n2 = len(ext2)
    diff = numpy.empty(n1+n2, dtype=numpy.float64)
    bary = numpy.empty(n1+n2, dtype=numpy.float64)
    for i, e1 in enumerate(ext1):
        di = Geom.getLength(e1)
        bi = Generator.barycenter(e1)
        for j, e2 in enumerate(ext2):
            dj = Geom.getLength(e2)
            bj = Generator.barycenter(e2)
            diff[i+n1*j] = abs(di-dj)
            bary[i+n1*j] = (bi[0]-bj[0])**2+(bi[1]-bj[1])**2+(bi[2]-bj[2])**2            
    
    ind = findBest(diff, bary)
    j = ind//n1; i = ind-j*n1
    ext1 = ext1[i]; ext2 = ext2[j]
    C.convertArrays2File([ext1,ext2], 'exts.plt')
    _forceMatch3(a1, a2, ext1, ext2)    
    return None

# force match avec deux courbes en entree
def _forceMatch3(a1, a2, ext1, ext2):
    import Post; import KCore; import Geom; import Transform; import Generator
                
    # Get pos
    posx1 = KCore.isCoordinateXPresent(a1)
    posy1 = KCore.isCoordinateYPresent(a1)
    posz1 = KCore.isCoordinateZPresent(a1)
    posx2 = KCore.isCoordinateXPresent(a2)
    posy2 = KCore.isCoordinateYPresent(a2)
    posz2 = KCore.isCoordinateZPresent(a2)
    
    # identifie ext1 sur a1
    hook = C.createHook(a1, function='nodes')
    indices1 = C.identifyNodes(hook, ext1)
    
    # identifie ext2 sur a2
    hook = C.createHook(a2, function='nodes')        
    indices2 = C.identifyNodes(hook, ext2)
    
    # match ext1 sur ext2
    hook = C.createHook(ext2, function='nodes')
    nodes,dist = C.nearestNodes(hook, ext1)
        
    ext1[1][posx1,:] = ext2[1][posx2,nodes[:]-1]
    ext1[1][posy1,:] = ext2[1][posy2,nodes[:]-1]
    ext1[1][posz1,:] = ext2[1][posz2,nodes[:]-1]
    a1[1][posx1,indices1[:]-1] = ext2[1][posx2,nodes[:]-1]
    a1[1][posy1,indices1[:]-1] = ext2[1][posy2,nodes[:]-1]
    a1[1][posz1,indices1[:]-1] = ext2[1][posz2,nodes[:]-1]
    
    # match ext2 sur new ext1
    hook = C.createHook(ext1, function='nodes')
    nodes,dist = C.nearestNodes(hook, ext2)
    a2[1][posx2,indices2[:]-1] = ext1[1][posx1,nodes[:]-1]
    a2[1][posy2,indices2[:]-1] = ext1[1][posy1,nodes[:]-1]
    a2[1][posz2,indices2[:]-1] = ext1[1][posz1,nodes[:]-1]
        
    return None

# Pour contour interne a a
def _forceMatch4(a, ext1, ext2):
    _forceMatch3(a, a, ext1, ext2)
    return None

def _forceMatch(a1, a2=None, P1=None, P2=None, C1=None, C2=None, tol=-1):
    if a2 is not None and P1 is not None and P2 is not None:
        _forceMatch2(a1, a2, P1, P2)
    elif a2 is not None and C1 is not None and C2 is not None:
        _forceMatch3(a1, a2, C1, C2)
    elif a2 is None and C1 is not None and C2 is not None:
        _forceMatch4(a1, C1, C2)
    elif a2 is not None: 
        _forceMatch1(a1, a2, tol)
    return None
    
# addnormalLayers pour une liste d'arrays structures
def addNormalLayersStruct__(surfaces, distrib, check=0, niterType=0, niter=0, niterK=[], 
                            smoothType=0, eps=0.4, nitLocal=3, 
                            kappaType=0, kappaS=[0.2,1.6], blanking=False, cellNs=[],
                            algo=0):
    import KCore
    try: import Transform as T; import Generator as G
    except: raise ImportError("addNormalLayers: requires Converter, Transform modules.")   
    kmax = distrib[1].shape[1] # nb of layers in the normal direction

    if kmax < 2: raise ValueError("addNormalLayers: distribution must contain at least 2 points.")

    vect = ['sx','sy','sz']
    # verifications 
    for nos, surfs in enumerate(surfaces):
        if surfs[4] != 1: raise ValueError("addNormalLayers: structured surface must be k=1.")
        if surfs[3] == 1: surfaces[nos] = T.addkplane(surfs)

    surfu = C.convertArray2Hexa(surfaces)
    surfu = T.join(surfu); surfu = close(surfu)
    surfu = T.reorder(surfu, (1,))
    listOfIndices = KCore.indiceStruct2Unstr2(surfaces, surfu, 1.e-14)
    
    listOfCoords = []
    for surfs in surfaces:
        npts = surfs[1].shape[1]
        nfld = surfs[1].shape[0]
        coords = C.array(surfs[0], surfs[2], surfs[3], kmax)
        coords[1][:,0:npts] = surfs[1][:,0:npts]
        listOfCoords.append(coords)

    nzones = len(listOfCoords)
    hmin = distrib[1][0,1]-distrib[1][0,0]
    hmax = distrib[1][0,kmax-1]-distrib[1][0,kmax-2]
    hmean = (distrib[1][0,kmax-1]-distrib[1][0,0])/(kmax-1)

    # determination de kb1,kb2
    kb1 = -1; kb2 = -1
    for k1 in range(kmax-1):
        if distrib[1][0,k1+1] >= 0.1*hmean and kb1 == -1: kb1 = k1
        elif distrib[1][0,k1+1] >= 1.*hmean and kb2 == -1: kb2 = k1
    kb2 = max(kb2, kb1+2)
    imax = surfu[1].shape[1]
    coordsloc = C.array(surfu[0],imax,2,1) # pour le check
    coordsloc[1][0,0:imax] = surfu[1][0,0:imax]
    coordsloc[1][1,0:imax] = surfu[1][1,0:imax]
    coordsloc[1][2,0:imax] = surfu[1][2,0:imax]
    k1 = 0
    stop = 0
    epsl = None # champ pour le lissage defini sur surfu
    cellNp = [None]*nzones # champ de blanking
    for i in range(nzones): cellNs.append(None) # Champ complet de blanking
    cellNu = None # champ de cellN defini sur surfu
    if blanking:
        cellNu = C.array('cellN', surfu[1].shape[1], 1, 1)
        cellNu = C.initVars(cellNu, 'cellN', 1.)

    while k1 < kmax-1 and stop == 0:
        hloc = distrib[1][0,k1+1]-distrib[1][0,k1]
        if algo == 0: # algo=0, lissage partout, hauteur 1
            if niter == 0:
                n = getNormalMap(surfu)
                n = C.normalize(n, vect)
                n = C.center2Node(n)
                n = C.normalize(n, vect)
                n, epsl = modifyNormalWithMetric(surfu, n, algo=0, eps=eps)
            else: # lissage
                if hmin < 0.01*hmean: # presence de couche limite
                    if k1 < kb1: # pas de lissage ds la couche limite
                        n = getSmoothNormalMap(surfu, niter=0, eps=eps)
                        np, epsl = modifyNormalWithMetric(surfu, n, algo=0, eps=eps)
                        n[1] = np[1]
                    elif k1 < kb2:
                        beta0 = (float(k1-kb1))/float(kb2-1-kb1)
                        n0 = getSmoothNormalMap(surfu, niter=0, eps=eps)
                        n0, epsl = modifyNormalWithMetric(surfu, n0, algo=0, eps=eps)
                        n = getSmoothNormalMap(surfu, niter=niter, eps=eps)
                        np, epsl = modifyNormalWithMetric(surfu, n, algo=0, eps=eps)
                        n[1] = (1-beta0)*n0[1] + beta0*np[1]
                    else: # lissage a fond
                        n = getSmoothNormalMap(surfu, niter=niter, eps=eps)
                        np, epsl = modifyNormalWithMetric(surfu, n, algo=0, eps=eps)
                        beta0 = float((kmax-2-k1))/float(kmax-2)
                        beta0 = beta0*beta0
                        n[1] = (1-beta0)*n[1] + beta0*np[1]
                else: # pas de couche limite
                    n = getSmoothNormalMap(surfu, niter=niter, eps=eps)
                    np, epsl = modifyNormalWithMetric(surfu, n, algo=0, eps=eps)
                    if kmax == 2: beta0 = 0.1
                    else: beta0 = float((kmax-2-k1))/float(kmax-2); beta0 = beta0*beta0
                    n[1] = (1-beta0)*n[1] + beta0*np[1]
        else: # algo=1, lissage ponctuel, hauteur 2
            if k1 == 0 and kappaType == 2:
                vol0 = getVolumeMap(surfu)
                vol0 = C.center2Node(vol0)
                h0 = distrib[1][0,1]-distrib[1][0,0]
            if kappaType == 2: 
                vol = getVolumeMap(surfu)
                vol = C.center2Node(vol)
            if niter == 0:
                n = getNormalMap(surfu)
                n = C.normalize(n, vect)
                n = C.center2Node(n)
                n = C.normalize(n, vect)
                n, epsl = modifyNormalWithMetric(surfu, n, algo=1, smoothType=smoothType, eps=eps, nitLocal=nitLocal, kappaType=kappaType, kappaS=kappaS)
                if kappaType == 2:
                    n[1][0,:] = n[1][0,:]*h0*vol0[1][0,:]/(vol[1][0,:]*hloc)
                    n[1][1,:] = n[1][1,:]*h0*vol0[1][0,:]/(vol[1][0,:]*hloc)
                    n[1][2,:] = n[1][2,:]*h0*vol0[1][0,:]/(vol[1][0,:]*hloc)
                if cellNu is not None: epsl[1][:] *= cellNu[1][:]
            else:
                if epsl is None:
                    n = getNormalMap(surfu)
                    n = C.normalize(n, vect)
                    n = C.center2Node(n)
                    n = C.normalize(n, vect)
                    n, epsl = modifyNormalWithMetric(surfu, n, algo=1, smoothType=smoothType, eps=eps, nitLocal=nitLocal, kappaType=kappaType, kappaS=kappaS)
                    if niterType == 0: niterl = niter
                    elif niterType == 1: niterl = niter*(1.+C.getMaxValue(epsl, 'hl')/eps) 
                    else:
                        if k1 < niterK[0]: niterl = k1*1./niterK[0]*niter
                        elif k1 < niterK[1]: niterl = niter
                        else: niterl = max(niter + niter*(niterK[1]-k1)/(niterK[2]-niterK[1]),0) 
                    print('%d: niter=%d'%(k1,niterl))
                    n = getSmoothNormalMap(surfu, niter=niterl, eps=epsl[1], algo=1)
                else:
                    if niterType == 0: niterl = niter
                    elif niterType == 1: niterl = niter*(1.+C.getMaxValue(epsl, 'hl')/eps)
                    else:
                        if k1 < niterK[0]: niterl = k1*1./niterK[0]*niter
                        elif k1 < niterK[1]: niterl = niter
                        else: niterl = max(niter + niter*(niterK[1]-k1)/(niterK[2]-niterK[1]),0) 
                    print('%d: niter=%d'%(k1,niterl))
                    n = getSmoothNormalMap(surfu, niter=niterl, eps=epsl[1], algo=1)
                n, epsl = modifyNormalWithMetric(surfu, n, algo=1, smoothType=smoothType, eps=eps, nitLocal=nitLocal, kappaType=kappaType, kappaS=kappaS)
                if cellNu is not None: epsl[1][:] *= cellNu[1][:]
                if kappaType == 2:
                    n[1][0,:] = n[1][0,:]*h0*vol0[1][0,:]/(vol[1][0,:]*hloc)
                    n[1][1,:] = n[1][1,:]*h0*vol0[1][0,:]/(vol[1][0,:]*hloc)
                    n[1][2,:] = n[1][2,:]*h0*vol0[1][0,:]/(vol[1][0,:]*hloc)

        n[1] = hloc*n[1]

        # modification eventuelle de la hauteur globale
        if KCore.isNamePresent(surfu, 'hsize') != -1:
            hsize = C.extractVars(surfu,['hsize'])
            n[1][0,:] = n[1][0,:] * hsize[1][0,:]
            n[1][1,:] = n[1][1,:] * hsize[1][0,:]
            n[1][2,:] = n[1][2,:] * hsize[1][0,:]
    
        surfu = C.addVars([surfu,n])
        surfu = T.deform(surfu, ['sx','sy','sz'])
        surfu = C.rmVars(surfu, ['sx','sy','sz'])

        kminout = kmax
        for noz in range(nzones):
            coords = listOfCoords[noz]
            ni = coords[2]; nj = coords[3]
            ninj = ni*nj
            indicesU = listOfIndices[noz]
            shift = (k1+1)*ninj
            for ind in range(ninj):
                indu = indicesU[ind]; inds = ind + shift
                coords[1][:,inds] = surfu[1][:,indu]
            listOfCoords[noz] = coords

            if blanking:
                # Cree et propage un cellN aux centres
                subc = T.subzone(coords,(1,1,k1+1),(ni,nj,k1+2))
                vol = getVolumeMap(subc)
                cellN = C.array('cellN',vol[2],vol[3],vol[4])
                cellN[1][0,:] = 1.
                indic = (vol[1][0,:] > -1.e-10)
                cellN[1][0,:] = cellN[1][0,:] * indic[:]
                if cellNp[noz] is not None: cellN[1] = numpy.minimum(cellNp[noz][1], cellN[1])
                else: cellNp[noz] = C.copy(cellN)
                cellNp[noz][1][0,:] = cellN[1][0,:]
                if cellNs[noz] is None: cellNs[noz] = C.copy(cellN)
                else:
                    #print(cellNs[noz])
                    #print(cellN)
                    cellNs[noz] = G.stack(cellNs[noz], cellN)
                    #print(cellNs[noz])
                # modification du lissage pour les points masques
                ni1 = ni-1
                ni2 = max(ni-2,0)
                nj2 = max(nj-2,0)
                for ind in range(ninj):
                    indu = indicesU[ind]
                    # in vertex
                    ii = ind//ni; jj = ind - ii*ni
                    # in centers
                    i0 = min(ii,ni2)
                    i1 = min(ii+1,ni2)
                    j0 = min(jj,nj2)
                    j1 = min(jj+1,nj2)
                    cellNu[1][0,indu] *= 0.25*(cellN[1][0,i0+ni1*j0]+
                                               cellN[1][0,i1+ni1*j0]+
                                               cellN[1][0,i0+ni1*j1]+
                                               cellN[1][0,i1+ni1*j1])
                    cellNu[1][0,indu] = 0.
            elif check == 1:
                subc = T.subzone(coords,(1,1,k1+1),(ni,nj,k1+2))
                vol = getVolumeMap(subc)
                if C.getMinValue(vol,'vol') <= -1.e-10:
                    print("Warning: addNormalLayers: only %d layers created."%(k1))
                    stop = 1; kminout = min(kminout,k1+1)

        if check == 1 and stop == 1:
            for noz in range(nzones):
                coords = listOfCoords[noz]
                ni = coords[2]; nj = coords[3]
                coords = T.subzone(coords,(1,1,1),(ni,nj,kminout))
                listOfCoords[noz] = coords
            return listOfCoords

        k1 += 1
    return listOfCoords

# addNormalLayers pour un array non structure
def addNormalLayersUnstr__(surface, distrib, check=0, niterType=0, niter=0, niterK=[], 
                           smoothType=0, eps=0.4, nitLocal=3, 
                           kappaType=0, kappaS=[0.2,1.6], blanking=False, cellNs=[], algo=0):
    try: import Transform as T; import KCore
    except: raise ImportError("addNormalLayers: requires Converter, Transform modules.")    
    if isinstance(surface[0], list): surf = T.join(surface)
    else: surf = surface
    surf = close(surf); surf = T.reorder(surf, (1,))
    kmax = distrib[1].shape[1] # nb of layers in the normal direction

    if kmax < 2: raise ValueError("addNormalLayers: distribution must contain at least 2 points.")

    vect = ['sx','sy','sz']
    hmin = distrib[1][0,1]-distrib[1][0,0]
    hmax = distrib[1][0,kmax-1]-distrib[1][0,kmax-2]
    hmean = (distrib[1][0,kmax-1]-distrib[1][0,0])/(kmax-1)
    # determination de kb1,kb2
    kb1 =-1; kb2 =-1
    for k1 in range(kmax-1):
        if distrib[1][0,k1+1] >= 0.1*hmean and kb1 == -1: kb1 = k1
        elif distrib[1][0,k1+1] >= 1.*hmean and kb2 == -1: kb2 = k1
    kb2 = max(kb2, kb1+2)
    epsl = None # champ pour le lissage
    cellN = None # champ pour le blanking

    for k1 in range(kmax-1):
        #print("Generating layer %d"%k1)
        hloc = distrib[1][0,k1+1]-distrib[1][0,k1]
        if algo == 0: # algo=0, lissage partout, hauteur 1
            if niter == 0:
                n = getNormalMap(surf)
                n = C.normalize(n, vect)
                n = C.center2Node(n)
                n = C.normalize(n, vect)
                n, epsl = modifyNormalWithMetric(surf, n, algo=0, eps=eps)
            else:
                # ancienne version
                if hmin < 0.01*hmean: # presence de couche limite
                    if k1 < kb1: # pas de lissage ds la couche limite
                        n = getSmoothNormalMap(surf, niter=0, eps=eps)
                        np, epsl = modifyNormalWithMetric(surf, n, algo=0, eps=eps)
                        n[1] = np[1]
                    elif k1 < kb2:
                        beta0 = (float(k1-kb1))/float(kb2-1-kb1)
                        n0 = getSmoothNormalMap(surf,niter=0, eps=eps)
                        n0, epsl = modifyNormalWithMetric(surf, n0, algo=0, eps=eps)
                        n = getSmoothNormalMap(surf, niter=niter, eps=eps)
                        np, epsl = modifyNormalWithMetric(surf, n, algo=0, eps=eps)
                        n[1] = (1-beta0)*n0[1] + beta0*np[1]
                    else: # lissage a fond
                        n = getSmoothNormalMap(surf, niter=niter, eps=eps)
                        np, epsl = modifyNormalWithMetric(surf, n, algo=0, eps=eps)
                        beta0 = float((kmax-2-k1))/float(kmax-2)
                        beta0 = beta0*beta0
                        n[1] = (1-beta0)*n[1] + beta0 *np[1]
                else: # pas de couche limite
                    n = getSmoothNormalMap(surf, niter=niter, eps=eps)
                    np, epsl = modifyNormalWithMetric(surf, n, algo=0, eps=eps)
                    if kmax == 2: beta0 = 0.1
                    else: beta0 = float((kmax-2-k1))/float(kmax-2); beta0 = beta0*beta0
                    n[1] = (1-beta0)*n[1] + beta0 *np[1]           
        else: # algo=1 (nouvelle version)
            if k1 == 0 and kappaType == 2: vol = getVolumeMap(surf)
            if niter == 0:
                n = getNormalMap(surf)
                n = C.normalize(n, vect)
                
                # Add cellN to n
                if cellN is not None and blanking:
                    fake = ['cellN',cellN[1],n[2],n[3]]
                    n = C.addVars([n, fake])
                    generator.extrapWithCellN(surf, n)
                    n = C.extractVars(n, ['sx','sy','sz'])
                    #C.convertArrays2File(surf, 'surf.plt')
                n = C.center2Node(n)
                n = C.normalize(n, vect)
                n, epsl = modifyNormalWithMetric(surf, n, algo=1, smoothType=smoothType, eps=eps, nitLocal=nitLocal, kappaType=kappaType, kappaS=kappaS)
                #n = getSmoothNormalMap(surf, niter=niter, eps=epsl[1], algo=1)
                #n, epsl = modifyNormalWithMetric(surf, n, algo=1, smoothType=smoothType, eps=eps, nitLocal=nitLocal, kappaType=kappaType, kappaS=kappaS)
            else:
                if epsl is None:
                    n = getNormalMap(surf)
                    n = C.normalize(n, vect)
                    n = C.center2Node(n)
                    n = C.normalize(n, vect)
                    n, epsl = modifyNormalWithMetric(surf, n, algo=1, smoothType=smoothType, eps=eps, nitLocal=nitLocal, kappaType=kappaType, kappaS=kappaS)
                    if niterType == 0: niterl = niter
                    elif niterType == 1: niterl = niter*(1.+C.getMaxValue(epsl, 'hl')/eps)
                    else:
                        if k1 < niterK[0]: niterl = k1*1./niterK[0]*niter
                        elif k1 < niterK[1]: niterl = niter
                        else: niterl = max(niter + niter*(niterK[1]-k1)/(niterK[2]-niterK[1]),0) 
                    print('%d: niter=%d'%(k1,niterl))
                    n = getSmoothNormalMap(surf, niter=niterl, eps=epsl[1], cellN=cellN, algo=1)
                else:
                    if niterType == 0: niterl = niter
                    elif niterType == 1: niterl = niter*(1.+C.getMaxValue(epsl, 'hl')/eps)
                    else:
                        if k1 < niterK[0]: niterl = k1*1./niterK[0]*niter
                        elif k1 < niterK[1]: niterl = niter
                        else: niterl = max(niter + niter*(niterK[1]-k1)/(niterK[2]-niterK[1]),0) 
                    print('%d: niter=%d'%(k1,niterl))
                    n = getSmoothNormalMap(surf, niter=niterl, eps=epsl[1], cellN=cellN, algo=1)
                n, epsl = modifyNormalWithMetric(surf, n, algo=1, smoothType=smoothType, eps=eps, nitLocal=nitLocal, kappaType=kappaType, kappaS=kappaS)

        n[1] = hloc*n[1]
        a = grow(surf, n)
        surf = close(surf); surf = T.reorder(surf, (1,))

        if check >= 1:
            vol = getVolumeMap(a)
            volmin = C.getMinValue(vol, 'vol')
            if volmin <= -1.e-10:
                if k1 != 0:
                    print("Warning: addNormalLayers: only %d layers created."%(k1))
                    return m
                else: raise ValueError("addNormalLayers: no layer created.")

        if blanking:
            # Cree un cellN aux centres sur la layer a
            if cellN is None:
                cellN = C.node2Center(a)
                cellN = C.initVars(cellN, 'cellN=1.')
                cellN = C.extractVars(cellN, ['cellN'])
            vol = getVolumeMap(a)
            cellN = C.addVars([cellN, vol])
            cellN = C.initVars(cellN, '{cellN}=minimum({vol}>0, {cellN})')
            cellN = C.extractVars(cellN, ['cellN'])
            generator.blankFirst(a, cellN)
            generator.blankSelf(a, cellN)
            generator.blankSelf(a, cellN)
            generator.blankSelf(a, cellN)
            
        n[0] = 'sx0,sy0,sz0'
        surf = C.addVars([surf,n])
        vectn2 = ['sx0','sy0','sz0']
        surf = T.deform(surf, vectn2)
        surf = C.rmVars(surf, vectn2)

        if k1 == 0: m = a
        else: m = T.join(m, a)
        
        if blanking:
            if k1 == 0: cellNs.append(cellN) 
            else:
                op = ['cellN', None, m[2], cellN[3]]
                op[1] = numpy.concatenate((cellNs[0][1], cellN[1]), axis=1) 
                cellNs[0] = op
                #C.convertArrays2File([a,m], 'check.plt')
                #generator.blankPrev(a, cellN, m, op)
                
    return m

# Fonction retournant la carte d'orthogonalite d'une grille
def getOrthogonalityMap(array):
    """Return the orthogonality map in an array.
    Usage: getOrthogonalityMap(array)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            b.append(generator.getOrthogonalityMap(i))
        return b
    else:
        return generator.getOrthogonalityMap(array)

# Fonction retournant la carte de regularite d'une grille
def getRegularityMap(array):
    """Return the regularity map in an array.
    Usage: getRegularityMap(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(generator.getRegularityMap(i))
        return b
    else:
        return generator.getRegularityMap(array)

def getAngleRegularityMap(array):
    """Return the regularity map in an array.
    Usage: getAngleRegularityMap(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(generator.getAngleRegularityMap(i))
        return b
    else:
        return generator.getAngleRegularityMap(array)
        
# Fonction retournant la carte de qualite d'une maillage TRI 
# (0. pour un triangle degenere. pour un triangle equilateral)
def getTriQualityMap(array):
    """Return a TRI quality measure map in an array.
    Usage: getTriQualityMap(array)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            b.append(generator.getTriQualityMap(i))
        return b
    else:
        return generator.getTriQualityMap(array)
        
# Fonction retournant les qualites min, max et mean d'un maillage TRI 
# (0. pour un triangle degenere, 1. pour un triangle equilateral)
def getTriQualityStat(array):
    """Return the orthogonality stats (min, max and mean) in aTRI mesh.
    Usage: getTriQualityStat(array)"""
    card = getTriQualityMap(array)
    meanv10 = C.getMeanRangeValue(card, 'quality', 0., 0.1)
    meanv90 = C.getMeanRangeValue(card, 'quality', 0.9, 1.)
    meanv = C.getMeanValue(card, 'quality')
    return (meanv10, meanv, meanv90)

#------------------------------------------------------------------------------
# Genere des pyramides ayant pour base les QUAD d'une surface donnee
#------------------------------------------------------------------------------
def quad2Pyra(array, hratio = 1.):
     """Create a set of pyramids from a set of quads.
     Usage: quad2Pyra(array, hratio)"""
     return generator.quad2Pyra(array, hratio)

def getMeshFieldInfo(array, field, critValue, verbose):
    fmin  = 1.e32
    fsum  = 0
    fmax  = -1.
    fcrit = 0
    size  = 0
    info = 'INFO %s: min = %1.2e, max = %1.2e, mean = %1.2e, crit(%s %s %s) = %s cells out of %s | %2.2f%% (%s)'

    DictFunction = {'vol':getVolumeMap, 'orthogonality':getOrthogonalityMap, 'regularity':getRegularityMap, 'regularityAngle':getAngleRegularityMap}

    for cpt, m in enumerate(array):
        f = DictFunction[field](m)[1]

        size_loc  = numpy.size(f)
        fcrit_loc = numpy.count_nonzero(f<critValue) if field == 'vol' else numpy.count_nonzero(f>critValue) 
        fmin_loc  = numpy.min(f)
        fmax_loc  = numpy.max(f)
        fsum_loc  = numpy.sum(f)

        fmin   = min(fmin_loc, fmin)
        fmax   = max(fmax_loc, fmax)
        fsum  += fsum_loc
        fcrit += fcrit_loc
        size  += size_loc

        if verbose == 2 or (verbose == 1 and fcrit_loc > 0):
            print(info%(field.upper(),fmin_loc,fmax_loc,fsum_loc/float(size_loc),field,'<' if field == 'vol' else '>',critValue,fcrit_loc,size_loc,fcrit_loc/float(size_loc)*100,"Zone %d"%(cpt)))


    if verbose == 2 or (verbose == 1 and fcrit_loc > 0):
        print('#'*(len(field)+7))
        print(info%(field.upper(),fmin,fmax,fsum/float(size),field,'<' if field == 'vol' else '>',critValue,fcrit,size,fcrit/float(size)*100,'GLOBAL'))
        print('#'*(len(field)+7)+'\n')

    return fmin, fmax, fsum/float(size), fcrit

def checkMesh(array, critVol=0., critOrtho=15., critReg=0.1, critAngReg=15., addGC=False, verbose=0):
    """Return information on mesh quality."""
    if not isinstance(array[0], list): array = [array] 

    #addGC: dummy argument to match the pyTree function

    vmin,vmax,vmean,vcrit = getMeshFieldInfo(array, 'vol', critVol, verbose)
    omin,omax,omean,ocrit = getMeshFieldInfo(array, 'orthogonality', critOrtho, verbose)
    rmin,rmax,rmean,rcrit = getMeshFieldInfo(array, 'regularity', critReg, verbose)
    amin,amax,amean,acrit = getMeshFieldInfo(array, 'regularityAngle', critAngReg, verbose)

    return {'vmin':vmin,'vmax':vmax,'vmean':vmean,'vcrit':vcrit,
            'rmin':rmin,'rmax':rmax,'rmean':rmean,'rcrit':rcrit,
            'amin':amin,'amax':amax,'amean':amean,'acrit':acrit,
            'omin':omin,'omax':omax,'omean':omean,'ocrit':ocrit}
