"""Grid generation module.
"""
__version__ = '2.9'
__author__ = "Stephanie Peron, Sam Landier, Christophe Benoit, Gaelle Jeanfaivre, Pascal Raud, Luis Bernardos"
# 
# Python Interface to create arrays defining meshes
#
from . import generator
import numpy
def cart(Xo, H, N, api=1):
    """Create a cartesian mesh defined by a structured array.
    Usage: cart((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    return generator.cart(Xo, H, N, api)

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
    if len(arrays) == 0:
        raise ValueError("bbox: array list is empty.")
    if not isinstance(arrays[0], list): ars = [arrays]
    else: ars = arrays
    try: import Converter as C
    except: raise ImportError("bbox: requires Converter module.")
    xmin = 1.e256; ymin = 1.e256; zmin = 1.e256
    xmax =-1.e256; ymax =-1.e256; zmax =-1.e256
    for a in ars:
        varx = KCore.isNamePresent(a, 'CoordinateX')
        if varx == -1:
            varx = KCore.isNamePresent(a, 'x')
            if varx == -1:
                raise ValueError("bbox: x-coordinate not found.")
                return []
            else: varx = 'x'
        else: varx = 'CoordinateX'

        vary = KCore.isNamePresent(a, 'CoordinateY')
        if vary == -1:
            vary = KCore.isNamePresent(a, 'y')
            if vary == -1:
                raise ValueError("bbox: y-coordinate not found.")
                return []
            else: vary = 'y'
        else: vary = 'CoordinateY'

        varz = KCore.isNamePresent(a, 'CoordinateZ')
        if varz == -1:
            varz = KCore.isNamePresent(a, 'z')
            if varz == -1 :
                raise ValueError("bbox: z-coordinate not found.")
                return []
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

def BB(array, method='AABB', weighting=0):
    """Return the axis-aligned or oriented bounding box of an array as an array.
    Usage: b = BB(a, method='AABB', weighting=0)"""
    try: import Converter as C
    except: raise ImportError("BB: requires Converter module.")
    if isinstance(array[0], list):
        out = []
        for a in array:
            if method == 'AABB':  # Computes AABB
                sbb = bbox(a)
                ar = C.array('x,y,z', 2, 2, 2)
                C.setValue(ar, (1,1,1), [sbb[0], sbb[1], sbb[2]])
                C.setValue(ar, (2,1,1), [sbb[3], sbb[1], sbb[2]])
                C.setValue(ar, (1,2,1), [sbb[0], sbb[4], sbb[2]])
                C.setValue(ar, (1,1,2), [sbb[0], sbb[1], sbb[5]])
                C.setValue(ar, (2,2,1), [sbb[3], sbb[4], sbb[2]])
                C.setValue(ar, (1,2,2), [sbb[0], sbb[4], sbb[5]])
                C.setValue(ar, (2,1,2), [sbb[3], sbb[1], sbb[5]])
                C.setValue(ar, (2,2,2), [sbb[3], sbb[4], sbb[5]])
                out.append(ar)
            elif method == 'OBB':  # Computes OBB
                out.append(generator.obbox(a,weighting))
            else:
                print('BB: Warning, method=%s not implemented, making an OBB.'%method)
                out.append(generator.obbox(a,weighting))
        return out
    else:
        if method == 'AABB':  # Computes AABB
            sbb = bbox(array)
            ar = C.array('x,y,z', 2, 2, 2)
            C.setValue(ar, (1,1,1), [sbb[0], sbb[1], sbb[2]])
            C.setValue(ar, (2,1,1), [sbb[3], sbb[1], sbb[2]])
            C.setValue(ar, (1,2,1), [sbb[0], sbb[4], sbb[2]])
            C.setValue(ar, (2,2,1), [sbb[3], sbb[4], sbb[2]])
            C.setValue(ar, (1,1,2), [sbb[0], sbb[1], sbb[5]])
            C.setValue(ar, (2,1,2), [sbb[3], sbb[1], sbb[5]])
            C.setValue(ar, (1,2,2), [sbb[0], sbb[4], sbb[5]])
            C.setValue(ar, (2,2,2), [sbb[3], sbb[4], sbb[5]])
        elif method == 'OBB':  # Computes OBB
            ar = generator.obbox(array,weighting)
        else:
            print('BB: Warning, method=%s not implemented, making an OBB.'%method)
            ar = generator.obbox(array,weighting)
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
            if isinstance(a[1], list): # array2
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
        return generator.bboxIntersection(m1, m2)

def checkPointInCEBB(array, P):
    """Check if point P is in the Cartesian Elements Bounding Box of
    array."""
    return generator.checkPointInCEBB(array, P)

def enforceX(array, x0, enforcedh, N, add=0):
    """Enforce a x0-centered line in a distribution defined by an array.
    Usage: enforceX(array, x0, enforcedh, supp, add) -or-
    Usage: enforceX(array, x0, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforceX(array, x0, enforcedh, N)
    else:
        return generator.enforce(array, "enforceX", x0, enforcedh, N, add)
    
def enforceY(array, y0, enforcedh, N, add=0):
    """Enforce a j line in a distribution defined by an array.
    Usage: enforceY(array, y0, enforcedh, supp, add) -or-
    Usage: enforceY(array, y0, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforceY(array, y0, enforcedh, N)
    else:
        return generator.enforce(array, "enforceY", y0, enforcedh, N, add)
    
def enforceZ(array, z0, enforcedh, N, add=0):
    """Enforce a k line in a distribution defined by an array.
    Usage: enforceZ(array, z0, enforcedh, supp, add)"""
    return generator.enforce(array, "enforceZ", z0, enforcedh, N, add)
    
def enforcePlusX(array, enforcedh, N, add=0):
    """Enforce the first X-line in a distribution defined by an array.
    (one sided distribution, right).
    Usage: enforcePlusX(array, enforcedh, supp, add) -or-
    Usage: enforcePlusX(array, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforcePlusX(array, enforcedh, N)
    else:
        return generator.enforce(array, "enforcePlusX", 0., enforcedh, N, add)
    
def enforcePlusY(array, enforcedh, N, add=0):
    """Enforce a j line in a distribution  defined by an array.
    (one sided distribution, top).
    Usage: enforcePlusY(array, enforcedh, supp, add) -or-
    Usage: enforcePlusY(array, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforcePlusY(array, enforcedh, N)
    else:
        return generator.enforce(array, "enforcePlusY", 0., enforcedh, N, add)

def enforcePlusZ(array, enforcedh, N, add=0):
    """Enforce a k line in a distribution  defined by an array.
    (one sided distribution, top).
    Usage: enforcePlusZ(array, enforcedh, supp, add)"""
    return generator.enforce(array, "enforcePlusZ", 0., enforcedh, N, add)
    
def enforceMoinsX(array, enforcedh, N, add=0):
    """Enforce the last X-line in a distribution (one sided, left).
    Usage: enforceMoinsX(array, enforcedh, supp, add) -or-
    Usage: enforceMoinsX(array, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforceMoinsX(array, enforcedh, N)
    else:
        return generator.enforce(array, "enforceMoinsX", 0., enforcedh, N, add)
    
def enforceMoinsY(array, enforcedh, N, add=0):
    """Enforce a j line in a distribution (one sided distribution, bottom).
    Usage: enforceMoinsY(array, enforcedh, supp, add) -or-
    Usage: enforceMoinsY(array, enforcedh, (supp,add))"""
    if isinstance(N, tuple):
        return generator.enforceMoinsY(array, enforcedh, N)
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
    try: import Geom as D; import math; import Converter as C; import KCore
    except: raise ImportError("enforceCurvature2: requires Converter and Geom modules.")
    
    tol = 1.e-12 # tolerance sur les pts confondues(=close)
    loop = 0 # le contour est il une boucle

    posx = KCore.isNamePresent(arrayC, 'x')
    if (posx != -1): xt = C.extractVars(arrayC, ['x'])[1]
    else:
        posx = KCore.isNamePresent(arrayC, 'CoordinateX')
        if (posx != -1): xt = C.extractVars(arrayC, ['CoordinateX'])[1]
        else: raise ValueError("enforceCurvature2: coordinates must be present in array.")
    posy = KCore.isNamePresent(arrayC, 'y')
    if (posy != -1): yt = C.extractVars(arrayC, ['y'])[1]
    else:
        posy = KCore.isNamePresent(arrayC, 'CoordinateY')
        if (posy != -1): yt = C.extractVars(arrayC, ['CoordinateY'])[1]
        else: raise ValueError("enforceCurvature2: coordinates must be present in array.")
    posz = KCore.isNamePresent(arrayC, 'z')
    if (posz != -1): zt = C.extractVars(arrayC, ['z'])[1]
    else:
        posz = KCore.isNamePresent(arrayC, 'CoordinateZ')
        if (posz != -1): zt = C.extractVars(arrayC, ['CoordinateZ'])[1]
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

def map(array, d, dir=0):
    """Map a distribution on a curve or a surface.
    Usage: map(array, d, dir)"""
    if (len(d) == 5 and d[3] != 1 and d[4] == 1 and dir == 0):
        return map2d(array, d)
    elif (len(d) == 5 and dir != 0):
        return map1dpl(array, d, dir)
    else: return map1d(array, d)
    
# map sur une courbe
def map1d(array, d):
    return generator.map(array, d)

# map par lignes dans la direction dir
def map1dpl(array, d, dir):
    try: import Transform as T; import Converter as C
    except:
        raise ImportError("map: requires Transform and Converter modules.")
    if (dir == 2): m = T.reorder(array, (2,1,3))
    elif (dir == 3): m = T.reorder(array, (3,2,1))
    elif (dir == 1): m = array
    ni = m[2]; nj = m[3]; nk = m[4]; ndi = d[2]; ndi2 = ndi*nj
    a = C.array('x,y,z', ndi, nj, nk)
    for k in range(nk):
        for j in range(nj):
            l = T.subzone(m, (1,j+1,k+1), (ni,j+1,k+1))
            am = map1d(l, d)
            ind = j*ndi+k*ndi2
            a[1][:,ind:ndi+ind] = am[1][:,0:ndi]
    if (dir == 2): a = T.reorder(a, (2,1,3))
    elif (dir == 3): a = T.reorder(a, (3,2,1))
    return a

# map sur une surface
def map2d(array, d):
    try: import Transform as T; import Converter as C
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
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(mapCurvature___(i, N, power, dir))
        return b
    else:
        return mapCurvature___(array, N, power, dir)

def mapCurvature___(array, N, power, dir):
    try: import Transform as T; import Converter as C
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
        import Transform as T; import Converter as C
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

#==============================================================================
# Densifie le maillage d'un i-array
# IN: array: i-array
# IN: h: pas de discretisation
# OUT: i-array avec la nouvelle discretisation
#==============================================================================
def densify(array, h):
    if isinstance(array[0], list):
        l = []
        for i in array:
            l.append(generator.densify(i, h))
        return l
    else:
        return generator.densify(array, h)

def hyper2D(array, arrayd, type):
    """Generate an hyperbolic mesh. 
    Usage: hyper2D(array, arrayd, type)"""
    return generator.hyper2D(array, arrayd, type)

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

def close(array, tol=1.e-12):
    """Close a mesh defined by an array gathering points closer than tol.
    Usage: close(array, tol)"""
    if isinstance(array[0], list):
        return generator.closeAll(array, tol)
    else:
        return generator.close(array, tol)

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
    try: import Transform as T; import Converter as C
    except:
        raise ImportError("plaster: requires Converter, Transform modules.")

    c = C.convertArray2Tetra(contours); c = T.join(c)
    s = C.convertArray2Tetra(surfaces); s = T.join(s)
   
    bb = bbox(contours)
    bb[0] = bb[0] - 1.e-2; bb[1] = bb[1] - 1.e-2; bb[2] = bb[2] - 1.e-2
    bb[3] = bb[3] + 1.e-2; bb[4] = bb[4] + 1.e-2; bb[5] = bb[5] + 1.e-2
    lx = bb[3]-bb[0]; ly = bb[4]-bb[1]; lz = bb[5]-bb[2]
    s1 = lx*ly; s2 = lx*lz; s3 = ly*lz
    if (s1 >= s2 and s1 >= s3):
        if side == 0:
            p = cart( (bb[0],bb[1],bb[2]), (lx/(ni-1),ly/(nj-1),1), (ni,nj,1) )
        else:
            p = cart( (bb[0],bb[1],bb[5]), (lx/(ni-1),ly/(nj-1),1), (ni,nj,1) )
        dir = (0,0,1)
    elif (s2 >= s1 and s2 >= s3):
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

def TFITri(a1, a2, a3):
    """Generate a transfinite interpolation mesh from 3 input curves.
    Usage: TFITri(a1,a2,a3)"""
    from . import TFIs
    return TFIs.TFITri(a1, a2, a3)

def TFIO(a):
    """Generate a transfinite interpolation mesh for 1 input curve.
    Usage: TFIO(a1,a2,a3)"""
    from . import TFIs
    return TFIs.TFIO(a)

def TFIHalfO(a1, a2):
    """Generate a transfinite interpolation mesh for 2 input curves.
    Usage: TFIHalfO(a1,a2)"""
    from . import TFIs
    return TFIs.TFIHalfO(a1, a2)

def TFIMono(a1, a2):
    """Generate a transfinite interpolation mesh for 2 input curves.
    Usage: TFIMono(a1, a2)"""
    from . import TFIs
    return TFIs.TFIMono(a1, a2)

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

def getVolumeMap(array):
    """Return the volume map in an array.
    Usage: getVolumeMap(array)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            b.append(generator.getVolumeMap(i))
        return b
    else:
        return generator.getVolumeMap(array)

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
  
def getSmoothNormalMap(array, niter=2, eps=0.4):
    """Return the map of smoothed and non-normalized surface normals in 
    an array.
    Usage: getSmoothNormalMap(array, niter, eps)"""
    try: import Converter as C
    except: raise ImportError("getSmoothNormalMap: requires Converter module.")
    it = 1
    n = getNormalMap(array)
    n = C.normalize(n, ['sx','sy','sz'])
    n = C.center2Node(n)
    n = C.normalize(n, ['sx','sy','sz'])
    while (it < niter):
        np = C.node2Center(n)
        np = C.normalize(np, ['sx','sy','sz'])
        np = C.center2Node(np)
        np = C.normalize(np, ['sx','sy','sz'])
        it += 1
        n[1][:] += eps*np[1][:]
        n = C.normalize(n, ['sx','sy','sz'])
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

    infos = Collar.createCollarMesh__(surfaces1,surfaces2,distribj,distribk,\
                                      niterj, niterk, ext, alphaRef, type,\
                                      contour, constraints1, constraints2,toldist)
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
def getLocalStepFactor__(s, sn):
    import Converter as C
    sc = C.addVars([s,sn])
    sc = C.node2Center(sc)
    ht = generator.getLocalStepFactor(s, sc)
    niter = 10; eps = 0.01; it = 0
    while it < niter:
        ht = C.node2Center(ht)
        ht = C.center2Node(ht)
        it += 1
    return ht

#==============================================================================
def modifyNormalWithMetric(array, narray):
    """Correct the normals located at nodes with respect to metric."""
    try: import Converter as C
    except: raise ImportError("modifyNormalWithMetric: requires Converter module.")
    a = C.copy(array); n = C.copy(narray)
    if (len(a) == 5): a = C.convertArray2Hexa(a); n = C.convertArray2Hexa(n)
    ht = getLocalStepFactor__(a, n)
    nx = C.extractVars(n,['sx']); ny = C.extractVars(n,['sy']); nz = C.extractVars(n,['sz'])
    nx[1][0,:] = ht[1][0,:]* nx[1][0,:]
    ny[1][0,:] = ht[1][0,:]* ny[1][0,:]
    nz[1][0,:] = ht[1][0,:]* nz[1][0,:]    
    n = C.addVars([n, nx, ny, nz])
    return n

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
        for i in array:
            b.append(generator.getInCircleMap(i))
        return b
    else:
        return generator.getInCircleMap(array)

# -- Add normal layers to a surface --
def addNormalLayers(surface, distrib, check=0, niter=0, eps=0.4):
    """Generate N layers to a surface following normals. distrib is the 
    height of the layers.
    If niter=0, the normal are not smoothed; else niter is the number of
    smoothing iterations applied to normals.
    Usage: addNormalLayers(surface, distrib, check, niter)"""
    if isinstance(surface[0], list): # liste d'arrays
        if len(surface[0]) == 5: type = 0 # structured
        else: type = 1 # unstructured
        for i in range(1,len(surface)):
            if ((len(surface[i]) == 5 and type != 0) or (len(surface[i]) == 4 and type != 1)):
                raise ValueError("addNormalLayers: all the surfaces must be structured or unstructured.")
        if type == 0: return addNormalLayersStruct__(surface, distrib, check, niter, eps)
        else: # NS
            out = []
            for s in surface: out.append(addNormalLayersUnstr__(s, distrib, check, niter, eps))
            return out
    else: # 1 array
        if len(surface) == 5: return addNormalLayersStruct__([surface], distrib, check, niter, eps)[0]
        else: return addNormalLayersUnstr__(surface, distrib, check, niter, eps) # NS

#-----------------------------------------------------------------------------
# Generation de grilles cartesiennes multibloc a partir de:
# grilles de corps: bodies
# h: pas d'espace sur la grille la plus fine
# Dfar: distance d eloignement au corps
# nlvl: nb de points par niveau, sauf sur le niveau le plus fin
# nlvl[0]: niveau le plus grossier
#-----------------------------------------------------------------------------
def gencartmb(bodies, h, Dfar, nlvl):
    import KCore
    try: import Transform as T; import Converter as C
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
        xmax = ref[1][0,0] + (pmax[0]-1)*Href
        ymax = ref[1][1,0] + (pmax[1]-1)*Href
        zmax = ref[1][2,0] + (pmax[2]-1)*Href
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
        if (ref[4]- 2*nb < 3  or  ref[3]- 2*nb < 3 or ref[2]- 2*nb < 3):
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
    ref = cart( (xmin,ymin,zmin), (Href,Href,Href), (Ni,Nj,Nk))
    
    # Niveaux
    out = []
    c = 1
    for i in nlvl :
        lev,ref = createLevel(i, c, ref)
        c  = c + 1
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
        if (array[3] != 1 or array[4] != 1):
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
        import Converter as C; import Transform as T
        import Geom as D
    except:
        raise ImportError("mapSplit: requires Converter, Geom, Transform modules.")

    a = C.copy(array); d = C.copy(dist)

    # Compute total length of dist
    ldt = D.getLength(d); ldti = 1./ldt

    posx = KCore.isNamePresent(d, "x")
    posy = KCore.isNamePresent(d, "y")
    posz = KCore.isNamePresent(d, "z")
    if (posx == -1): posx = KCore.isNamePresent(a, "CoordinateX")
    if (posy == -1): posy = KCore.isNamePresent(a, "CoordinateY")
    if (posz == -1): posz = KCore.isNamePresent(a, "CoordinateZ")
    if (posx == -1 or posy == -1 or posz == -1):
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
    if (lt/stepmin > densMax): stepmin = lt/densMax
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
        if (d1 < d2): indsplit = ind1
        else: indsplit = ind2
        # Split distribution with "split points"
        if (indsplit_previous < indsplit):
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
        import Converter as C
        b = C.convertArray2Tetra(a); b = close(b)
        return generator.T3mesher2D(b, grading, triangulateOnly, metricInterpYype)
    except:
        return generator.T3mesher2D(a, grading, triangulateOnly, metricInterpType)

def tetraMesher(a, maxh=-1., grading=0.4, triangulateOnly=0, 
                remeshBoundaries=0, algo=1):
    """Create a TRI/TETRA mesh given a set of BAR or surfaces in a.
    Usage: tetraMesher(a, maxh, grading)"""
    try:
        import Converter as C; import Transform as T
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
        if algo == 0:
            if remeshBoundaries == 0:
                if maxh < 0: maxh = 1.e6
                return generator.netgen1(a, maxh, grading)
            else:
                try: import Converter as C
                except: raise ImportError("tetraMesher: requires Converter module.")
                C.convertArrays2File([a], 'netgen.0120.stl', 'fmt_stl')
                if maxh < 0: maxh = 1.e6
                return generator.netgen2(a, maxh, grading)
        else:
            # find holes coords
            holes = []
            try: 
                import Transform as T
                sp = T.splitConnexity(a)
                for s in sp:
                    px = C.isNamePresent(s, 'x')
                    if (px != -1): ext = C.extractVars(s, ['x','y','z'])
                    else: ext = C.extractVars(s, ['CoordinateX','CoordinateY','CoordinateZ'])
                    n = getNormalMap(s)
                    n = C.center2Node(n)
                    n = C.normalize(n, ['sx','sy','sz'])
                    pt = C.getValue(ext, 0)
                    n = C.getValue(n, 0)
                    eps = 1.e-10
                    holes.append([pt[0]+eps*n[0], pt[1]+eps*n[1], pt[2]+eps*n[2]])
            except: pass
            return generator.tetgen(a, maxh, grading, holes)
    else:
        raise TypeError("tetraMesher: requires BAR or TRI mesh.")

def fittingPlaster(contour, bumpFactor=0.):
    """Generate a structured points cloud over a BAR.
    Usage: fittingPlaster(contour, bumpFactor)"""
    try:
        import Converter as C
        contour = C.convertArray2Tetra(contour)
    except: pass
    contour = close(contour)
    return generator.fittingPlaster(contour, bumpFactor)

def gapfixer(contour, cloud, hardPoints=None, refine=1):
    """Fix a gap defined by a contour bar and a point cloud representing the gap surface.
    Some hard points can be specified to force the constructed surface to pass by.
    If the optional refine argument is set to 0, the resulting surface will be a contrained triangulation of the contour [and the additional hard nodes].
    Usage: gapFixer(contour, cloud, hardPoints, refine)"""
    try:
        import Converter as C
        contour = C.convertArray2Tetra(contour)
    except: pass
    contour = close(contour)
    return generator.gapfixer(contour, cloud, hardPoints, refine)

def gapsmanager(components, mode=0, refine=0, coplanar=0):
    """Fix a gap between several component surfaces (list of arrays).
    The mode sets the case (0 for POST with a center mesh, 1 for POST with a nodal mesh, 2 otherwise).
    The planar argument tells whether the components are coplanar or not (1 for coplanar case).
    Usage: gapsmanager(components, mode, refine, coplanar)"""
    try:
        import Converter as C
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
    try:
        import Converter as C; import Transform as T
    except:
        raise ImportError("snapFront: requires Converter, Transform module.")

    try:
        import Converter as C; import Transform as T
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

    #import Converter as C
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
# OUT: outList: list which contains :
#               - the surfaces joined in a unique array
#               - the contours (if not nul)
#               - the corners (if not nul)
#------------------------------------------------------------------------------
def refinedSharpEdges__(surfaces, step, angle):
    """Get refined sharp edges from a given surface. 
    Usage: snapSharpEdges(meshes, surfaces)"""
    try:
        import Post as P; import Geom as D; import Converter as C 
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
                                             ext=ext, optimized=optimized)
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

def octree(stlArrays, snearList=[], dfarList=[], dfar=-1., balancing=0, levelMax=1000, ratio=2, octant=None):
    """Generate an octree (or a quadtree) mesh starting from a list of TRI (or BAR) arrays defining bodies,
    a list of corresponding snears, and the extension dfar of the mesh."""
    try: import Converter as C; s = C.convertArray2Tetra(stlArrays)
    except: s = stlArrays
    if ratio == 2:
        o = generator.octree(s, snearList, dfarList, dfar, levelMax, octant)
        if balancing == 0: return o
        elif balancing==1: return balanceOctree__(o, 2, corners=0)
        elif balancing==2: return balanceOctree__(o, 2, corners=1)
    else: #27-tree
        o = generator.octree3(s, snearList, dfar, levelMax, octant)
        if balancing == 0: return o
        else: return balanceOctree__(o,3, corners=0)

def conformOctree3(octree):
    """Conformize an octree3.
    Usage: conformOctree3(octree)"""
    return generator.conformOctree3(octree)
    
def balanceOctree__(octree, ratio=2, corners=0):
    return generator.balanceOctree(octree, ratio, corners)

def extendOctreeGrids__(A, ext, optimized):
    """Extend grids with ext cells. If optimized is ext, the minimum overlapping is ensured.
    Usage: extendOctreeGrids__(cartGrids, ext, optimized)"""
    return generator.extendCartGrids(A, ext, optimized)

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
    try: import Converter as C
    except: raise ImportError("expandLayer: requires Converter module.")
    
    indic = C.node2Center(octreeHexa)
    indic = C.initVars(indic, 'indicator', 0.)
    indic = generator.modifyIndicToExpandLayer(octreeHexa, indic,
                                               level, corners)
    return adaptOctree(octreeHexa, indic, balancing)

# addnormallayers pour une liste d'arrays structures
def addNormalLayersStruct__(surfaces, distrib, check=0, niter=0, eps=0.4):
    import KCore
    try: import Converter as C; import Transform as T
    except: raise ImportError("addNormalLayers: requires Converter, Transform modules.")   
    kmax = distrib[1].shape[1] # nb of layers in the normal direction

    if kmax < 2: raise ValueError("addNormalLayers: distribution must contain at least 2 points.")

    vect = ['sx','sy','sz']
    # verifications 
    for nos in range(len(surfaces)):
        surfs = surfaces[nos]
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
    kb1 =-1; kb2 =-1
    for k1 in range(kmax-1):
        if (distrib[1][0,k1+1] >= 0.1*hmean and kb1 == -1): kb1 = k1
        elif (distrib[1][0,k1+1] >= 1.*hmean and kb2 == -1): kb2 = k1
    kb2 = max(kb2, kb1+2)
    imax = surfu[1].shape[1]
    coordsloc = C.array(surfu[0],imax,2,1) # pour le check
    coordsloc[1][0,0:imax] = surfu[1][0,0:imax]
    coordsloc[1][1,0:imax] = surfu[1][1,0:imax]
    coordsloc[1][2,0:imax] = surfu[1][2,0:imax]
    k1 = 0
    stop = 0
    while k1 < kmax-1 and stop == 0:
        hloc = distrib[1][0,k1+1]-distrib[1][0,k1]
        if niter == 0:
            n = getNormalMap(surfu)
            n = C.normalize(n, vect)
            n = C.center2Node(n)
            n = C.normalize(n, vect)
        else:
            if hmin < 0.01*hmean: # presence de couche limite
                if (k1 < kb1): # pas de lissage ds la couche limite
                    n = getSmoothNormalMap(surfu, niter=0, eps=eps)
                    np = modifyNormalWithMetric(surfu, n)
                    n[1] = np[1]
                elif (k1 < kb2):
                    beta0 = (float(k1-kb1))/float(kb2-1-kb1)
                    n0 = getSmoothNormalMap(surfu, niter=0, eps=eps)
                    n0 = modifyNormalWithMetric(surfu,n0)
                    n = getSmoothNormalMap(surfu, niter=niter, eps=eps)
                    np = modifyNormalWithMetric(surfu,n)
                    n[1] = (1-beta0)*n0[1] + beta0*np[1]
                else: # lissage a fond
                    n = getSmoothNormalMap(surfu, niter=niter, eps=eps)
                    np = modifyNormalWithMetric(surfu,n)
                    beta0 = float((kmax-2-k1))/float(kmax-2)
                    beta0 = beta0*beta0
                    n[1] = (1-beta0)*n[1] + beta0*np[1]
            else: # pas de couche limite
                n = getSmoothNormalMap(surfu, niter=niter, eps=eps)
                np = modifyNormalWithMetric(surfu, n)
                if (kmax == 2): beta0 = 0.1
                else: beta0 = float((kmax-2-k1))/float(kmax-2); beta0 = beta0*beta0
                n[1] = (1-beta0)*n[1] + beta0*np[1]
                
        n[1] = hloc*n[1]
        surfu = C.addVars([surfu,n])
        surfu = T.deform(surfu, ['sx','sy','sz'])
        surfu = C.rmVars(surfu, ['sx','sy','sz'])

        kminout = kmax
        for noz in range(nzones):
            coords = listOfCoords[noz]
            ni=coords[2]; nj=coords[3]
            ninj = ni*nj
            indicesU = listOfIndices[noz]
            shift = (k1+1)*ninj
            for ind in range(ninj):
                indu = indicesU[ind]; inds = ind + shift
                coords[1][:,inds] = surfu[1][:,indu]
            listOfCoords[noz] = coords

            if check == 1: 
                subc = T.subzone(coords,(1,1,k1+1),(ni,nj,k1+2))
                vol = getVolumeMap(subc)
                if C.getMinValue(vol,'vol') <= -1.e-10: 
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
def addNormalLayersUnstr__(surface, distrib, check=0, niter=0, eps=0.4):
    try: import Converter as C; import Transform as T; import KCore
    except: raise ImportError("addNormalLayers: requires Converter, Transform modules.")    
    if isinstance(surface[0], list): surf = T.join(surface)
    else: surf = surface
    surf = close(surf); surf = T.reorder(surf, (1,))
    kmax = distrib[1].shape[1] # nb of layers in the normal direction

    if (kmax < 2): raise ValueError("addNormalLayers: distribution must contain at least 2 points.")

    vect = ['sx','sy','sz']
    hmin = distrib[1][0,1]-distrib[1][0,0]
    hmax = distrib[1][0,kmax-1]-distrib[1][0,kmax-2]
    hmean =  (distrib[1][0,kmax-1]-distrib[1][0,0])/(kmax-1)
    # determination de kb1,kb2
    kb1 =-1; kb2 =-1
    for k1 in range(kmax-1):
        if (distrib[1][0,k1+1] >= 0.1*hmean and kb1 == -1): kb1 = k1
        elif (distrib[1][0,k1+1] >= 1.*hmean and kb2 == -1): kb2 = k1
    kb2 = max(kb2, kb1+2)  
    for k1 in range(kmax-1):
        hloc = distrib[1][0,k1+1]-distrib[1][0,k1]
        if (niter == 0):
            n = getNormalMap(surf)
            n = C.normalize(n, vect)
            n = C.center2Node(n)
            n = C.normalize(n, vect)
        else:
            if (hmin < 0.01*hmean): # presence de couche limite
                if (k1 < kb1): # pas de lissage ds la couche limite
                    n = getSmoothNormalMap(surf, niter=0, eps=eps)
                    np = modifyNormalWithMetric(surf, n)
                    n[1] = np[1]
                    
                elif (k1 < kb2):
                    beta0 = (float(k1-kb1))/float(kb2-1-kb1)
                    n0 = getSmoothNormalMap(surf,niter=0, eps=eps)
                    n0 = modifyNormalWithMetric(surf,n0)
                    n = getSmoothNormalMap(surf,niter=niter, eps=eps)
                    np = modifyNormalWithMetric(surf,n)
                    n[1] = (1-beta0)*n0[1] + beta0*np[1]
                else: # lissage a fond
                    n = getSmoothNormalMap(surf,niter=niter, eps=eps)
                    np = modifyNormalWithMetric(surf,n)
                    beta0 = float((kmax-2-k1))/float(kmax-2)
                    beta0 = beta0*beta0
                    n[1] = (1-beta0)*n[1] + beta0 *np[1]
            else: # pas de couche limite
                n = getSmoothNormalMap(surf, niter=niter, eps=eps)
                np = modifyNormalWithMetric(surf, n)
                if (kmax == 2): beta0 = 0.1
                else: beta0 = float((kmax-2-k1))/float(kmax-2); beta0 = beta0*beta0
                n[1] = (1-beta0)*n[1] + beta0 *np[1]           

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
        
        n[0]='sx0,sy0,sz0'
        surf = C.addVars([surf,n])
        vectn2 = ['sx0','sy0','sz0']
        surf = T.deform(surf, vectn2)
        surf = C.rmVars(surf, vectn2)

        if k1 == 0: m = a
        else: 
            m = T.join(m, a)
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
    try: import Converter as C
    except: raise ImportError("getTriQualityStat: requires Converter module.")
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
