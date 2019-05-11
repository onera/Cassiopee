"""Transformation of arrays.
"""
__version__ = '2.9'
__author__ = "Stephanie Peron, Christophe Benoit, Gaelle Jeanfaivre, Pascal Raud"
# 
# Python Interface to make basic transformations on arrays
#
from . import transform
try: import Converter 
except: raise ImportError("Transform: requires Converter module.")

try: range = xrange
except: pass

__all__ = ['_translate', 'translate', 'addkplane', 'breakElements', 'cart2Cyl', '_cart2Cyl', 'collapse', 
    'computeDeformationVector', '_contract', 'contract', 'cyl2Cart', '_cyl2Cart','deform', 'deformNormals', 'deformPoint', 
    'dual', '_homothety', 'homothety', 'join', 'makeCartesianXYZ', 'makeDirect', 'merge', 'mergeCart', 
    'mergeCartByRefinementLevel', 'oneovern', 'patch', 'perturbate', 'projectAllDirs', 'projectDir', 
    'projectOrtho', 'projectOrthoSmooth', 'projectRay', 'reorder', 'reorderAll', 'rotate', '_scale', 'scale', 
    'smooth', 'splitBAR', 'splitConnexity', 'splitCurvatureAngle', 'splitCurvatureRadius', 'splitManifold', 
    'splitMultiplePts', 'splitNParts', 'splitSharpEdges', 'splitSize', 'splitTBranches', 
    'splitTRI', 'subzone', '_symetrize', 'symetrize', 'deformMesh']

#========================================================================================
# Merge a set of cart grids in A for each refinement level
#========================================================================================
def mergeCartByRefinementLevel(A, sizeMax):
    dhmin = 1.e10 
    allDh = []
    for a in A:
        xt = a[1][0,:]
        dh = xt[1]-xt[0]
        dhmin = min(dhmin,dh)
        allDh.append(dh)

    out = []; ok = 0
    levels={}; level=0
    nzones = len(allDh)
    count = 0
    while ok == 0:
        found = 0
        for noc in range(nzones):
            dh = allDh[noc]
            if dh < 1.2*dhmin and dh > 0.8*dhmin: 
                if level in levels:
                    levels[level].append(noc)
                else:
                    levels[level] = [noc]
                count += 1
                found += 1
        print('Level %d: merging %d zones over %d (Total: %d).'%(level,found,nzones,count))
        if found > 0:
            res = []
            for i in levels[level]: res.append(A[i])            
            res = mergeCart(res, sizeMax)
            out += res
        if count == nzones: ok = 1
        
        dhmin = 2.*dhmin
        level += 1
    return out

def mergeCart(A, sizeMax=1000000000, tol=1.e-10):
    """Merge a list of Cartesian zones using the method of weakest descent.
    Usage: mergeCart(A, sizeMax, tol)"""
    return transform.mergeCart(A, sizeMax, tol)

def merge(A, Ac=[], sizeMax=1000000000, dir=0, tol=1.e-10, alphaRef=180.):
    """Merge a list of matching structured grids.
    Usage: merge(A, Ac, sizeMax, dir, tol, alphaRef)"""
    if len(Ac) != 0 and len(A) != len(Ac):
        raise ValueError("merge: node and center arrays must have the same length.")
    # Tri suivant les types
    STRUCTs = []; BARs = []; TRIs = []; QUADs = []; TETRAs = []
    PYRAs = []; PENTAs = []; HEXAs = []; NGONs = []
    for i in A:
        l = len(i)
        if l == 5: STRUCTs.append(i)
        else:
            elt = i[3]
            if elt == 'BAR': BARs.append(i)
            elif elt == 'TRI': TRIs.append(i)
            elif elt == 'QUAD': QUADs.append(i)
            elif elt == 'TETRA': TETRAs.append(i)
            elif elt == 'PYRA': PYRAs.append(i)
            elif elt == 'PENTA': PENTAs.append(i)
            elif elt == 'HEXA': HEXAs.append(i)
            elif elt == 'NGON': NGONs.append(i)
    STRUCTc = []; BARc = []; TRIc = []; QUADc = []; TETRAc = []
    PYRAc = []; PENTAc = []; HEXAc = []; NGONc = []
    for i in Ac:
        l = len(i)
        if l == 5: STRUCTc.append(i)
        else:
            elt = i[3]
            if elt == 'BAR*': BARc.append(i)
            elif elt == 'TRI*': TRIc.append(i)
            elif elt == 'QUAD*': QUADc.append(i)
            elif elt == 'TETRA*': TETRAc.append(i)
            elif elt == 'PYRA*': PYRAc.append(i)
            elif elt == 'PENTA*': PENTAc.append(i)
            elif elt == 'HEXA*': HEXAc.append(i)
            elif elt == 'NGON*': NGONc.append(i)
    # struct
    ret = []; retc = []
    if len(STRUCTs) > 0:
        r = transform.merge(STRUCTs, STRUCTc, sizeMax, dir, tol, alphaRef)
        if len(STRUCTc) == 0: ret += r
        else: ret += r[0]; retc += r[1]
    if len(BARs) > 0:
        if len(BARc) > 0: 
            r = join(BARs, arrayc=BARc)
            ret += [r[0]]; retc += [r[1]]
        else: ret += [join(BARs)]
    if len(TRIs) > 0:
        if len(TRIc) > 0:
            r = join(TRIs, arrayc=TRIc)
            ret += [r[0]]; retc += [r[1]]
        else: ret += [join(TRIs)]
    if len(QUADs) > 0:
        if len(QUADc) > 0:
            r = join(QUADs, arrayc=QUADc)
            ret += [r[0]]; retc += [r[1]]
        else: ret += [join(QUADs)]
    if len(TETRAs) > 0:
        if len(TETRAc) > 0:
            r = join(TETRAs, arrayc=TETRAc)
            ret += [r[0]]; retc += [r[1]]
        else: ret += [join(TETRAs)]
    if len(PYRAs) > 0:
        if len(PYRAc) > 0:
            r = join(PYRAs, arrayc=PYRAc)
            ret += [r[0]]; retc += [r[1]]
        else: ret += [join(PYRAs)]
    if len(PENTAs) > 0:
        if len(PENTAc) > 0:
            r = join(PENTAs, arrayc=PENTAc)
            ret += [r[0]]; retc += [r[1]]
        else: ret += [join(PENTAs)]
    if len(NGONs) > 0:
        if len(NGONc) > 0:
            r = join(NGONs, arrayc=NGONc)
            ret += [r[0]]; retc += [r[1]]
        else: ret += [join(NGONs)]
    if len(Ac) > 0: return (ret, retc)
    else: return ret

def cart2Cyl(a, center=(0,0,0), axis=(0,0,1)):
    """Transform a mesh defined in Cartesian coordinates into cylindrical coordinates.
    Usage: cart2Cyl(a, center, axis)"""
    b = Converter.copy(a)
    _cart2Cyl(b, center, axis)
    return b

def _cart2Cyl(a, center=(0,0,0), axis=(0,0,1)):
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(transform._cart2CylA(i, center, axis))
        return b
    else:
        return transform._cart2CylA(a, center, axis)

def cyl2Cart(a, center=(0,0,0), axis=(0,0,1)):
    """Transform a mesh defined in Cylindrical coordinates into cartesian coordinates.
    Usage: cyl2Cart(a, center, axis)"""
    b = Converter.copy(a)
    _cyl2Cart(b, center, axis)
    return b

def _cyl2Cart(a, center=(0,0,0), axis=(0,0,1)):
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(transform._cyl2CartA(i, center, axis))
        return b
    else:
        return transform._cyl2CartA(a, center, axis)

def collapse(a):
    """Collapse the smallest edge of each element for TRI arrays. Return a BAR.
    Usage: collapse(a)"""
    if isinstance(a[0], list): 
        b = []
        for i in a:
            b.append(transform.collapse(i))
        return b
    else:
        return transform.collapse(a)
    
def translate(a, transvect):
    """Translate a grid.
    Usage: translate(a, (v1,v2,v3))"""
    b = Converter.copy(a)
    _translate(b, transvect)
    return b

def _translate(a, transvect):
    if isinstance(a[0], list):
        for i in a:
            transform.translate(i, transvect)
    else:
        transform.translate(a, transvect)
    return None

def rotate(a, center, arg1, arg2=None, 
           vectors=[['VelocityX','VelocityY','VelocityZ'],['MomentumX','MomentumY','MomentumZ']]):
    """Rotate a grid."""
    if arg2 is None: # kind of euler angles
        return rotate3__(a, center, arg1, vectors)
    elif isinstance(arg2, float) or isinstance(arg2, int):
        return rotate1__(a, center, arg1, arg2, vectors)
    else: return rotate2__(a, center, arg1, arg2, vectors)
    
def rotate1__(array, center, rotvect, angle, vectors): # centre+axe+angle+vecteurs a modifier
    """Rotate a mesh defined by an array around vector n of center Xc
    and of angle teta.
    Usage: rotate(a, (xc,yc,zc), (nx,ny,nz), teta, vectors)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            b.append(transform.rotateA1(i, center, rotvect, angle, vectors))
        return b
    else:
        return transform.rotateA1(array, center, rotvect, angle, vectors)

def rotate2__(array, center, e1, e2, vectors): # centre+axe1->axe2
    """Rotate a mesh defined by an array of center Xc
    transforming a unitary vector e1 into e2.
    Usage: rotate2(a, (xc,yc,zc), (e1x,e1y,e1z), (e2x,e2y,e2z), vectors)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            b.append(transform.rotateA2(i, center, e1, e2, vectors))
        return b
    else:
        return transform.rotateA2(array, center, e1, e2, vectors)

def rotate3__(array, center, angles, vectors): # centre+3 angles+ champs vectoriels a modifier
    """Rotate a mesh defined by an array of center Xc
    and a set of euler angles.
    Usage: rotate3(a, (xc,yc,zc), (alpha,beta,gamma), vectors)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            b.append(transform.rotateA3(i, center, angles, vectors))
        return b
    else:
        return transform.rotateA3(array, center, angles, vectors)

def homothety(a, center, alpha):
    """Make for a mesh defined by an array an homothety of center Xc and
    of factor alpha.
    Usage: homothety(a, (xc,yc,zc), alpha)"""
    b = Converter.copy(a)
    _homothety(b, center, alpha)
    return b

def _homothety(a, center, alpha):
    if isinstance(a[0], list): 
        for i in a:
            transform.homothety(i, center, alpha)
    else:
        transform.homothety(a, center, alpha)
    return None

def contract(a, center, dir1, dir2, alpha):
    """Contract a mesh around a plane defined by (center, dir1, dir2) and of factor alpha.
    Usage: contract(a, (xc,yc,zc), dir1, dir2, alpha)"""
    b = Converter.copy(a)
    _contract(b, center, dir1, dir2, alpha)
    return b

def _contract(a, center, dir1, dir2, alpha):
    if isinstance(a[0], list): 
        for i in a:
            transform.contract(i, center, dir1, dir2, alpha)
    else:
        transform.contract(a, center, dir1, dir2, alpha)
    return None

def scale(a, factor=1.):
    """Scale a mesh following factor (constant) or (f1,f2,f3) following dir.
    Usage: scale(a, 1)"""
    b = Converter.copy(a)
    _scale(b, factor)
    return b

def _scale(a, factor=1.):
    X = (0,0,0)
    try: import Generator; X = Generator.barycenter(a)
    except: pass
    if isinstance(factor, list) or isinstance(factor, tuple):
        if len(factor) == 1:
            _homothety(a, X, factor)
        elif len(factor) == 3:
            axe1 = (1,0,0); axe2 = (0,1,0); axe3 = (0,0,1)
            _contract(a, X, axe2, axe3, factor[0])
            _contract(a, X, axe1, axe3, factor[1])
            _contract(a, X, axe1, axe2, factor[2])
    else:
        _homothety(a, X, factor)
    return None

def symetrize(a, point, vector1, vector2):
    """Make a symetry of mesh from plane passing by point and of director vector: vector1 and vector2.
    Usage: symetrize(a, (xc,yc,zc), (v1x,v1y,v1z), (v2x,v2y,v2z))"""
    b = Converter.copy(a)
    _symetrize(b, point, vector1, vector2)
    return b

def _symetrize(a, point, vector1, vector2):
    """Make a symetry of mesh from plane passing by point and of director vector: vector1 and vector2.
    Usage: symetrize(a, (xc,yc,zc), (v1x,v1y,v1z), (v2x,v2y,v2z))"""
    if isinstance(a[0], list):
        for i in a:
            transform.symetrize(i, point, vector1, vector2)
    else:
        return transform.symetrize(a, point, vector1, vector2)
    return None

def perturbate(a, radius, dim=3):
    """Perturbate a mesh randomly of radius
    Usage: perturbate(a, radius, dim)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(transform.perturbate(i, radius, dim))
        return b
    else:
        return transform.perturbate(a, radius, dim)

def smooth(a, eps=0.5, niter=4, type=0, fixedConstraints=[], 
           projConstraints=[], delta=1., point=(0,0,0), radius=-1.):
    """Smooth a mesh with a Laplacian.
    Usage: smooth(a, eps, niter, type, fixedConstraints, projConstraints, delta, (xR,yR,zR), radius)"""
    import KCore; import numpy
    try: import Generator as G
    except:
        raise ImportError("smooth: requires Converter and Generator modules.")
    if len(fixedConstraints) != 0:
        try:
            fixedConstraints = Converter.convertArray2Tetra(fixedConstraints)
        except: fixedConstraint = []
        else:
            if isinstance(fixedConstraints[0], list):
                fixedConstraint = join(fixedConstraints)
            else: fixedConstraint = fixedConstraints
    else: fixedConstraint = []

    if len(projConstraints) != 0:
        try:
            projConstraints = Converter.convertArray2Tetra(projConstraints)
        except: projConstraint = []
        else:
            if isinstance(projConstraints[0], list):
                projConstraint = join(projConstraints)
            else: projConstraint = projConstraints    
    else: projConstraint = []

    listeType = -1
    if isinstance(a[0], list):
        for i in a:
            if len(i) == 4 and listeType == -1: listeType = 0
            if len(i) == 4 and listeType != 0:
                raise TypeError("smooth: list of zones must be all structured or all unstructured.")
            if len(i) == 5 and listeType == -1: listeType = 1
            if len(i) == 5 and listeType != 1:
                raise TypeError("smooth: list of zones must be all structured or all unstructured.")

        if listeType == 1: # all struct    
            b = Converter.convertArray2Hexa(a); b = join(b); b = G.close(b)
            listOfIndices = KCore.indiceStruct2Unstr2(a, b, 1.e-14)
            c = transform.smooth(b, eps, niter, type, 
                                 fixedConstraint, projConstraint, delta,
                                 point, radius)
            listOfCoords = []
            for noz in range(len(a)):
                coords = Converter.copy(a[noz])
                ninjnk = coords[2]*coords[3]*coords[4]
                indicesU = listOfIndices[noz]
                coords[1][0:3,:] = c[1][0:3,indicesU[:]]
                listOfCoords.append(coords)
            return listOfCoords
        else: # all unstruct
            coords = []
            for i in a:
                coords.append(transform.smooth(i, eps, niter, type, 
                                               fixedConstraint, projConstraint,
                                               delta, point, radius))
            return coords
        
    elif len(a) == 5: # array structure    
        b = Converter.convertArray2Hexa(a); b = G.close(b)
        c = transform.smooth(b, eps, niter, type, fixedConstraint, 
                             projConstraint, delta, point, radius)
        listOfIndices = []
        listOfIndices = KCore.indiceStruct2Unstr2([a], b, 1.e-14)
        coords = Converter.copy(a)
        c = transform.smooth(b, eps, niter, type, fixedConstraint, 
                             projConstraint, delta, point, radius)
        ninjnk = coords[2]*coords[3]*coords[4]
        indicesU = listOfIndices[0]
        coords[1][0:3,:] = c[1][0:3,indicesU[:]]
        return coords
    else: # array non structure
        return transform.smooth(a, eps, niter, type, fixedConstraint, 
                                projConstraint, delta, point, radius)

def projectAllDirs(arrays, surfaces, vect=['nx','ny','nz'], oriented=0):
    """Project points defined in arrays to surfaces according to the direction provided by vect.
    Usage: projectAllDirs(arrays,surfaces,vect,oriented)
    """
    try:
        b = Converter.convertArray2Tetra(surfaces) 
    except: b = surfaces
    if not isinstance(arrays[0], list):
        return transform.projectAllDirs([arrays], b, vect, oriented)[0] 
    else: 
        return transform.projectAllDirs(arrays, b, vect, oriented)

def projectDir(surfaces, arrays, dir, smooth=0, oriented=0):
    """Project surfaces onto surface arrays following dir. 
    Usage: projectDir(surfaces, arrays, dir)"""
    try:
        b = Converter.convertArray2Tetra(arrays) 
        if isinstance(b[0], list): b = join(b)
    except: b = arrays[0]
    if isinstance(surfaces[0], list):
        if smooth == 0: return transform.projectDir(surfaces, b, dir, oriented)
        else: return transform.projectSmoothDir(surfaces, b, dir, oriented)
    else:
        if smooth == 0: return transform.projectDir([surfaces], b, dir, oriented)[0]
        else: return transform.projectSmoothDir([surfaces], b, dir, oriented)[0]

def projectOrtho(surfaces, arrays):
    """Project a list of zones surfaces onto surface arrays following normals. 
    Usage: projectOrtho(surfaces, arrays)"""
    try:
        b = Converter.convertArray2Tetra(arrays)
        if isinstance(b[0], list): b = join(b)
    except: b = arrays[0]
    if isinstance(surfaces[0], list):
        return transform.projectOrtho(surfaces, b)
    else:
        return transform.projectOrtho([surfaces], b)[0]

def projectOrthoSmooth(surfaces, arrays, niter=1):
    """Project a list of zones surfaces onto surface arrays following normals. 
    Usage: projectOrthoSmooth(surfaces, arrays)"""
    if isinstance(surfaces[0], list): surfs = surfaces
    else: surfs = [surfaces]

    # Projection orthogonale directe
    a = projectOrtho(surfs, arrays)
    # Calcul du vecteur normal
    for i in range(len(surfs)):
        a[i][1][:] = surfs[i][1][:]-a[i][1][:]
        a[i][0] = a[i][0].replace('x','nx')
        a[i][0] = a[i][0].replace('y','ny')
        a[i][0] = a[i][0].replace('z','nz')
        a[i][0] = a[i][0].replace('CoordinateX','nx')
        a[i][0] = a[i][0].replace('CoordinateY','ny')
        a[i][0] = a[i][0].replace('CoordinateZ','nz')
   
    # Lissage du vecteur
    n = a; vect = ['nx','ny','nz']
    for i in range(niter):
        #n = Cpnverter.normalize(n, vect)
        #for k in n:
        #    if len(k) == 5: transform.extrapInside(k) # dark hack
        n = Converter.node2ExtCenter(n)
        #n = Converter.normalize(n, vect)
        n = Converter.extCenter2Node(n)

    for i in range(len(surfs)):
        surfs[i] = Converter.addVars([surfs[i], n[i]])
    
    # Projection
    a = projectAllDirs(surfs, arrays, vect)
    a = Converter.rmVars(a, vect)
    if isinstance(surfaces[0], list): return a
    else: return a[0]

def projectRay(surfaces, arrays, P):
    """Project surfaces onto surface arrays using rays starting from P. 
    Usage: projectRay(surfaces, arrays, P)"""
    try:
        b = Converter.convertArray2Tetra(arrays)
        if isinstance(b[0], list): b = join(b)
    except: b = arrays[0]
    if isinstance(surfaces[0], list):
        return transform.projectRay(surfaces, b, P)
    else:
        return transform.projectRay([surfaces], b, P)[0]

def deform(a, vector=['dx','dy','dz']):
    """Deform surface by moving surface of the vector dx, dy, dz. 
    Usage: deform(a, vector)"""    
    if len(vector) != 3 or not isinstance(vector[1],str):
        a = Converter.addVars([a,vector])
        if isinstance(vector[0], list):
            vector = vector[0][0].split(',')
        else: vector = vector[0].split(',')
    else :
        if len(vector) != 3: 
            raise ValueError("deform: 3 variables are required.")

    if isinstance(a[0], list):
        out = []
        for i in a: out+=[transform.deform(i, vector)]
        return out
    else: return transform.deform(a, vector)
    
def deformNormals(array, alpha, niter=1):
    """Deform a a surface of alpha times the surface normals.
    Usage: deformNormals(array, alpha, niter)"""
    try: import Generator as G
    except: raise ImportError("deformNormals: requires Generator module.")   
    if isinstance(array[0], list) and isinstance(alpha[0], list): 
        if len(array) != len(alpha): raise ValueError("deformNormals: number of arrays in a and alpha must be equal.")
        b = []; noi = 0
        for i in array:
            npts = i[1].shape[1]
            alpi = alpha[noi][1][0]/niter
            if npts != len(alpi): raise ValueError("deformNormals: array and alpha must be of same length.")
            aloc = i
            for ite in range(niter):
                n = G.getSmoothNormalMap(aloc, niter=0)
                n[1][:,:] = n[1][:,:]*alpi[:]
                aloc = Converter.addVars([aloc,n])
                aloc = deform(aloc,['sx','sy','sz'])
                aloc = Converter.rmVars(aloc,['sx','sy','sz'])
            b.append(aloc)
            noi += 1
        return b
    
    elif (not isinstance(array[0], list) and not isinstance(alpha[0], list)):
        if array[1].shape[1] != alpha[1].shape[1]: raise ValueError("deformNormals: array and alpha must be of same length.")
        npts = array[1].shape[1]
        alp = alpha[1][0]/niter
        aloc = array
        for ite in range(niter):
            n = G.getSmoothNormalMap(aloc, niter=0)
            n[1][:,:] = n[1][:,:]*alp[:]
            aloc = Converter.addVars([aloc,n])
            aloc = deform(aloc,['sx','sy','sz'])
            aloc = Converter.rmVars(aloc,['sx','sy','sz'])            
        return aloc
    
    else: raise ValueError("deformNormals: array and alpha must be both an array or a list of arrays.")

def deformPoint(a, xyz, dxdydz, depth, width):
    """Deform mesh by moving point (x,y,z) of a vector (dx, dy, dz). depth
    is the depth along vector and width controls the width of deformation.
    Usage: deform(a, (x,y,z), (dx,dy,dz), depth, width)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(transform.deformPoint(i, xyz, dxdydz, depth, width))
        return b
    else:
        return transform.deformPoint(a, xyz, dxdydz, depth, width)

def deformMesh(a, surfDelta, beta=4., type='nearest'):
    """Deform a mesh wrt surfDelta defining surface grids and deformation vector on it.
    Usage: deformMesh(a, surfDelta, beta, type)"""
    if not isinstance(a[0], list):
        if len(a) == 5: 
            if type=='nearest': return deformMeshStruct1__(a, surfDelta, beta)
            elif type=='gridline': return deformMeshStruct2__(a, surfDelta, beta)
            else: raise TypeError("deformMesh: type not valid.")
        else: raise TypeError("deformMesh: not valid for unstructured arrays.")
    else:
        out = []
        for i in a:
            if len(i) == 5: 
                if type=='nearest': out.append(deformMeshStruct1__(i, surfDelta, beta))
                elif type =='gridline': out.append(deformMeshStruct2__(i, surfDelta, beta))
                else: raise TypeError("deformMesh: type not valid.")
            else: raise TypeError("deformMesh: not valid for unstructured arrays.")
        return out

def deformMeshStruct1__(arrayi, surfDelta, beta):
    try: import Generator as G
    except: raise ImportError("deformMesh: requires Converter and Generator modules.")
    array = Converter.copy(arrayi)
    surfDelta = Converter.convertArray2Tetra(surfDelta)

    ni = array[2]; nj = array[3]; nk = array[4]
    if ((ni == 1) and (nj == 1)) or ((ni == 1) and (nk == 1)) or ((nj == 1) and (nk == 1)):
        raise TypeError("deformMesh: not for 1D-arrays.")
        
    borders = []
    dim = 3
    if ni == 1 or nj == 1 or nk == 1: dim = 2
    if dim == 3:
        m1 = subzone(array, (1,1,1),(1,nj,nk)); m1[2] = nj; m1[3] = nk; m1[4] = 1 
        m2 = subzone(array,(ni,1,1),(ni,nj,nk)); m2[2] = nj; m2[3] = nk; m2[4] = 1 
        m3 = subzone(array, (1,1,1),(ni,1,nk)); m3[3] = nk; m3[4] = 1 
        m4 = subzone(array,(1,nj,1),(ni,nj,nk)); m4[3] = nk; m4[4] = 1 
        m5 = subzone(array,(1,1,1),(ni,nj,1))
        m6 = subzone(array,(1,1,nk),(ni,nj,nk))
        borders = [m1,m2,m3,m4,m5,m6]
    else:
        m1 = subzone(array, (1,1,1),(1,nj,nk)); m1[2] = nj; m1[3] = 1 
        m2 = subzone(array,(ni,1,1),(ni,nj,nk)); m2[2] = nj; m2[3] = 1
        m3 = subzone(array, (1,1,1),(ni,1,nk))
        m4 = subzone(array,(1,nj,1),(ni,nj,nk))
        borders = [m1,m2,m3,m4]
    res = computeDeformationVector(borders, surfDelta, beta)
    borders = Converter.addVars([borders,res])
    delta = G.TFI(borders)
    delta = Converter.extractVars(delta, ['dx','dy','dz'])
    dim = delta[1].shape[0] 
    array[1][:dim,:] += delta[1][:dim,:]
    return array

def deformMeshStruct2__(arrayi, surfDelta, beta):
    try: import Generator as G
    except: raise ImportError("deformMesh: requires Converter and Generator modules.")
    array = Converter.copy(arrayi)
    surfDelta = Converter.convertArray2Tetra(surfDelta)
    if not isinstance(surfDelta[0], list): surfDelta = [surfDelta]

    ni = array[2]; nj = array[3]; nk = array[4]
    if (ni==1 and nj==1) or (ni==1 and nk==1) or (nj==1 and nk==1):
        raise TypeError("deformMesh: not for 1D-arrays.")
    
    res = transform.deformMeshStruct(array, surfDelta, beta)
    res = Converter.addVars([array,res])

    borders = []
    dim = 3
    if ni == 1 or nj == 1 or nk == 1: dim = 2
    if dim == 3:
        m1 = subzone(res,(1,1,1),(1,nj,nk)); m1[2]=nj; m1[3]=nk; m1[4]=1 
        m2 = subzone(res,(ni,1,1),(ni,nj,nk)); m2[2]=nj; m2[3]=nk; m2[4]=1 
        m3 = subzone(res,(1,1,1),(ni,1,nk)); m3[3]=nk; m3[4]=1 
        m4 = subzone(res,(1,nj,1),(ni,nj,nk)); m4[3]=nk; m4[4]=1 
        m5 = subzone(res,(1,1,1),(ni,nj,1))
        m6 = subzone(res,(1,1,nk),(ni,nj,nk))
        borders = [m1,m2,m3,m4,m5,m6]
    else:
        m1 = subzone(res,(1,1,1),(1,nj,nk)); m1[2]=nj; m1[3]=1 
        m2 = subzone(res,(ni,1,1),(ni,nj,nk)); m2[2]=nj; m2[3]=1
        m3 = subzone(res,(1,1,1),(ni,1,nk))
        m4 = subzone(res,(1,nj,1),(ni,nj,nk))
        borders = [m1,m2,m3,m4]

    del res
    delta = G.TFI(borders)
    delta = Converter.extractVars(delta, ['dx','dy','dz'])
    dim = delta[1].shape[0] 
    array[1][:dim,:] += delta[1][:dim,:]
    return array

def computeDeformationVector(array, surfDelta, beta=4.):
    """Computes a deformation vector for each border of a mesh
    Usage: computeDeformationVector(array, delta)"""
    try:
        surfDelta = Converter.convertArray2Tetra(surfDelta)
    except: pass
    if not isinstance(surfDelta[0], list): surfDelta = [surfDelta]
    if isinstance(array[0], list):
        return transform.computeDeformationVector(array, surfDelta, beta)
    else: return transform.computeDeformationVector([array], surfDelta, beta)  

def join(array, array2=[], arrayc=[], arrayc2=[], tol=1.e-10):
    """Join two arrays in one or join a list of arrays in one. 
    Usage: join(array, array2) or join(arrays)"""
    if arrayc == []:
        if array2 == []: return joing__(array, tol)
        else: return joins__(array, array2, tol)
    else:
        if array2 == []: return joingb__(array, arrayc, tol)
        else: return joinsb__(array, array2, arrayc, arrayc2, tol)
            
def joing__(arrays, tol):
    if len(arrays) > 1: a = arrays[0]
    elif len(arrays) == 1: return arrays[0]
    else: return []
    if len(a) == 4 and a[3] != 'NGON': return transform.joinAll(arrays, tol)
    pool = arrays[:]
    pool.pop(0)
    while len(pool) > 0:
        success = 0
        c = 0
        for i in pool:
            try:
                a = joins__(a, i, tol)
                pool.pop(c)
                success = 1
                break
            except: pass
            c += 1
        if success == 0: raise ValueError("join: cannot join!")
    return a

def joingb__(arrays, arraysc, tol):
    if len(arrays) != len(arraysc): raise ValueError("join: arrays and arraysc must be of same length.")
    if len(arrays) > 1: a = arrays[0]; ac = arraysc[0]
    elif len(arrays) == 1: return arrays[0],arraysc[0]
    else: return []
    if len(a) == 4 and a[3] != 'NGON': return transform.joinAllBoth(arrays, arraysc, tol)
    pool = arrays[:]; poolc = arraysc[:]
    pool.pop(0); poolc.pop(0)
    while len(pool) > 0:
        success = 0; c = 0
        for noi in range(len(pool)):
            try:
                a,ac = joinsb__(a, pool[noi], ac, poolc[noi], tol)
                pool.pop(c); poolc.pop(c)
                success = 1
                break
            except: pass
            c += 1
        if success == 0: raise ValueError("join: cannot join!")
    return [a,ac]

def joinsb__(array1, array2, arrayc1, arrayc2, tol):
    if len(array1) == 5 and len(array2) == 5:
        return transform.joinBoth( array1, array2, arrayc1, arrayc2, tol)
    
    elif len(array1) == 4 and len(array2) == 5:
        if array1[3] == "NGON":
            a = Converter.convertArray2NGon(array2)
            ac = Converter.convertArray2NGon(arrayc2)
        else: 
            a = Converter.convertArray2Hexa(array2)
            ac = Converter.convertArray2Hexa(arrayc2)
        cn = a[2]; ac = [ac[0], ac[1], cn, ac[3]+'*']
        return transform.joinBoth(array1, a, arrayc1, ac, tol)

    elif len(array1) == 5 and len(array2) == 4:
        if array2[3] == "NGON":
            a = Converter.convertArray2NGon(array1)
            ac = Converter.convertArray2NGon(arrayc1)
        else: 
            a = Converter.convertArray2Hexa(array1)
            ac = Converter.convertArray2Hexa(arrayc1)
        cn = a[2]; ac = [ac[0], ac[1], cn, ac[3]+'*']# la conversion d'un array structure en centres en non structure ne donne pas un array de type elt*
        return transform.joinBoth(a, array2, ac, arrayc2, tol)

    else:
        if array1[3] == "NGON" and array2[3] != "NGON":
            a = Converter.convertArray2NGon(array2)
            ac = Converter.convertArray2NGon(arrayc2)
            cn = a[2]; ac = [ac[0], ac[1], cn, ac[3]+'*']# la conversion d un array structure en centres en non structure ne donne pas un array de type elt*
            return transform.joinBoth(array1, a, arrayc1, ac, tol)
      
        elif array1[3] != "NGON" and array2[3] == "NGON":
            a = Converter.convertArray2NGon(array1)
            ac = Converter.convertArray2NGon(arrayc1)
            cn = a[2]; ac = [ac[0], ac[1], cn, ac[3]+'*']# la conversion d'un array structure en centres en non structure ne donne pas un array de type elt*
            return transform.join(a, array2, tol)
        else: return transform.joinBoth(array1, array2, arrayc1, arrayc2, tol)

def joins__(array1, array2, tol):
    if len(array1) == 5 and len(array2) == 5:
        return transform.join(array1, array2, tol)

    elif len(array1) == 4 and len(array2) == 5:
        if array1[3] == "NGON": a = Converter.convertArray2NGon(array2)
        else: a = Converter.convertArray2Hexa(array2)
        return transform.join(array1, a, tol)

    elif len(array1) == 5 and len(array2) == 4:
        if array2[3] == "NGON": a = Converter.convertArray2NGon(array1)
        else: a = Converter.convertArray2Hexa(array1)
        return transform.join(a, array2, tol)

    else: # Unstruct           
        if array1[3] == "NGON" and array2[3] != "NGON":
            a = Converter.convertArray2NGon(array2)
            return transform.join(array1, a, tol)
        elif (array1[3] != "NGON" and array2[3] == "NGON"):
            a = Converter.convertArray2NGon(array1)
            return transform.join(a, array2, tol)
        else: return transform.join(array1, array2, tol)

def patch(a1, a2, position=None, nodes=None):
    """Patch mesh2 defined by a2 in mesh1 defined by a1 at position (i,j,k).
    Usage: patch(a1, a2, (i,j,k))"""
    import numpy
    if (isinstance(a1[0], list) or isinstance(a2[0], list)):
        raise TypeError("patch: not for a list of arrays.")
    if position is None and nodes is None:
        raise TypeError("patch: either position or nodes must be defined.")
    if position is not None and nodes is not None:
        raise TypeError("patch: position and nodes can not be both defined.")
    if position is not None:
        return transform.patch(a2, a1, position)
    elif nodes is not None:
        if isinstance(nodes, list):
            nodes = numpy.asarray(nodes,numpy.int32, order='F')
        return transform.patch2(a2, a1, nodes)

def oneovern(a, N, add=1):
    """Take one over N points from mesh.
    Usage: oneovern(a, (Ni,Nj,Nk))"""
    if isinstance(a[0], list): 
        b = []
        for i in a:
            b.append(transform.oneovern(i, N, add))
        return b
    else:
        return transform.oneovern(a, N, add)

def subzone(array, minIndex, maxIndex=None, type=None):
    """Take a subzone of mesh.
    Usage: subzone(array, (imin,jmin,kmin), (imax,jmax,kmax))"""
    if maxIndex is None: # non structure
        if type == 'elements':
            if len(array) == 5: 
                raise TypeError("subzone with a list of elements not yet implemented for structured arrays.")
            return transform.subzoneElements(array, minIndex)
        elif type == 'faces':
            if len(array) == 5:
                return transform.subzoneStructInt(array,minIndex)
            else:
                return transform.subzoneFaces(array, minIndex)
        elif type == 'nodes':
            if len(array) == 5: 
                raise TypeError("subzone with a list of nodes not yet implemented for structured arrays.")
            return  transform.subzoneUnstruct(array, minIndex)
        else: 
            if len(array) == 5:
                raise TypeError("subzone with a list of nodes not yet implemented for structured arrays.")
            return transform.subzoneUnstruct(array, minIndex)
    else: # structure (subzone par range)
        if len(array) == 4:
            raise TypeError("subzone with two ranges is not valid for unstructured arrays.")
        return transform.subzoneStruct(array, minIndex, maxIndex)

def reorder(a, order):
    """Reorder the numerotation of mesh.
    Usage: reorder(a, (2,1,-3))"""
    if isinstance(a[0], list): 
        b = []
        for i in a:
            b.append(transform.reorder(i, order))
        return b
    else:
        return transform.reorder(a, order)     

def reorderAll(arrays, dir=1):
    """Orientate normals of all surface blocks consistently in one direction (1) or the opposite (-1).
    For unstructured inputs, when dir is set to 1(-1), it means outward(inward).
    Usage: reorderAll(arrays, dir)"""
    btype = 0
    arrs = []
    if isinstance(arrays[0], list): arrs = arrays
    else: arrs.append(arrays)
    arrays = arrs
    
    for i in arrays:
      if btype == 0: btype = len(i)
      elif btype != len(i): raise TypeError("reorderAll: mixed blocks (structured or unstructured) is not supported.")

    if btype == 5: return transform.reorderAll(arrays, dir)
    elif btype == 4: return transform.reorderAllUnstr(arrays, dir)
    else: raise TypeError("reorderAll: blocks types are not supported.")

def makeCartesianXYZ(a):
    """Reorder a Cartesian mesh in order to get i,j,k aligned with X,Y,Z."""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(makeCartesianXYZ__(i))
        return b
    else:
        return makeCartesianXYZ__(a)

def makeCartesianXYZ__(z):
    if len(z) == 5:
        import KCore
        ni = z[2]; nj = z[3]; nk = z[4]
        ind = 0; indi = ind+1; indj = ind+ni; indk = ind+ni*nj

        posx = KCore.isCoordinateXPresent(z)
        posy = KCore.isCoordinateYPresent(z)
        posz = KCore.isCoordinateZPresent(z)
        valind = Converter.getValue(z,ind)
        valindi = Converter.getValue(z,indi)
        valindj = Converter.getValue(z,indj)
        valindk = Converter.getValue(z,indk)

        dx_i = valindi[posx]-valind[posx]
        dy_i = valindj[posx]-valind[posx]
        dz_i = valindk[posx]-valind[posx]
        diri = 1; dirj = 2; dirk = 3
        if abs(dx_i) > 0.: diri = 1
        elif abs(dy_i) > 0.: diri = 2
        elif abs(dz_i) > 0.: diri = 3
        dx_j = valindi[posy]-valind[posy]
        dy_j = valindj[posy]-valind[posy]
        dz_j = valindk[posy]-valind[posy]
        if abs(dx_j) > 0.: dirj = 1
        elif abs(dy_j) > 0.: dirj = 2
        elif abs(dz_j) > 0.: dirj = 3
        dx_k = valindi[posz]-valind[posz]
        dy_k = valindj[posz]-valind[posz]
        dz_k = valindk[posz]-valind[posz]
        if abs(dx_k) > 0.: dirk = 1
        elif abs(dy_k) > 0.: dirk = 2
        elif abs(dz_k) > 0.: dirk = 3
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
        z = reorder(z,(dirs[0], dirs[1], dirs[2]))
        ni = z[2]; nj = z[3]; nk = z[4]
        ind = 0; indi = ind+1; indj = ind+ni; indk = ind+ni*nj
        diri = 1; dirj = 1; dirk = 1
        valind = Converter.getValue(z,ind)
        valindi = Converter.getValue(z,indi)
        valindj = Converter.getValue(z,indj)
        valindk = Converter.getValue(z,indk)
        dx_i = valindi[posx]-valind[posx]
        ok = 0
        diri = 1; dirj = 2; dirk = 3
        if dx_i < 0.: diri =-1; ok = 1
        dy_j = valindj[posy]-valind[posy] 
        if dy_j < 0.: dirj =-2; ok = 1
        dz_k = valindk[posz]-valind[posz] 
        if dz_k < 0.: dirk =-3; ok = 1
        if ok == 1: z = reorder(z,(diri,dirj,dirk))
    return z

def makeDirect(a):
    """Reorder a structured mesh to make it direct."""
    if isinstance(a[0], list): 
        b = []
        for i in a:
            b.append(makeDirect__(i))
        return b
    else:
        return makeDirect__(a)     
  
def makeDirect__(a):
    import KCore
    import KCore.Vector as Vector
    l = len(a)
    if len(a) == 4: return Converter.copy(a)
    px = KCore.isNamePresent(a, 'x')
    py = KCore.isNamePresent(a, 'y')
    pz = KCore.isNamePresent(a, 'z')
    if px == -1 or py == -1 or pz == -1: return Converter.copy(a)
    ni = a[2]; nj = a[3]; nk = a[4]; p = a[1]
    i = ni/2; j = nj/2; k = nk/2
    ip1 = max(i+1,ni-1); jp1 = max(j+1,nj-1); kp1 = max(k+1,nk-1)
    ind = int(i + j*ni + k*ni*nj)
    P0 = [ p[px,ind], p[py,ind], p[pz,ind] ]
    ind = int(ip1 + j*ni + k*ni*nj)
    P1 = [ p[px,ind], p[py,ind], p[pz,ind] ]
    ind = int(i + jp1*ni + k*ni*nj)
    P2 = [ p[px,ind], p[py,ind], p[pz,ind] ]
    ind = int(i + j*ni + kp1*ni*nj)
    P3 = [ p[px,ind], p[py,ind], p[pz,ind] ]
    l1 = Vector.sub(P1,P0); ln1 = Vector.norm2(l1)
    l2 = Vector.sub(P2,P0); ln2 = Vector.norm2(l2)
    l3 = Vector.sub(P3,P0); ln3 = Vector.norm2(l3)
    if ln1 > 0 and ln2 > 0 and ln3 > 0:
        c = Vector.cross(l1,l2)
        c = Vector.dot(c,l3)
        if c < 0: b = reorder(a, (1,2,-3)); return b
    return Converter.copy(a)

def addkplane(a, N=1):
    """Add N k-plane(s) to a mesh.
    Usage: addkplane(a, N)"""
    if isinstance(a[0], list): 
        b = []
        for i in a:
            c = addkplane__(i, N)
            b.append(c)
        return b
    else:
        return addkplane__(a, N)
 
# IN: a: surface array
def addkplane__(a, N):
    if len(a) == 5: # structure
        res = a
        for j in range(N): res = transform.addkplane(res)
        return res
    else: # non structure
        try: import Generator as G
        except: return transform.addkplane(a)
        res = []
        for j in range(N):
            b = translate(a, (0,0,j*1.))
            b = transform.addkplane(b)
            if res == []: res = b
            else: res = join(res, b); res = G.close(res)
        return res

def addkplaneCenters(arrayC, arrayK, N=1):
    if isinstance(arrayC[0], list): 
        b = []
        for noi in range(len(arrayC)):
            c = transform.addkplaneCenters(arrayC[noi], arrayK[noi],N)
            b.append(c)
        return b
    else:
        return transform.addkplaneCenters(arrayC, arrayK,N)

# Essaie de couper en 2 en respectant level niveaux de multigrille
def findMGSplit__(n, level):
    ns = (n+1)//2
    if level == 0: return ns
    power = 2**level
    if ((ns-1)%power == 0 and (n-ns)%power == 0): return ns
    if ((ns-2)%power == 0 and (n-ns+1)%power == 0): return ns-1
    if ((ns-3)%power == 0 and (n-ns+2)%power == 0): return ns-2
    if ((ns-4)%power == 0 and (n-ns+3)%power == 0): return ns-3
    return -1

# Fait un split nv et le reste en respectant le multigrille
def findMGSplitUp__(n, nv, level):
    ns = nv
    if ns < 4: ns = 4
    if level == 0: return ns
    power = 2**level
    if ((ns-1)%power == 0 and (n-ns)%power == 0): return ns
    if ((ns-2)%power == 0 and (n-ns+1)%power == 0): return ns-1
    if ((ns-3)%power == 0 and (n-ns+2)%power == 0): return ns-2
    if ((ns-4)%power == 0 and (n-ns+3)%power == 0): return ns-3
    return -1

# Get split dir
def getSplitDir__(ni, nj, nk, dirs):
    dirl = 1
    if ni >= nj and ni >= nk:
        dirl = 1
        if 1 in dirs: dirl = 1
        else:
            if nj >= nk:
                if 2 in dirs: dirl = 2
                elif 3 in dirs: dirl = 3
            else:
                if 3 in dirs: dirl = 3
                elif 2 in dirs: dirl = 2
    elif nj >= ni and nj >= nk:
        dirl = 2
        if 2 in dirs: dirl = 2
        else:
            if ni >= nk:
                if 1 in dirs: dirl = 1
                elif 3 in dirs: dirl = 3
            else:
                if 3 in dirs: dirl = 3
                elif 1 in dirs: dirl = 1
    elif nk >= ni and nk >= nj:
        dirl = 3
        if 3 in dirs: dirl = 3
        else:
            if ni >= nj:
                if 1 in dirs: dirl = 1
                elif 2 in dirs: dirl = 2
            else:
                if 2 in dirs: dirl = 2
                elif 1 in dirs: dirl = 1
    return dirl

# Split size au milieu
def splitSize__(a, N, multigrid, dirs):
    if len(a) == 4: # unstructured
        print('Warning: splitSize: unstructured array not treated.')
        return [a]
    if len(a) == 5: # structured
        ni = a[2]; nj = a[3]; nk = a[4]
        if ni*nj*nk > N:
            dirl = getSplitDir__(ni, nj, nk, dirs)

            if dirl == 1:
                ns = findMGSplit__(ni, level=multigrid)
                if ns > 0: 
                    a1 = subzone(a, (1,1,1), (ns,nj,nk))
                    a2 = subzone(a, (ns,1,1), (ni,nj,nk))
                else: return [a]
            elif dirl == 2:
                ns = findMGSplit__(nj, level=multigrid)
                if ns > 0:
                    a1 = subzone(a, (1,1,1), (ni,ns,nk))
                    a2 = subzone(a, (1,ns,1), (ni,nj,nk))
                else: return [a]
            elif dirl == 3:
                ns = findMGSplit__(nk, level=multigrid)
                if ns > 0:
                    a1 = subzone(a, (1,1,1), (ni,nj,ns))
                    a2 = subzone(a, (1,1,ns), (ni,nj,nk))
                else: return [a]
            else:
                ns = findMGSplit__(ni, level=multigrid)
                if ns > 0:
                    a1 = subzone(a, (1,1,1), (ns,nj,nk))
                    a2 = subzone(a, (ns,1,1), (ni,nj,nk))
                else: return [a]
            l1 = splitSize__(a1, N, multigrid, dirs)
            l2 = splitSize__(a2, N, multigrid, dirs)
            return l1+l2
        else: return [a]

# Split size decentre
def splitSizeUp__(a, N, multigrid, dirs):
    if len(a) == 4: # unstructured
        print('Warning: splitSize: unstructured zone not treated.')
        return [a]
    if len(a) == 5: # structured
        ni = a[2]; nj = a[3]; nk = a[4]
        nij = ni*nj; nik = ni*nk; njk = nj*nk
        if ni*nj*nk > N:
            dirl = getSplitDir__(ni, nj, nk, dirs)
            if dirl == 1:
                nc = N//njk
                ns = findMGSplitUp__(ni, nc, level=multigrid)
                if ns > 0: 
                    a1 = subzone(a, (1,1,1), (ns,nj,nk))
                    a2 = subzone(a, (ns,1,1), (ni,nj,nk))
                else: return [a]
            elif dirl == 2:
                nc = N//nik
                ns = findMGSplitUp__(nj, nc, level=multigrid)
                if ns > 0:
                    a1 = subzone(a, (1,1,1), (ni,ns,nk))
                    a2 = subzone(a, (1,ns,1), (ni,nj,nk))
                else: return [a]
            elif dirl == 3:
                nc = N//nik
                ns = findMGSplitUp__(nk, nc, level=multigrid)
                if ns > 0:
                    a1 = subzone(a, (1,1,1), (ni,nj,ns))
                    a2 = subzone(a, (1,1,ns), (ni,nj,nk))
                else: return [a]
            else:
                nc = N//njk
                ns = findMGSplitUp__(ni, nc, level=multigrid)
                if ns > 0:
                    a1 = subzone(a, (1,1,1), (ns,nj,nk))
                    a2 = subzone(a, (ns,1,1), (ni,nj,nk))
                else: return [a]
            l1 = splitSizeUp__(a1, N, multigrid, dirs)
            l2 = splitSizeUp__(a2, N, multigrid, dirs)
            return l1+l2
        else: return [a]
    
# Split size decentre avec ressources
def splitSizeUpR__(A, N, R, multigrid, dirs, minPtsPerDir):
    # Cree le pool
    SP = []; Nl = 0
    for i in A:
        if len(i) == 5: # structure
            SP.append((Converter.getNCells(i),i)); Nl += Converter.getNCells(i)
    if N == 0: N = Nl*1. / R
    #print 'average cells ', N
    from operator import itemgetter
    
    # Init le vecteur des ressources
    Rs = [0]*R
    mins = minPtsPerDir-1 # nbre de cellules mini des blocs

    out = []
    while len(SP) > 0:
        SP = sorted(SP, key=itemgetter(0), reverse=True)
        Rs = sorted(Rs)
        #print 'ress', Rs[0], Converter.getNCells(SP[0][1])
        a = SP[0][1] # le plus gros
        ni = a[2]; nj = a[3]; nk = a[4]
        ni1 = max(1, ni-1); nj1 = max(1, nj-1); nk1 = max(1, nk-1)
        nik = ni1*nk1; njk = nj1*nk1; nij = ni1*nj1
        Nr = min(N, N-Rs[0])
        ncells = Converter.getNCells(a)
        if ncells > Nr:
            # Calcul le meilleur split
            nc = int(round(Nr*1./njk,0))+1
            ns = findMGSplitUp__(ni, nc, level=multigrid)
            if ns-1 < mins: ns = 5
            delta1 = ns-1
            delta2 = ni-ns
            if delta2 < mins: delta2 -= 1.e6
            delta3 = abs((ns-1)*njk - Nr)/njk
            deltai = delta3-delta1-delta2
            nc = int(round(Nr*1./nik,0))+1
            ns = findMGSplitUp__(nj, nc, level=multigrid)
            if ns-1 < mins: ns = 5
            delta1 = ns-1
            delta2 = nj-ns
            if delta2 < mins: delta2 -= 1.e6
            delta3 = abs(ni1*(ns-1)*nk1 - Nr)/nik
            deltaj = delta3-delta1-delta2
            nc = int(round(Nr*1./nij,0))+1
            ns = findMGSplitUp__(nk, nc, level=multigrid)
            if (ns-1 < mins): ns = 5
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
                ns = findMGSplitUp__(ni, nc, level=multigrid)
                if ns-1 >= mins and ni-ns >= mins:
                    a1 = subzone(a, (1,1,1), (ns,nj,nk))
                    a2 = subzone(a, (ns,1,1), (ni,nj,nk))
                    SP[0] = (Converter.getNCells(a2), a2)
                    out += [a1]; Rs[0] += Converter.getNCells(a1)
                    trynext = 0
            elif dirl == 2:
                nc = int(round(Nr*1./nik,0))+1
                ns = findMGSplitUp__(nj, nc, level=multigrid)
                if ns-1 >= mins and nj-ns >= mins:
                    a1 = subzone(a, (1,1,1), (ni,ns,nk))
                    a2 = subzone(a, (1,ns,1), (ni,nj,nk))
                    SP[0] = ( Converter.getNCells(a2), a2 )
                    out += [a1]; Rs[0] += Converter.getNCells(a1)
                    trynext = 0
            elif dirl == 3:
                nc = int(round(Nr*1./nij,0))+1
                ns = findMGSplitUp__(nk, nc, level=multigrid)
                if ns-1 >= mins and nk-ns >= mins:
                    a1 = subzone(a, (1,1,1), (ni,nj,ns))
                    a2 = subzone(a, (1,1,ns), (ni,nj,nk))
                    SP[0] = ( Converter.getNCells(a2), a2 )
                    out += [a1]; Rs[0] += Converter.getNCells(a1)
                    trynext = 0
            if trynext == 1:
                out += [a]; Rs[0] += Converter.getNCells(a); del SP[0]
        else:
            out += [a]; Rs[0] += Converter.getNCells(a); del SP[0]
    #print 'ress:', Rs
    #Tot = 0
    #for i in Rs: Tot += i
    #print 'Tot', Tot
    return out

def splitSize(array, N=0, multigrid=0, dirs=[1,2,3], type=0, R=None, 
              minPtsPerDir=5):
    """Split a block until it has less than N points.
    Usage: splitSize(array, N, multigrid=0, dirs=[1,2,3], type=0)"""
    minPtsPerDir = max(minPtsPerDir, 2**(multigrid+1)+1)
    if R is not None: type = 2
    if type == 0: # middle split
        if isinstance(array[0], list):
            b = []
            for i in array: b += splitSize__(i, N, multigrid, dirs)
            return b
        else: return splitSize__(array, N, multigrid, dirs)
    elif type == 1: # upwind split
        if isinstance(array[0], list):
            b = []
            for i in array: b += splitSizeUp__(i, N, multigrid, dirs)
            return b
        else: return splitSizeUp__(array, N, multigrid, dirs)
    else: # decentre avec ressources (greedy)
        if isinstance(array[0], list):
            return splitSizeUpR__(array, N, R, multigrid, dirs, minPtsPerDir)
        else: return splitSizeUpR__([array], N, R, multigrid, dirs)
        
#==============================================================================
# find splits pour splitNParts
# IN: ni,nj,nk de l'array a decouper
# IN: N nbre de blocs a obtenir pour cet array
# IN: dirs: directions autorisees pour la decoupe de ce bloc
# IN: multigrid: niveau de multigrille a respecter
# OUT: liste des splits a effectuer [(dir, n)] de taille N
#==============================================================================
def findSplits__(ni, nj, nk, N, dirs, multigrid):
    ldir = len(dirs)
    # Passage en multigrille
    plev = 2**multigrid
    nig = (ni-1)//plev+1; njg = (nj-1)//plev+1; nkg = (nk-1)//plev+1
    out = []
    if ldir == 1: # pas le choix
        if dirs[0] == 1:
            #ns = nig/N
            ns = round(nig*1./N, 0); ns = int(ns)
            r = (ni-N*ns*plev)/plev
            #print ns, ni-(N-1)*ns*plev
            b1 = 1
            for j in range(N):
                if r > 0: b2 = b1+plev*(ns+1); r -= 1
                elif r < 0: b2 = b1+plev*(ns-1); r += 1
                else: b2 = b1+plev*ns
                if j == N-1: b2 = ni
                out.append((b1,b2,1,nj,1,nk))
                b1 = b2
  
        elif dirs[0] == 2: 
            ns = round(njg/N, 0); ns = int(ns)
            r = (nj-N*ns*plev)/plev
            b1 = 1
            for j in range(N):
                if r > 0: b2 = b1+plev*(ns+1); r -= 1
                elif r < 0: b2 = b1+plev*(ns-1); r += 1
                else: b2 = b1+plev*ns
                if j == N-1: b2 = nj
                out.append((1,ni,b1,b2,1,nk))
                b1 = b2
        else: 
            ns = round(nkg/N, 0); ns = int(ns)
            r = (nk-N*ns*plev)/plev
            b1 = 1
            for j in range(N):
                if r > 0: b2 = b1+plev*(ns+1); r -= 1
                elif r < 0: b2 = b1+plev*(ns-1); r += 1
                else: b2 = b1+plev*ns
                if j == N-1: b2 = nk
                out.append((1,ni,1,nj,b1,b2))
                b1 = b2
                
    elif ldir == 2:
        if dirs[0] == 1: ns1 = ni; bs1 = ni; ng1 = (ni-1)//plev+1 
        elif dirs[0] == 2: ns1 = nj; bs1 = nj; ng1 = (nj-1)//plev+1 
        else: ns1 = nk; bs1 = nk; ng1 = (nk-1)//plev+1 
        if dirs[1] == 1: ns2 = ni; bs2 = ni; ng2 = (ni-1)//plev+1 
        elif dirs[1] == 2: ns2 = nj; bs2 = nj; ng2 = (nj-1)//plev+1 
        else: ns2 = nk; bs2 = nk; ng2 = (nk-1)//plev+1 
        best = [1,1,1]
        size = -1
        for N1 in range(1,N+1):
            for N2 in range(1,N+1):
                if N1*N2 == N:
                    ns1 = ng1//N1; ns2 = ng2//N2
                    s = min(ns1, ns2)
                    if s > size: best = [N1,N2]; size = s
        N1 = best[0]; N2 = best[1]
        #ns1 = ng1/N1; ns2 = ng2/N2
        ns1 = round(ng1*1./N1, 0); ns1 = int(ns1)
        ns2 = round(ng2*1./N2, 0); ns2 = int(ns2)
        r1 = (ni-N1*ns1*plev)/plev
        r2 = (nj-N2*ns2*plev)/plev
        i1 = 1; i2 = ni
        j1 = 1; j2 = nj
        k1 = 1; k2 = nk
        b1 = 1
        for i in range(N1): # tous les splits en 1
            if r1 > 0: b2 = b1+plev*(ns1+1); r1 -= 1
            elif r1 < 0: b2 = b1+plev*(ns1-1); r1 += 1
            else: b2 = b1+plev*ns1
            if dirs[0] == 1: 
                i1 = b1; i2 = b2; 
                if i == N1-1: i2 = ni
            elif dirs[0] == 2: 
                j1 = b1; j2 = b2
                if i == N1-1: j2 = nj
            else: 
                k1 = b1; k2 = b2     
                if i == N1-1: k2 = nk

            r2 = (ni-N2*ns2*plev)/plev
            c1 = 1
            for j in range(N2): # tous les splits en 2
                if (r2 > 0): c2 = c1+plev*(ns2+1); r2 -= 1
                elif (r2 < 0): c2 = c1+plev*(ns2-1); r2 += 1
                else: c2 = c1+plev*ns2
                if (dirs[1] == 1): 
                    i1 = c1; i2 = c2
                    if j == N2-1: i2 = ni
                elif (dirs[1] == 2): 
                    j1 = c1; j2 = c2
                    if j == N2-1: j2 = nj
                else: 
                    k1 = c1; k2 = c2
                    if j == N2-1: k2 = nk
                out.append( (i1,i2,j1,j2,k1,k2) )
                c1 = c2
            b1 = b2
            
    else: # ldir == 3
        ns1 = ni; bs1 = ni; ng1 = (ni-1)/plev+1
        ns2 = nj; bs2 = nj; ng2 = (nj-1)/plev+1
        ns3 = nk; bs3 = nk; ng3 = (nk-1)/plev+1
        best = [1,1,1]
        size = -1
        for N1 in range(1,N+1):
            for N2 in range(1,N+1):
                for N3 in range(1,N+1):
                    if N1*N2*N3 == N:          
                        ns1 = ng1/N1; ns2 = ng2/N2; ns3 = ng3/N3
                        s = min(ns1, ns2, ns3)
                        if s > size:
                            best = [N1,N2,N3]; size = s
                        elif s == size:
                            if N1 == best[0]: # discrimine suivant 2/3
                                sl = min(ns2, ns3)
                                sb = min(ng2/best[1], ng3/best[2])
                                if sl > sb: best = [N1,N2,N3]
                            elif (N2 == best[1]): # discrimine suivant 1/3
                                sl = min(ns1, ns3)
                                sb = min(ng1/best[0], ng3/best[2])
                                if (sl > sb): best = [N1,N2,N3]
                            else:  # discrimine suivant 1/2
                                sl = min(ns1, ns2)
                                sb = min(ng1/best[0], ng2/best[1])
                                if (sl > sb): best = [N1,N2,N3]

        N1 = best[0]; N2 = best[1]; N3 = best[2]
        #ns1 = ng1/N1; ns2 = ng2/N2; ns3 = ng3/N3
        ns1 = round(ng1*1./N1, 0); ns1 = int(ns1)
        ns2 = round(ng2*1./N2, 0); ns2 = int(ns2)
        ns3 = round(ng3*1./N3, 0); ns3 = int(ns3)
        r1 = (ni-N1*ns1*plev)/plev
        r2 = (nj-N2*ns2*plev)/plev
        r3 = (nk-N3*ns3*plev)/plev
        b1 = 1
        for i in range(N1): # tous les splits en i
            if r1 > 0: b2 = b1+plev*(ns1+1); r1 -= 1
            elif r1 < 0: b2 = b1+plev*(ns1-1); r1 += 1
            else: b2 = b1+plev*ns1
            i1 = b1; i2 = b2
            if i == N1-1: i2 = ni
            r2 = (nj-N2*ns2*plev)/plev
            c1 = 1
            for j in range(N2): # tous les splits en j
                if r2 > 0: c2 = c1+plev*(ns2+1); r2 -= 1
                elif r2 < 0: c2 = c1+plev*(ns2-1); r2 += 1
                else: c2 = c1+plev*ns2
                j1 = c1; j2 = c2
                if j == N2-1: j2 = nj
                r3 = (nk-N3*ns3*plev)/plev
                d1 = 1
                for k in range(N3): # tous les splits en k                    
                    if r3 > 0: d2 = d1+plev*(ns3+1); r3 -= 1
                    elif r3 < 0: d2 = d1+plev*(ns3-1); r3 += 1
                    else: d2 = d1+plev*ns3
                    k1 = d1; k2 = d2
                    if k == N3-1: k2 = nk
                    out.append( (i1,i2,j1,j2,k1,k2) )
                    d1 = d2
                c1 = c2
            b1 = b2
    return out

#==============================================================================
# IN: l: nbre de zones
# IN: N: nbre total de blocs voulu
# IN: Np: nbre de pts pour chaque zone
# OUT: Ns: nbre de blocs a obtenir par zones
#==============================================================================
def findNsi__(l, N, Np):
    Sum = 0.
    for i in range(l): Sum += Np[i]
    Nm = Sum *1. / N
    Nm = max(Nm, 1.)
    Ns = [0]*l # nbre de splits a effectuer
    Er = [0]*l # Erreur de split
    for i in range(l): Ns[i] = Np[i]*1. / Nm
    # Passage en entier
    for i in range(l):
        val = round(Ns[i], 0)
        if val == 0: Ns[i] = 1
        else: Ns[i] = int(val)
        Er[i] = (Ns[i]*Nm - Np[i],i)

    # Tri suivant Er
    from operator import itemgetter
    Er = sorted(Er, key=itemgetter(0))

    # Check for N
    ND = 0
    for i in range(l): ND += Ns[i]
    while ND != N:
        #print 'Round ', ND, N
        if ND < N: # pas assez de blocs
            # On cherche a augmenter les splits des plus grands Er
            for i in range(N-ND):
                e = Er[i]
                no = e[1]
                Ns[no] += 1 #print Ns[no]*Nm-Np[no]
            pass
        elif ND > N: # trop de blocs
            # On cherche a diminuer les splits des plus petits Er 
            for i in range(ND-N):
                e = Er[l-i-1]
                no = e[1]
                Ns[no] -= 1 #print Ns[no]*Nm-Np[no]

        ND = 0
        for i in range(l): ND += Ns[i]
        #print 'Final Round ', ND, N
    return Ns

# split une liste d'arrays structures en N parties a peu pres egales
def splitNParts(arrays, N, multigrid=0, dirs=[1,2,3]):
    """Split blocks in N blocks."""
    if not isinstance(arrays[0], list): arrays = [arrays]
    # Fait des paquets de zones structurees et NGON
    arraysS = []; arraysN = []
    NpS = []; NpN = [] # nbre de points
    NeS = []; NeN = [] # nbre de cellules
    outO = []
    for a in arrays:
        if len(a) == 5:
            arraysS.append(a)
            NpS.append(a[2]*a[3]*a[4])
            NeS.append(max(a[2]-1,1)*max(a[3]-1,1)*max(a[4]-1,1))
        elif a[3] == 'NGON': 
            arraysN.append(a)
            NpN.append(a[1].size)
            c = a[2]; NeN.append(c[0,c[0,1]+2])
        else:
            arraysN.append(a)
            NpN.append(a[1].size)
            c = a[2]; NeN.append(c.shape[1])

    SumS = 0.; SumN = 0.
    for i in NeS: SumS += i
    for i in NeN: SumN += i
    if SumS+SumN < 0.01: return outO
    alpha = N*1./(SumS+SumN)
    NbN = len(NeN) # nbre de grilles non structurees
    NPart = [0]*(NbN+1); Nt = 0
    for i in range(NbN): NPart[i] = max(int(alpha*NeN[i]),1); Nt += NPart[i]
    if SumS != 0: NPart[NbN] = max(N-Nt,1)
    else: NPart[NbN-1] = max(N-Nt+NPart[NbN-1],1)

    # Blocs non structures
    outN = []
    for i in range(len(arraysN)):
        a = arraysN[i]
        if NPart[i] > 1:
            if a[3] == 'NGON': 
                elts = transform.splitNGon(a, NPart[i])
            else:
                elts = transform.splitElement(a, NPart[i])
            for e in elts: outN.append(subzone(a, e, type='elements'))
        else: outN.append(a)

    # Blocs structures
    l = len(arraysS)
    if l == 0: return outN+outO
    NPa = NPart[NbN]
    Ns = findNsi__(l, NPa, NpS)

    outS = []
    for i in range(l):
        a = arraysS[i]
        ni = a[2]; nj = a[3]; nk = a[4]
        splits = findSplits__(ni, nj, nk, Ns[i], dirs, multigrid)
        for j in splits:
            a1 = subzone(a, (j[0],j[2],j[4]), (j[1],j[3],j[5]))
            outS.append(a1)
    return outS+outN+outO

def splitCurvatureRadius__(array, Rs):
    out = transform.splitCurvatureRadius(array, Rs)
    n = len(out)-1
    if (n >= 1):
        try:
            f = join(out[n], out[0])
            ret = transform.splitCurvatureRadius(f, Rs)
            if (len(ret) == 1):
                out[0] = f; del out[n]
            return out
        except: return out
    else: return out
    
def splitCurvatureRadius(a, Rs=100.):
    """Return the indices of the array where the curvature radius is low.
    Usage: splitCurvatureRadius(a, Rs)"""
    if isinstance(a[0], list):
        out = []
        for i in a:
            out += splitCurvatureRadius__(i, Rs)
        return out
    else: return splitCurvatureRadius__(a, Rs)
    
def splitConnexity(a):
    """Split array into connex zones.
    Usage: splitConnexity(a)"""
    if isinstance(a[0], list):
        out = []
        for i in a:
            out += transform.splitConnexity(i)
        return out
    return transform.splitConnexity(a)

def breakElements(a):
    """Break an array (in general NGON) in a set of arrays 
    of BAR, TRI, ... elements
    Usage: breakElements(a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            if len(i) == 4 and i[3] == 'NGON':
                b += transform.breakElements(i)
            else: b.append(i)
        return b
    else:
        if len(a) == 4 and a[3] == 'NGON':
            return transform.breakElements(a)
        else: return [a]

def dual(array, extraPoints=1):
    """Returns the dual mesh of a conformal mesh.
    Usage: dual(array, extraPoints)"""
    try:
        import Generator as G
        a = Converter.convertArray2NGon(array)
        a = G.close(a)
    except: # NODE forcement
        return array
    
    if isinstance(array[0], list):
        out = []
        for i in array:
            out.append(transform.dualNGon(i, extraPoints))
        return out
    else: return transform.dualNGon(a, extraPoints)
    
def splitSharpEdges(array, alphaRef=30.):
    """Split array into smooth zones (angles between elements are less than
    alphaRef).
    Usage: splitSharpEdges(array, alphaRef)"""
    try: import Generator as G
    except: raise ImportError("splitSharpEdges: requires Converter, Generator modules.")
    if isinstance(array[0], list):
        out = []
        for i in array:
            if len(i) == 5: # structured
                out.append(Converter.convertArray2Hexa(i))
            else: out.append(i)
        out = join(out)
    else:
        if len(array) == 5: out = Converter.convertArray2Hexa(array)
        else: out = array
    if out[3] == "TRI" or out[3] == "QUAD":
        out = G.close(out)
        out = reorder(out, (1,))
    return transform.splitSharpEdges(out, alphaRef)

#-----------------------------------------------------------------------------
# Decoupage d'une courbe 1D en fonction de l'angle de courbure
#-----------------------------------------------------------------------------
def splitCurvatureAngle(array, sensibility):
    """Split a line following curvature angle.
    Usage: splitCurvatureAngle(array, sensibility)"""
    out = []; array3 = array; ispl = 1 
    if len(array) != 5:
        raise TypeError("splitCurvatureAngle: defined for a i-array only.")

    while ispl > 0:
        im = array3[2]; jm = array3[3]; km = array3[4]
        ispl = transform.splitCurvatureAngle( array3, sensibility )
        if (ispl == 0 or ispl == im):
            try:
                f = join(array3, out[0])
                ispl = transform.splitCurvatureAngle( f, sensibility )
                if (ispl == 0 or ispl == f[2]): out[0] = f
                else: out.append(array3)
            except:
                out.append(array3)
            break
        arrayL = subzone(array3, (1, 1, 1), (ispl, jm, km))
        out.append(arrayL)
        array3 = subzone(array3, (ispl, 1, 1), (im, jm, km))
    return out

#-----------------------------------------------------------------------------
# Split blocks at borders joined to multiple blocks
#-----------------------------------------------------------------------------
def splitMultiplePts2D__(A):
    restart = 0
    nzones = len(A)
    for noz in range(nzones):
        z = A[noz]
        taga = Converter.extractVars(z,['definedBC'])
        ni = taga[2]; nj = taga[3]; nk = taga[4]; ninj = ni*nj
        isplit = -1; jsplit = -1
        # detecte si un pt interieur est de tag > 1
        # fenetre i = 1
        for j in range(1,nj-1): 
            if taga[1][0,j*ni] > 1.:
                isplit = 1; jsplit = j+1
                z1 = subzone(z,(1,1,1),(ni,jsplit,nk))
                z2 = subzone(z,(1,jsplit,1),(ni,nj,nk))
                del A[noz]; A+= [z1,z2]
                restart = 1
                return A, restart
             
        # fenetre i = ni
        for j in range(1,nj-1): 
            if taga[1][0,ni-1+j*ni] > 1.:
                isplit = ni; jsplit = j+1
                z1 = subzone(z,(1,1,1),(ni,jsplit,nk))
                z2 = subzone(z,(1,jsplit,1),(ni,nj,nk))
                del A[noz]; A+= [z1,z2]
                restart = 1
                return A, restart
                
        # fenetre j = 1
        for i in range(1,ni-1): 
            if taga[1][0,i] > 1.:
                isplit = i+1; jsplit = 1
                z1 = subzone(z,(1,1,1),(isplit,nj,nk))
                z2 = subzone(z,(isplit,1,1),(ni,nj,nk))
                del A[noz]; A+= [z1,z2]
                restart = 1
                return A, restart
            
        # fenetre j = nj
        for i in range(1,ni-1): 
            if taga[1][0,i+(nj-1)*ni] > 1.:
                isplit = i+1; jsplit = nj
                z1 = subzone(z,(1,1,1),(isplit,nj,nk))
                z2 = subzone(z,(isplit,1,1),(ni,nj,nk))
                del A[noz]; A+= [z1,z2]
                restart = 1
                return A, restart   
    return A, restart

def splitMultiplePts3D__(A):
    restart = 0
    nzones = len(A)
    for noz in range(nzones):
        z = A[noz]
        taga = Converter.extractVars(z,['definedBC'])
        ni = taga[2]; nj = taga[3]; nk = taga[4]; ninj = ni*nj
        isplit = -1; jsplit = -1; ksplit = -1
        # detecte si un pt interieur est de tag > 1
        ni1 = max(2,ni-1); nj1 = max(2,nj-1); nk1 = max(2,nk-1)
        # fenetre i = 1
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
                                del A[noz]; A+=[z1,z2]; restart = 1
                                return A,restart

                        # split en j ?
                        indm = ind-ninj; indp = ind+ninj
                        if indp > ninj*nk-1: indp = ind
                        if taga[1][0,indm] > 1. and taga[1][0,indp] > 1.:
                            z1 = subzone(z,(1,1,1),(ni,jsplit,nk))
                            z2 = subzone(z,(1,jsplit,1),(ni,nj,nk))
                            del A[noz]; A+=[z1,z2]; restart = 1
                            return A,restart
     
        # fenetre j = 1
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
                                del A[noz]; A+=[z1,z2]; restart = 1
                                return A,restart

                        # split en i ?
                        indm = ind-ninj; indp = ind+ninj;
                        if indp > ninj*nk-1: indp = ind
                        if taga[1][0,indm] > 1. and taga[1][0,indp] > 1.:
                            z1 = subzone(z,(1,1,1),(isplit,nj,nk))
                            z2 = subzone(z,(isplit,1,1),(ni,nj,nk))
                            del A[noz]; A+=[z1,z2]; restart = 1
                            return A,restart
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
                                del A[noz]; A+=[z1,z2]; restart = 1
                                return A,restart
                            # split en j ?
                            indm = ind-1; indp = ind+1
                            if taga[1][0,indm] > 1. and taga[1][0,indp] > 1.:
                                z1 = subzone(z,(1,1,1),(ni,jsplit,nk))
                                z2 = subzone(z,(1,jsplit,1),(ni,nj,nk))
                                del A[noz]; A+=[z1,z2]; restart = 1
                                return A,restart               
    return A, restart

def splitMultiplePts__(A,dim=3):
    try: import Generator as G
    except: raise ImportError("splitMultiplePts requires Converter and Generator modules.")
    nzones = len(A)
    allWins =[]
    tags = Converter.addVars(A, 'definedBC')
    tags = Converter.extractVars(tags, ['definedBC'])

    for noz1 in range(nzones):
        z = A[noz1]; ni=z[2]; nj=z[3]; nk=z[4] 
        winp=subzone(z,(1,1,1),(1,nj,nk));allWins.append(winp)     
        winp=subzone(z,(ni,1,1),(ni,nj,nk));allWins.append(winp)             
        winp=subzone(z,(1,1,1),(ni,1,nk));allWins.append(winp)             
        winp=subzone(z,(1,nj,1),(ni,nj,nk));allWins.append(winp)             
        if dim == 3:
            winp=subzone(z,(1,1,1),(ni,nj,1));allWins.append(winp)             
            winp=subzone(z,(1,1,nk),(ni,nj,nk));allWins.append(winp)             
    globWin = Converter.convertArray2Hexa(allWins); globWin = join(globWin); globWin = G.close(globWin)
    hook = Converter.createHook(globWin,function='nodes')
    tagG = [-1]*globWin[1].shape[1]
    for noz1 in range(nzones):
        res = Converter.identifyNodes(hook,A[noz1])
        for ind in range(A[noz1][1].shape[1]):
            if res[ind] != -1: 
                indg = res[ind]-1; tagG[indg]+=1
                    
    for noz1 in range(nzones):
        tag1 = tags[noz1]
        res = Converter.identifyNodes(hook,A[noz1])
        for ind in range(A[noz1][1].shape[1]):
            if res[ind] != -1: 
                indg = res[ind]-1
                tag1[1][0,ind] = tagG[indg]
    A = Converter.addVars([A,tags])
    Converter.freeHook(hook)

    # split des zones si definedBC est superieur a 1 en un pt (2D) ou sur une ligne i,j ou k (3D)
    split = 1; count = 0
    while split == 1:
        count += 1
        if dim == 2: A, split = splitMultiplePts2D__(A)
        else: A, split = splitMultiplePts3D__(A)
    A = Converter.rmVars(A, 'definedBC')
    return A, count

def splitMultiplePts(A, dim=3):
    """Split any zone of A if it is connected to several blocks at a given border.
    Usage: splitMultiplePts(A, dim)"""
    count = 2
    while count > 1:
        A,count = splitMultiplePts__(A, dim)
    return A

def splitBAR(array, N):
    """Split BAR at index N (start 0).
    Usage: splitBAR(array, N)"""
    a = transform.splitBAR(array, N)
    A = splitConnexity(a)
    return A

def splitTBranches(array, tol=1.e-13):
    """Split a BAR into a set of BARS at vertices where T branches exist.
    Usage: splitTBranches(array, tol=1.e-13)"""
    if isinstance(array[0], list): 
        b = []
        for i in array:
            res = transform.splitTBranches(i, tol)
            if not isinstance(res[0], list): res = [res]
            b += res
        return b
    else:
        res = transform.splitTBranches(array, tol)
        if not isinstance(res[0], list): res = [res] 
        return res
    
def splitTRI(array, idxList):
    """Split a TRI into several TRIs delimited by the input poly line defined by the lists of indices idxList.
    Usage: splitTRI(array, idxList)"""
    return transform.splitTRI(array, idxList)
    
def splitManifold(array):
    """Split an unstructured mesh (only TRI or BAR currently) into several manifold pieces.
    Usage: splitManifold(array)"""
    return transform.splitManifold(array)
