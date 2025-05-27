"""Post-processing of solutions."""
__version__ = '4.1'
__author__ = "Stephanie Peron, Christophe Benoit, Gaelle Jeanfaivre, Pascal Raud, Christelle Wervaecke, Xavier Juvigny"

from . import post
import numpy
try: import Converter
except ImportError: raise ImportError("Post: requires Converter module.")

try: range = xrange
except: pass

## [AJ - KEEP FOR NOW - FROM MASTER]
__all__ = ['coarsen', 'computeCurl', 'computeDiff', 'computeExtraVariable',
           'computeGrad', 'computeGrad2', 'computeGradLSQ',
           'computeDiv', 'computeDiv2', 'computeIndicatorField',
           'computeIndicatorFieldForBounds', 'computeIndicatorValue',
           'computeNormCurl', 'computeNormGrad', 'computeVariables',
           'computeVariables2', '_computeVariables2',
           'enforceIndicatorForCoarsestLevel', 'enforceIndicatorForFinestLevel',
           'enforceIndicatorNearBodies', 'exteriorElts', 'exteriorEltsStructured',
           'exteriorFaces', 'exteriorFacesStructured', 'extractMesh', 'extractPlane',
           'extractPoint', 'frontFaces', 'integ', 'integMoment', 'integMomentNorm',
           'integNorm', 'integNormProduct', 'interiorFaces', 'isoLine', 'isoSurf',
           'isoSurfMC', 'perlinNoise', 'projectCloudSolution',
           'refine', 'renameVars', 'selectCells', 'selectCells2', 'selectCells3',
           'sharpEdges', 'silhouette', 'slice', 'streamLine', 'streamLine2',
           'streamRibbon', 'streamRibbon2', 'streamSurf', 'usurp', 'zip', 'zipper',
           'growOfEps__','computeIndicatorField_AMR']

#==============================================================================
# Add two layers to surface arrays
# if planarity=True: height is max(eps, cellPlanarity)
# else height is eps
# retourne 1 si tous les arrays ont ete etendus
#==============================================================================
def extrudeLayer__(i, nlayers, planarity, eps, dplus, dmoins):
    import Generator; import Transform
    if planarity:
        p = Generator.getCellPlanarity(i)
        epsmax = max(eps, 2*Converter.getMaxValue(p, 'dist'))
    else: epsmax = eps
    if i[3] == 'BAR' or (i[3] == 1 and i[4] == 1): # 1D
        for k in range(nlayers+1): dplus[1][0,k] = k*epsmax; dmoins[1][0,k] =-k*epsmax
        b = Generator.addNormalLayers(i, dplus)
        c = Generator.addNormalLayers(i, dmoins)
        b = Converter.convertArray2Tetra(b)
        c = Converter.convertArray2Tetra(c)
        p = Transform.join(b, c); p = Generator.close(p)
    else: # other
        j = Converter.convertArray2Tetra(i)
        for k in range(nlayers+1): dplus[1][0,k] = k*epsmax; dmoins[1][0,k] =-k*epsmax
        j = Transform.reorder(j, (1,))
        b = Generator.addNormalLayers(j, dplus)
        j = Transform.reorder(j, (-1,))
        c = Generator.addNormalLayers(j, dplus)
        p = Transform.join(b, c); p = Generator.close(p)
        p = Converter.convertArray2Tetra(p)

    if p[3] == 'TRI': # une BAR au depart
        p = Transform.reorder(p, (1,))
        b = Generator.addNormalLayers(p, dplus)
        p = Transform.reorder(p, (-1,))
        c = Generator.addNormalLayers(p, dplus)
        p = Transform.join(b, c); p = Generator.close(p)
        p = Converter.convertArray2Tetra(p)
    return p

def growOfEps__(arrays, eps, nlayers=1, planarity=True):
    inl = []
    dplus = Converter.array('h', nlayers+1, 1, 1)
    dmoins = Converter.array('h', nlayers+1, 1, 1)
    modified = 0
    for i in arrays:
        if len(i) == 5: # structure
            if i[2] == 1 or i[3] == 1 or i[4] == 1: # 2D : to be extruded
                p = extrudeLayer__(i, nlayers, planarity, eps, dplus, dmoins)
                modified += 1
                inl.append(p)
            else: # 3D OK for ADT
                inl.append(i)
        elif len(i) == 4: # non-structure
            if i[3] == 'TRI' or i[3] == 'QUAD':
                p = extrudeLayer__(i, nlayers, planarity, eps, dplus, dmoins)
                modified += 1
                inl.append(p)
            else:
                i = Converter.convertArray2Tetra(i)
                if i[3] == 'TRI': # if TRI : extrude to be ok for ADT
                    p = extrudeLayer__(i, nlayers, planarity, eps, dplus, dmoins)
                    modified += 1
                    inl.append(p)
                else: # NGON is now a TETRA : ok for ADT
                    inl.append(i)
    if modified == len(arrays): modified = 1
    else: modified = 0
    #Converter.convertArrays2File(inl, 'ext.plt')
    return inl, modified

def extractPoint(arrays, Pts, order=2, extrapOrder=1,
                 constraint=40., tol=1.e-6, hook=None):
    """Extract the solution in one or more points.
    Usage: extractPoint( arrays, Pts, order,extrapOrder,constraint,hook)"""
    if arrays[0][1].shape[0] < 4: return [] # nofield
    inl, modified = growOfEps__(arrays, tol, nlayers=2, planarity=False)
    if isinstance(Pts,list): res = post.extractPoint(inl, Pts, order, extrapOrder, constraint, hook)
    else: res = post.extractPoint(inl, [Pts], order, extrapOrder, constraint, hook)
    out = []
    npts = res[1].shape[1]; nfld = res[1].shape[0]
    for i in range(npts):
        outi = []
        for nof in range(nfld): outi.append(res[1][nof,i])
        out.append(outi)
    if isinstance(Pts, list): return out
    else: return out[0]

def extractPlane(arrays, T, order=2, tol=1.e-6):
    """Slice solution with a plane.
    Usage: extractPlane(arrays, (coefa, coefb, coefc, coefd), order)"""
    try: import Generator; import Transform
    except ImportError:
        return post.extractPlane(arrays,T, order)
    inl, modified = growOfEps__(arrays, tol, nlayers=2, planarity=False)
    ret1 = post.extractPlane(inl, T, order)
    if modified == 0: return ret1
    else:
        ret1 = Generator.close(ret1)
        ret1 = Transform.collapse(ret1)
        return ret1

def extractMesh(arrays, extractArray, order=2, extrapOrder=1,
                constraint=40., tol=1.e-6, hook=None):
    """Extract the solution on a given mesh.
    Usage: extractMesh(arrays, extractArray, order, hook)"""
    inl, modified = growOfEps__(arrays, tol, nlayers=2, planarity=False)
    if isinstance(extractArray[0], list):
        return post.extractMesh(inl, extractArray, order, extrapOrder,
                                constraint, hook)
    else:
        return post.extractMesh(inl, [extractArray], order, extrapOrder,
                                constraint, hook)[0]

def slice(a, type=None, eq=None):
    """Extract a slice of different shapes."""
    if eq is not None:
        eq = eq.replace(' ', '')
        eq = eq.replace('=0.', '')
        eq = eq.replace('=0', '')
        eq = '{eslice}='+eq

    if type == 'cone':
        if eq is None: raise ValueError('slice: equation is needed.')
        a = Converter.initVars(a, '{r}=sqrt({z}*{z}+{y}*{y})')
        a = Converter.initVars(a, eq)
        p = isoSurfMC(a, 'eslice', 0.)
        return p
    elif type == 'cone_struct':
        raise ValueError('slice: Structured slices not already implemented.')
    elif type == 'plane':
        a = Converter.initVars(a, eq)
        p = isoSurfMC(a, 'eslice', 0.)
        return p
    return a

def projectCloudSolution(cloudArray, surfArray, dim=3, loc='nodes', ibm=False, old=False):
    """Project the solution defined on a set of points to a TRI surface."""
    surfArray = Converter.convertArray2Tetra(surfArray)
    if isinstance(surfArray[0], list):
        try:
            import Transform
            surfArray = Transform.join(surfArray)
        except: pass
    return projectCloudSolution__(cloudArray,surfArray,dim=dim,loc=loc,ibm=ibm, old=old)

# Much lighter than projectCloudSolution: no conversion to TETRA and no join.
def projectCloudSolution__(cloudArray, surfArray, dim=3, loc='nodes', ibm=False, old=False):
    """Project the solution defined on a set of points to a TRI surface."""
    cloudArray = Converter.convertArray2Node(cloudArray)
    if loc == 'centers': surfArray = Converter.node2Center(surfArray)
    return post.projectCloudSolution2Triangle(cloudArray,surfArray,dim,int(ibm), int(old))

def projectCloudSolutionWithInterpData(cloudArray, surfArray, offset, interpDonor, interpCoef, dim=3, loc='nodes'):
    """Project the solution defined on a set of points to a TRI surface using pre-calculated interpolation data."""
    surfArray = Converter.convertArray2Tetra(surfArray)
    if isinstance(surfArray[0], list):
        try:
            import Transform
            surfArray = Transform.join(surfArray)
        except: pass
    return projectCloudSolutionWithInterpData__(cloudArray, surfArray, offset, interpDonor, interpCoef, dim=dim, loc=loc)

# Much lighter than projectCloudSolution: no conversion to TETRA and no join.
def projectCloudSolutionWithInterpData__(cloudArray, surfArray, offset, interpDonor, interpCoef, dim=3, loc='nodes'):
    """Project the solution defined on a set of points to a TRI surface using pre-calculated interpolation data."""
    cloudArray = Converter.convertArray2Node(cloudArray)
    if loc == 'centers': surfArray = Converter.node2Center(surfArray)
    return post.projectCloudSolution2TriangleWithInterpData(cloudArray, surfArray, offset, interpDonor, interpCoef, dim)

def prepareProjectCloudSolution(cloudArray, surfArray, dim=3, loc='nodes', ibm=False):
    """Compute the MLS interpolation data for projectCloudSolutionWithInterpData."""
    surfArray = Converter.convertArray2Tetra(surfArray)
    if isinstance(surfArray[0], list):
        try:
            import Transform
            surfArray = Transform.join(surfArray)
        except: pass
    return prepareProjectCloudSolution__(cloudArray,surfArray,dim=dim,loc=loc,ibm=ibm)

# Much lighter than projectCloudSolution: no conversion to TETRA and no join.
def prepareProjectCloudSolution__(cloudArray, surfArray, dim=3, loc='nodes', ibm=False):
    """Compute the MLS interpolation data for projectCloudSolutionWithInterpData."""
    cloudArray = Converter.convertArray2Node(cloudArray)
    if loc == 'centers': surfArray = Converter.node2Center(surfArray)
    return post.prepareProjectCloudSolution2Triangle(cloudArray,surfArray,dim,int(ibm))

def coarsen(a, indic, argqual=0.1, tol=1.e6):
    """Coarsen a surface TRI-type mesh given a coarsening indicator for each
    element.
    Usage: coarsen(a, indic, argqual, tol)"""
    if isinstance(a[0], list):
        b = []
        for i in range(len(a)):
            b.append(post.coarsen(a[i], indic[i], argqual, tol))
        return b
    else:
        return post.coarsen(a, indic, argqual, tol)

def refine(a, indic=None, w=-1):
    """Refine a surface TRI-type mesh given an indicator for each
    element.
    Usage: refine(a, indic, w)"""
    if isinstance(a[0], list):
        b = []
        for i in range(len(a)):
            b.append(refine__(a[i], indic[i], w))
        return b
    else:
        return refine__(a, indic, w)

def refine__(a, indic, w):
    if indic is not None: return post.refine(a, indic) # linear
    else: return post.refineButterfly(a, w) # butterfly

def interiorFaces(a, strict=0):
    """Interior faces of an array a. The argument strict equal to 1 means
    that interior faces with only interior nodes are taken into account.
    Usage: interiorFaces(array, strict)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(post.interiorFaces(i, strict))
        return b
    else:
        return post.interiorFaces(a, strict)

#==============================================================================
# Return a list of structured arrays defining exterior faces for a structured
# array
#==============================================================================
def exteriorFacesStructured(a):
    """Return the list of exterior faces for a structured array
    Usage: exteriorFacesStructured(a,indices)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(post.exteriorFacesStructured(i))
        return b
    else:
        return post.exteriorFacesStructured(a)

#==============================================================================
# Return an unstructured array of exterior faces
#==============================================================================
def exteriorFaces(a, indices=None):
    """Exterior faces of an array.
    Usage: exteriorFaces(a, indices)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(exteriorFacesForOneArray__(i, indices))
        return b
    else:
        return exteriorFacesForOneArray__(a, indices)

def exteriorFacesForOneArray__(a, indices):
    # To be commented in next release
    if len(a) == 4 and (a[3] == 'PENTA' or a[3] == 'PYRA'):
        try:
            import Generator
            a = Converter.convertArray2NGon(a)
            a = Generator.close(a)
        except: pass
    return post.exteriorFaces(a, indices)

def exteriorElts(array):
    """Exterior elements of an array.
    Usage: exteriorElts(a)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            if len(i) == 5: i = Converter.convertArray2Hexa(i)
            b.append(post.exteriorElts(i))
        return b
    else:
        a = array
        if len(a) == 5: a = Converter.convertArray2Hexa(a)
        return post.exteriorElts(a)

def exteriorEltsStructured(array, depth=1):
    """Exterior elements of an array as a structured grid.
    Usage: exteriorEltsStructured(a, depth)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(post.exteriorEltsStructured(i, depth))
        return b
    else:
        return post.exteriorEltsStructured(array, depth)

def integ(coordArrays, FArrays, ratioArrays):
    """Integral of fields.
    Usage: integ(coordArrays, FArrays, ratioArrays)"""
    return post.integ(coordArrays, FArrays, ratioArrays)

def integNorm(coordArrays, FArrays, ratioArrays):
    """Integral of fields times normal.
    Usage: integNorm(coordArrays, FArrays, ratioArrays)"""
    return post.integNorm(coordArrays, FArrays, ratioArrays)

def integNormProduct(coordArrays, FArrays, ratioArrays):
    """Integral of scalar product fields times normal.
    Usage: integNormProduct( coordArrays, FArrays, ratioArrays )"""
    return post.integNormProduct(coordArrays, FArrays, ratioArrays)

def integMoment(coordArrays, FArrays, ratioArrays, center):
    """Integral of moments.
    Usage: integMoment( coordArrays, FArrays, ratioArrays, (xc, yc, zc) )"""
    return post.integMoment(coordArrays, FArrays, ratioArrays, center)

def integMomentNorm(coordArrays, FArrays, ratioArrays, center):
    """Integral of moments (OM^f.vect(n)).
    Usage:
    integMomentNorm( coordArrays, FArrays, ratioArrays, (xc, yc, zc) )"""
    return post.integMomentNorm(coordArrays, FArrays, ratioArrays, center)

def zipper(arrays, options=[]):
    """Extract Chimera surface as an unique unstructured surface.
    Usage: zipper(arrays, options)"""
    return post.zipper(arrays, options)

def zip(arrays,  options=[]):
    """Extract Chimera surface as an unique unstructured surface.
    Usage: zip(arrays, options)"""
    return post.zip(arrays, options)

def usurp(blkArrays, ibArrays=[]):
    """Extract unique surfaces using ranked polygons.
    Usage: usurp(blkArrays, ibArrays)"""
    return post.usurp(blkArrays, ibArrays)

def computeVariables(array, varname,
                     gamma=1.4, rgp=287.053, s0=0., betas=1.458e-6,
                     Cs=110.4, mus=1.76e-5, Ts=273.15):
    """Compute the variables defined in varname for array.
    Usage: computeVariables(array, varname, gamma=1.4, rgp=287.053, s0=0., betas=1.458e-6, Cs=110.4, mus=0., Ts=0.)"""
    if varname == []:
        #print('Warning: computeVariables: varname list is empty.')
        return array
    if isinstance(varname, str): varname = [varname]
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(post.computeVariables(i, varname, gamma, rgp, s0, betas, Cs, mus, Ts))
        return b
    else:
        return post.computeVariables(array, varname, gamma, rgp, s0, betas, Cs, mus, Ts)


def computeVariables2(array, varname,
                      gamma=1.4, rgp=287.053, s0=0., betas=1.458e-6,
                      Cs=110.4, mus=1.76e-5, Ts=273.15):
    """In place compute variable2"""
    b = Converter.copy(array)
    _computeVariables2(b, varname, gamma, rgp, s0, betas, Cs, mus, Ts)
    return b

def _computeVariables2(array, varname,
                       gamma=1.4, rgp=287.053, s0=0., betas=1.458e-6,
                       Cs=110.4, mus=1.76e-5, Ts=274.05):
    if varname == []:
        #print('Warning: computeVariables: varname list is empty.')
        return array
    if isinstance(varname, str): varname = [varname]
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(post.computeVariables2(i, varname, gamma, rgp, s0, betas, Cs, mus, Ts))
        return b
    else:
        return post.computeVariables2(array, varname, gamma, rgp, s0, betas, Cs, mus, Ts)

def computeExtraVariable(array, varname, gamma=1.4, rgp=287.53,
                         Cs=110.4, mus=1.76e-5, Ts=273.15):
    """Compute variables that require a change of location."""
    from . import extraVariables
    if varname == 'Vorticity':
        return extraVariables.computeVorticity(array)
    elif varname == 'VorticityMagnitude':
        return extraVariables.computeVorticityMagnitude(array)
    elif varname == 'QCriterion':
        return extraVariables.computeQCriterion(array)
    elif varname == 'ShearStress':
        return extraVariables.computeShearStress(array, gamma, rgp, Cs,
                                                 mus, Ts)
    elif varname == 'SkinFriction':
        return extraVariables.computeSkinFriction(array, tangent=0)
    elif varname == 'SkinFrictionTangential':
        return extraVariables.computeSkinFriction(array, tangent=1)
    else:
        print('Warning: computeExtraVariable: unknown variable: %s.'%varname)

def perlinNoise(array, alpha=2., beta=2., n=8):
    """Compute perlin noise for array.
    Usage: perlinNoise(array, alpha, beta, n)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(post.perlinNoise(i, alpha, beta, n))
        return b
    else:
        return post.perlinNoise(array, alpha, beta, n)

def computeGrad(array, varname):
    """Compute the gradient of the field varname defined in array.
    Usage: computeGrad(array, varname) """
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(post.computeGrad(i, varname))
        return b
    else:
        return post.computeGrad(array, varname)

def computeGrad2(array, arrayc, vol=None, cellN=None, indices=None, BCField=None):
    """Compute the gradient of a field defined on centers."""
    if isinstance(array[0], list):
        raise ValueError("computeGrad2: input must be a single zone.")
    if len(array) == 4:
        if array[3] == 'NGON' and arrayc[3] == 'NGON*':
            return post.computeGrad2NGon(array, arrayc, vol, cellN, indices, BCField)
        else:
            raise ValueError("computeGrad2: only valid for NGon unstructured zones.")
    else:
        return post.computeGrad2Struct(array, arrayc, cellN, indices, BCField)

def computeNormGrad(array, varname):
    """Compute the norm of gradient of field varname defined in array.
    Usage: computeNormGrad(array, varname) """
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(post.computeNormGrad(i, varname))
        return b
    else:
        return post.computeNormGrad(array, varname)

def computeGradLSQ(array, arrayc):
    """Compute gradient using least-squares method."""
    return post.computeGradLSQ(array, arrayc)

def computeDiv(array, vector):
    """Compute the divergence of the given vector, whose components are defined in array
    using the computeGrad method for gradients.
    Usage: computeDiv(array, vector) """
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(post.computeDiv(i, vector))
        return b
    else:
        return post.computeDiv(array, vector)

def computeDiv2(array, arrayc, vol=None, cellN=None, indices=None, BCFieldX=None, BCFieldY=None, BCFieldZ=None):
    """Compute the divergence of the field varname, whose components are defined in array
    using the computeGrad2 method for gradients.
    Usage: computeDiv2(array, arrayc, indices, BCFieldX, BCFieldY, BCFieldZ) """
    if isinstance(array[0], list):
        raise ValueError("computeDiv2: input must be a single zone.")
    if len(array) == 4:
        if array[3] == 'NGON' and arrayc[3] == 'NGON*':
            return post.computeDiv2NGon(array, arrayc, vol, cellN, indices, BCFieldX, BCFieldY, BCFieldZ)
        else:
            raise ValueError("computeDiv2: only valid for NGon unstructured zones.")
    else:
        return post.computeDiv2Struct(array, arrayc, cellN, indices, BCFieldX, BCFieldY, BCFieldZ)

def computeCurl(array, vector):
    """Compute the curl of the 3D-field defined in array.
    vector defines the name of the 3 components of the vector on which the curl is applied.
    Usage: computeCurl(array, vector) """
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(post.computeCurl(i, vector))
        return b
    else:
        return post.computeCurl(array, vector)

def computeNormCurl(array, vector):
    """Compute the norm of the curl of the 3D-field defined in array.
    vector defines the name of the 3 components of the vector on which the curl is applied.
    Usage: computeNormCurl(array, vector) """
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(post.computeNormCurl(i, vector))
        return b
    else:
        return post.computeNormCurl(array, vector)

def computeDiff(array, varname):
    """Compute the difference of the field varname defined in array.
    Usage: computeDiff(array, varname) """
    import KCore
    if isinstance(array[0], list):
        A = []
        for a in array:
            b =  Converter.extractVars(a,[varname])
            posc = KCore.isNamePresent(a, 'cellN')
            if posc != -1:
                celln = Converter.extractVars(a, ['cellN'])
                b = Converter.addVars([b,celln])
            A.append(b)

        b = []
        for i in A:
            b0 = post.computeDiff(i, varname); b0[0] = 'diff'+b0[0]
            b.append(b0)
        return b
    else:
        b = Converter.extractVars(array,[varname])
        posc = KCore.isNamePresent(array, 'cellN')
        if posc != -1:
            celln = Converter.extractVars(array, ['cellN'] )
            b = Converter.addVars([b,celln])
        b = post.computeDiff(b, varname)
        b[0] = 'diff'+b[0]
        return b

def streamLine(arrays, X0, vector, N=2000, dir=2):
    """Compute a streamline starting from (x0,y0,z0) given
    a list of arrays containing 'vector' information.
    Usage: streamLine(arrays, (x0,y0,z0), vector, N, dir)"""
    # get an (unstructured) array containing all  2D-surface arrays
    surf = []
    for a in arrays:
        elt='None'
        ni = 2; nj=2; nk=2
        if len(a) == 5: # structure
            ni = a[2]; nj = a[3]; nk = a[4]
        else: elt = a[3]
        mult = (ni - 1)*(nj - 1)*(nk - 1)
        add = (ni - 1)*(nj - 1) + (ni - 1)*(nk - 1) + (nj - 1)*(nk - 1)
        if ((mult == 0) and (add != 0)) or (elt == 'QUAD') or (elt == 'TRI'):
            a = Converter.convertArray2Tetra(a)
            surf.append(a)
    if surf != []:
        try: import Transform, Generator
        except: raise ImportError("streamLine: requires Transform and Generator modules.")
        surf = Transform.join(surf)
        surf = Generator.close(surf)
    # grow 2D arrays
    if surf != []:
        tol = 1.
        inl, modified = growOfEps__(arrays, tol, nlayers=2, planarity=False)
        arrays = inl

    if dir == +1:
        return post.compStreamLine(arrays, surf, X0, vector, 1., N)
    elif dir == -1:
        return post.compStreamLine(arrays, surf, X0, vector, -1., N)
    else:
        try: a = post.compStreamLine(arrays, surf, X0, vector, 1., N)
        except: a = 0
        try: b = post.compStreamLine(arrays, surf, X0, vector, -1., N)
        except: b = 0
        if a != 0 and b != 0:
            try: import Transform
            except: return a
            b = Transform.reorder(b, (-1,2,3))
            c = Transform.join(b, a)
            return c
        elif b == 0 and a != 0: return a
        elif a == 0 and b != 0: return b
        else: raise ValueError('Empty streamline.')

# IN: arrays: coords + solution
# IN: X0: (x,y,z) ou liste de tuples
# IN: vector: ['vx','vy','vz']
# IN: eps: hauteur extrusion quand on nous donne une surface
def streamLine2(arrays, X0, vector, N=2000, dir=2, eps=1.e-2):
    """Compute a streamline starting from (x0,y0,z0) given
    a list of arrays containing 'vector' information.
    Usage: streamLine(arrays, (x0,y0,z0), vector, N, dir)"""
    # get an (unstructured) array containing all  2D-surface arrays
    surf = []
    for a in arrays:
        elt = 'None'
        ni = 2; nj = 2; nk = 2
        if len(a) == 5: # structure
            ni = a[2]; nj = a[3]; nk = a[4]
        else: elt = a[3]
        mult = (ni - 1)*(nj - 1)*(nk - 1)
        add = (ni - 1)*(nj - 1) + (ni - 1)*(nk - 1) + (nj - 1)*(nk - 1)
        if (mult == 0 and add != 0) or elt == 'QUAD' or elt == 'TRI':
            a = Converter.convertArray2Tetra(a)
            surf.append(a)
    if surf != []:
        try: import Transform, Generator
        except: raise ImportError("streamLine: requires Transform and Generator modules.")
        surf = Transform.join(surf)
        surf = Generator.close(surf)
    # grow 2D arrays
    if surf != []:
        inl, modified = growOfEps__(arrays, eps, nlayers=2, planarity=False)
        arrays = inl
    rets = post.comp_stream_line(arrays, surf, X0, vector, dir, N)
    return rets

def streamRibbon(arrays, X0, N0, vector, N=2000, dir=2):
    """Compute a streamribbon starting from (x0,y0,z0) given
    a list of arrays containing 'vector' information. The width and orientation
    of the ribbon is given by N0 the normal vector to the streamline starting
    from X0.
    Usage: streamRibbon(arrays, (x0,y0,z0), (n0x,n0y,n0z), N, dir)"""
    if dir == +1:
        return post.compStreamRibbon(arrays, X0, N0, vector, 1., N)
    elif dir == -1:
        return post.compStreamRibbon(arrays, X0, N0, vector, -1., N)
    else:
        try: a = post.compStreamRibbon(arrays, X0, N0, vector, 1., N)
        except: a = 0
        try: b = post.compStreamRibbon(arrays, X0, N0, vector, -1., N)
        except: b = 0

        if a != 0 and b != 0:
            try: import Transform
            except: return a
            c = Transform.join(b,a)
            return c
        elif b == 0 and a != 0: return a
        elif a == 0 and b != 0: return b
        else: raise

def streamRibbon2(arrays, X0, vector, N=2000, dir=1, width=1.):
    """Compute a streamribbon starting from (x0,y0,z0) given
    a list of arrays containing 'vector' information. The width could be given (default 1).
    Usage: streamRibbon2(arrays, (x0,y0,z0), vector, dir)"""
    r = post.comp_stream_ribbon(arrays, X0, vector, dir, N, width)
    return r


def streamSurf(arrays, b, vector, N=2000, dir=1):
    """Compute a stream surface."""
    b = Converter.convertArray2Hexa(b)
    if b[3] != 'BAR': raise TypeError("streamSurf: b must be a BAR.")
    coord = b[1]
    for no in range(len(arrays)):
        if len(arrays[no])!= 5:
            arrays[no] = Converter.convertArray2Tetra(arrays[no])
    return post.compStreamSurf(arrays, b, vector, dir, N)

#------------------------------------------------------------------------------
# Construit le tag pour les formules
#------------------------------------------------------------------------------
def buildTag2__(array, F):
    a = Converter.copy(array)

    # Extrait les variables de a
    varstring = a[0]
    vars = varstring.split(',')

    eq = F.replace('centers:', '')

    # Instantiation de la formule
    ap = a[1]
    loc = eq
    c = 0
    for v in vars:
        loc = loc.replace('{%s}'%v, 'ap[%d,:]'%c); c += 1

    # Evaluation
    formula = eval(loc) # est un numpy bool
    tag = numpy.zeros(ap.shape[1], numpy.float64)
    tag[:] = formula[:]
    tag = tag.reshape(1, tag.size)
    if len(array) == 5: out = ['__tag__', tag, array[2], array[3], array[4]]
    else: out = ['__tag__', tag, array[2], array[3]]
    return out

#------------------------------------------------------------------------------
# construit le tag pour une fonction F
#------------------------------------------------------------------------------
def buildTag1__(array, F, varStrings):
    import KCore
    pos = []
    for i in varStrings:
        i = i.replace('centers:', '')
        p = KCore.isNamePresent(array, i)+1
        if p == 0:
            raise ValueError("selectCells: can't find %s in array."%i)
        else:
            pos.append(p)
    n = array[1]
    nsize = n.shape[1]
    l = len(varStrings)
    tag = numpy.zeros(nsize, numpy.float64)

    if l == 0:
        if F(): tag[:] = 1
    else:
        for i in range(nsize):
            x = [n[pos[j]-1,i] for j in range(l)]
            if F(*x): tag[i] = 1
    tag = tag.reshape(1, tag.size)
    if len(array) == 5: out = ['__tag__', tag, array[2], array[3], array[4]]
    else: out = ['__tag__', tag, array[2], array[3]]
    return out


def selectCells__(arrayNodes, F, arrayCenters, varStrings, strict, F2E, cleanConnectivity):
    if varStrings == []: tag = buildTag2__(arrayNodes, F)
    else: tag = buildTag1__(arrayNodes, F, varStrings)

    if arrayCenters != []:
        return post.selectCellsBoth(arrayNodes, arrayCenters, tag, strict, F2E, cleanConnectivity)
    else:
        return post.selectCells(arrayNodes, tag, strict, F2E, cleanConnectivity)

def selectCells(arrayNodes, F, arrayCenters=[], varStrings=[], strict=0, F2E=None, cleanConnectivity=True):
    """Select cells in a given array.
    Usage: selectCells(array, F, varStrings, strict)"""
    if isinstance(arrayNodes[0], list):
        b = []
        if arrayCenters != []:
            if len(arrayNodes) != len(arrayCenters): raise ValueError("selectCells: Nodes and Centers arrays have different size.")

        for i, an in enumerate(arrayNodes):
            if arrayCenters != []:
                ret = selectCells__(an, F, arrayCenters[i], varStrings, strict, F2E)
            else:
                ret = selectCells__(an, F, [], varStrings, strict, F2E)
            if F2E is None and arrayCenters == []: b.append(ret[0])
            else: b.append(ret)
        return b
    else:
        ret = selectCells__(arrayNodes, F, arrayCenters, varStrings, strict, F2E, cleanConnectivity)

        if F2E is None and arrayCenters == []: return ret[0]
        else: return ret

#------------------------------------------------------------------------------
# Select cells where tag=1
# loc=-1: unknown, 1: centers, 0: nodes (forced)
#------------------------------------------------------------------------------
# selectCells preserving center flow field solutions
# an : coordinates and fields in nodes
# ac : fields in centers
def selectCells2(an, tag, ac=[], strict=0, loc=-1, F2E=None, cleanConnectivity=True):
    """Select cells in a given array following tag.
    Usage: selectCells2(arrayN, arrayC, tag, strict)"""
    if isinstance(an[0], list):
        b = []
        lenan = len(an)
        lenac = len(ac)
        if ac != []:
            if lenan != lenac: raise ValueError("selectCells2: Nodes and Centers arrays have different size.")
        for i in range(lenan):
            sizetag = tag[i][1].shape[1]
            sizean  = an[i][1].shape[1]
            if sizetag != sizean or loc == 1: # centers
                if F2E is not None:
                    (PE2, retn, retc) = post.selectCellCenters(an[i], ac[i], tag[i], F2E, cleanConnectivity)
                else:
                    (retn, retc) = post.selectCellCenters(an[i], ac[i], tag[i], cleanConnectivity)
            else:
                if ac == []:
                    if F2E is not None:
                        (PE2, retn) = post.selectCells(an[i], tag[i], strict, F2E, cleanConnectivity)
                    else:
                        retn = post.selectCells(an[i], tag[i], strict, None, cleanConnectivity)[0]
                else:
                    if F2E is not None:
                        (PE2, retn, retc) = post.selectCellsBoth(an[i], ac[i], tag[i], strict, F2E, cleanConnectivity)
                    else:
                        (retn, retc) = post.selectCellsBoth(an[i], ac[i], tag[i], strict, None, cleanConnectivity)

            if ac == []:
                if F2E is None: b.append(retn)
                else: b.append(PE2, retn)
            else:
                if F2E is None: b.append((retn, retc))
                else: b.append((PE2,retn, retc))

        return b

    else:
        sizetag = tag[1].shape[1]
        sizean  = an[1].shape[1]
        if sizean != sizetag or loc == 1: # centers
            if ac == []:
                if F2E is not None:
                    (PE2, retn) = post.selectCellCenters(an, tag, F2E, cleanConnectivity)
                else:
                    retn = post.selectCellCenters(an, tag, None, cleanConnectivity)[0]
            else:
                if F2E is not None:
                    (PE2, retn, retc) = post.selectCellCentersBoth(an, ac, tag, F2E, cleanConnectivity)
                else:
                    (retn, retc) = post.selectCellCentersBoth(an, ac, tag, None, cleanConnectivity)
        else:
            if ac == []:
                if F2E is not None:
                    (PE2, retn)  = post.selectCells(an, tag, strict, F2E, cleanConnectivity)
                else:
                    retn = post.selectCells(an, tag, strict, None, cleanConnectivity)[0]
            else:
                if F2E is not None:
                    (PE2, retn, retc) = post.selectCellsBoth(an, ac, tag, strict, F2E, cleanConnectivity)
                else:
                    (retn, retc) = post.selectCellsBoth(an, ac, tag, strict, None, cleanConnectivity)

        if ac != []:
            if F2E is None: return (retn,retc)
            else: return (PE2, retn, retc)
        else:
            if F2E is None: return retn
            else: return (PE2, retn)

#==============================================================================
def selectCells3(a, tag):
    try: import Transform as T
    except: raise ImportError('selectCells: Transform module is required.')
    if isinstance(a[0], list):
        b = []
        lena = len(a)
        for i in tag:
            ret = post.selectCells3(i, 0)
            ret = T.subzone(a, ret, type='elements')
            b.append(ret)
        return b
    else:
        ret = post.selectCells3(tag, 0)
        ret = T.subzone(a, ret, type='elements')
        return ret

#==============================================================================
def frontFaces(a, tag):
    """Select faces located in the front of tag=0 and tag=1.
    Usage: frontFaces(a, tag)"""
    if isinstance(a[0], list):
        b = []
        lena = len(a)
        for i in range(lena):
            try:
                if len(a[i]) == 5: # structure
                    a[i] = Converter.convertArray2Hexa(a[i])
            except: pass
            b.append(post.frontFaces(a[i], tag[i]))
        return b
    else:
        try:
            if len(a) == 5: # structure
                a = Converter.convertArray2Hexa(a)
        except: pass
        return post.frontFaces(a, tag)

#==============================================================================
def isoLine(array, var, value):
    """Compute an isoLine correponding to value of field 'var' on arrays.
    Usage: isoLine(array, 'Density', 1.2)"""
    try:
        b = Converter.convertArray2Tetra(array)
        if isinstance(b[0], list):
            for i in range(len(b)):
                if b[i][3] == 'TETRA': b[i] = exteriorFaces(b[i])
        else:
            if b[3] == 'TETRA': b = exteriorFaces(b)
    except: b = array

    if isinstance(value, list): values = value
    else: values = [value]

    out = []
    for v in values:
        if isinstance(b[0], list):
            ret = []
            for i in b:
                try:
                    i = post.isoLine(i, var, v)
                    ret.append(i)
                except: pass
            try:
                import Transform
                if ret != []:
                    ret = Transform.join(ret); out.append(ret)
            except:
                if ret != []: out.append(ret[0])
        else:
            try:
                i = post.isoLine(b, var, v)
                out.append(i)
            except: pass

    if out == []: raise ValueError("isoLine: isoline is empty.")
    if isinstance(value, list) and len(value) == 1: return out[0]
    if not isinstance(value, list): return out[0]
    return out

#==============================================================================
def isoSurf(array, var, value, split='simple'):
    """Compute an isoSurf corresponding to value of field 'var' in
    volume arrays.
    Usage: isoSurf(array, 'Density', 1.2)"""
    try: import Transform
    except: return []

    if isinstance(array[0], list):
        ret = []
        for i in array:
            try:
                dim = 3
                if i[3] == 'NGON':
                    ap = i[2].ravel('K')
                    if ap[2] == 1: dim = 1
                    if ap[2] == 2: dim = 2
                if i[3] != 'NGON' or dim != 3:
                    i = Converter.convertArray2Tetra(i, split=split)
            except: pass
            try:
                if i[3] == 'TRI' or i[3] == 'QUAD' or i[3] == 'BAR':
                    i = post.isoLine(i, var, value)
                    ret.append(i)
                elif i[3] == 'NGON':
                    i = post.isoSurfNGon(i, var, value)
                    i = Transform.reorder(i, (1,))
                    ret.append(i)
                else:
                    i = post.isoSurf(i, var, value)
                    i = Transform.reorder(i, (1,))
                    ret.append(i)
            except: pass
        return ret
    else:
        try:
            dim = 3
            if array[3] == 'NGON':
                ap = array[2].ravel('K')
                if ap[2] == 1: dim = 1
                if ap[2] == 2: dim = 2
            if array[3] != 'NGON' or dim != 3: b = Converter.convertArray2Tetra(array, split=split)
            else: b = array
        except: b = array
        try:
            if b[3] == 'TRI' or b[3] == 'QUAD' or b[3] == 'BAR':
                b = post.isoLine(b, var, value)
            elif b[3] == 'NGON':
                b = post.isoSurfNGon(b, var, value)
                b = Transform.reorder(b, (1,))
            else:
                b = post.isoSurf(b, var, value)
                b = Transform.reorder(b, (1,))
            return [b]
        except: return []

#==============================================================================
def isoSurfMC(array, var, value, split='simple'):
    """Compute an isoSurf correponding to value of field 'var' in
    volume arrays.
    Usage: isoSurfMC(array, 'Density', 1.2)"""
    try: import Transform
    except: return []

    if isinstance(array[0], list):
        ret = []; pool = []
        for i in array:
            if len(i) == 5: i = Converter.convertArray2Hexa(i)
            if i[3] == 'HEXA':
                try:
                    i = post.isoSurfMC(i, var, value)
                    i = Transform.reorder(i, (1,))
                    ret.append(i)
                except: pass
            else: pool.append(i)

        if pool != []: ret2 = isoSurf(pool, var, value, split)
        else: ret2 = []
        return ret + ret2
    else:
        if len(array) == 5: array = Converter.convertArray2Hexa(array)
        if array[3] == 'HEXA':
            try:
                b = post.isoSurfMC(array, var, value)
                b = Transform.reorder(b, (1,))
                return [b]
            except: return []
        else: return isoSurf(array, var, value, split)

#==============================================================================
def enforceIndicatorNearBodies(indicator, octreeHexa, bodies):
    """Enforce the refinement indicator to -1000 near bodies."""
    return post.enforceIndicatorNearBodies(indicator, octreeHexa, bodies)

def enforceIndicatorForFinestLevel(indicator, octreeHexa):
    """Enforce the indicator field for the octree octreeHexa for the finest
    level: the indicator is set to -2000 for the elements of finest level."""
    return post.enforceIndicatorForFinestLevel(indicator, octreeHexa)

def enforceIndicatorForCoarsestLevel(indicator, octreeHexa):
    """Enforce the indicator field for the octree octreeHexa for the coarsest
    level: the indicator is set to -3000 for the elements of coarsest level."""
    return post.enforceIndicatorForCoarsestLevel(indicator, octreeHexa)

def computeIndicatorFieldForBounds(indicator, indicatorValues, valMin, valMax, isAMR=False):
    """Return the modified indicator field in order to obtain a number of
    points controlled by bounds epsMin and epsMax. Indicator values greater
    than epsMax are refined (indicator=1), those lower than epsMin are
    coarsened (indicator=-1), and 0 between bounds."""
    valt = indicatorValues[1]; nelts = indicator[1][0].shape[0]
    indicator2 = Converter.copy(indicator)
    indict = indicator2[1]
    for i in range(nelts):
        if indict[0,i] == -1000.: # enforce near bodies
            if isAMR:
                indict[0,i] = 0.
            else:
                if valt[0,i] >= valMax: indict[0,i] = 1.
                else: indict[0,i] = 0.
        elif indict[0,i] == -2000.: # no refinement of finest level
            indict[0,i] = 0.
            if isAMR:
                if valt[0,i] <= valMin: indict[0,i] = -1.
        elif indict[0,i] == -3000.: # coarsening of coarsest level
            indict[0,i] = -1.
            if isAMR:
                if valt[0,i] >= valMax and valt[0,i] <=1: indict[0,i] = 1.
        else:
            if isAMR:
                if valt[0,i] >= valMax and valt[0,i] <=1: indict[0,i] = 1.
                elif valt[0,i] <= valMin: indict[0,i] = -1.
                else: indict[0,i] = 0.
            else:
                if valt[0,i] >= valMax: indict[0,i] = 1.
                elif valt[0,i] <= valMin: indict[0,i] = -1.
                else: indict[0,i] = 0.
    return indicator2

def computeIndicatorValue(octreeHexa, zones, indicField):
    """Computes the indicator value on the octree mesh based on the maximum
    value of the absolute value of indicField. Returns the projected field located
    at elements of the octree mesh."""
    return post.computeIndicatorValue(octreeHexa, zones, indicField)

def computeIndicatorField(octreeHexa, indicVal, nbTargetPts=-1, bodies=[],
                          refineFinestLevel=1, coarsenCoarsestLevel=1):
    """Compute the indicator -1, 0 or 1 for each element of the HEXA octree
    with respect to the indicatorValue field located at element centers. The
    bodies fix the indicator to 0 in the vicinity of bodies. nbTargetPts
    controls the number of points after adaptation.
    If refineFinestLevel=1, the finest levels are refined.
    If coarsenCoarsestLevel=1, the coarsest levels are coarsened wherever possible.
    Returns the indicator field.
    Usage: computeIndicatorField(octreeHexa, indicVal, nbTargetPts, bodies, refinestLevel, coarsenCoarsestLevel)"""
    try: import Generator as G
    except: raise ImportError("computeIndicatorField: requires Generator module.")
    npts = octreeHexa[1][0].shape[0]
    if nbTargetPts == -1: nbTargetPts = npts

    valName = indicVal[0]; nelts = indicVal[2]
    indicVal[1] = numpy.absolute(indicVal[1])
    indicator = Converter.initVars(indicVal, 'indicator', 0.)
    indicator = Converter.extractVars(indicator, ['indicator'])

    if bodies != []:
        bodies = Converter.convertArray2Tetra(bodies)
        indicator = post.enforceIndicatorNearBodies(indicator, octreeHexa, bodies)
    if refineFinestLevel == 0:
        indicator = post.enforceIndicatorForFinestLevel(indicator, octreeHexa)
    if coarsenCoarsestLevel == 1:
        indicator = post.enforceIndicatorForCoarsestLevel(indicator, octreeHexa)

    valMin = Converter.getMinValue(indicVal,valName); epsInf = valMin/4.
    valMax = Converter.getMaxValue(indicVal,valName); epsSup = valMax*4.
    # calcul de l'indicateur : tous les pts sont raffines (low)
    indicator1 = computeIndicatorFieldForBounds(indicator, indicVal,
                                                epsInf/4., 4.*epsInf)
    res = G.adaptOctree(octreeHexa, indicator1)
    nptsfin = len(res[1][0])
    print('Number of points for low bound value %g is %d (targetPts=%d)'%(epsInf, nptsfin, nbTargetPts))
    if nptsfin < nbTargetPts: return indicator1, epsInf/4., epsInf*4.

    # calcul de l'indicateur : ts les pts sont deraffines
    indicator1 = computeIndicatorFieldForBounds(indicator, indicVal, \
                                                epsSup/4., epsSup*4.)
    res = G.adaptOctree(octreeHexa, indicator1)
    nptsfin = len(res[1][0])
    print('Number of points for high bound value %g is %d (targetPts=%d)'%(epsSup, nptsfin, nbTargetPts))
    if nptsfin > nbTargetPts:
        #print('Warning: computeIndicator: the number of final points cannot be lower than the target.')
        return indicator1, epsSup/4., epsSup*4.

    # dichotomie
    count = 0; Delta = nbTargetPts
    diffmax = 1.e-8*nbTargetPts/max(Delta,1e-6); diff = diffmax+1.
    while count < 100 and Delta > 0.02*nbTargetPts and diff > diffmax:
        eps = 0.5*(epsInf+epsSup)
        #print('epsInf =', epsInf, ' | epsSup = ', epsSup)
        indicator1 = computeIndicatorFieldForBounds(indicator, indicVal, eps/4., 4.*eps)
        res = G.adaptOctree(octreeHexa, indicator1)
        nptsfin = len(res[1][0])
        if nptsfin > nbTargetPts: epsInf = eps
        else: epsSup = eps
        Delta = abs(nbTargetPts-nptsfin)
        diffmax = 1.e-8*nbTargetPts/max(Delta, 1e-6)
        diff = abs(epsSup-epsInf)
        count += 1
        print('Number of points for bound value %g is %d (targetPts=%d)'%(eps, nptsfin, nbTargetPts))
    return indicator1, eps/4., eps*4.

#==============================================================================
def sharpEdges(array, alphaRef=30.):
    """Detect sharp edges between adjacent cells of a surface. Angle out of
    [180-alpharef,180+alpharef] are considered as sharp.
    Usage: sharpEdges(a, alpharef)"""
    try: import Generator as G; import Transform as T
    except: raise ImportError("sharpEdges: requires Generator and Transform modules.")
    if isinstance(array[0], list): # liste d'arrays
        out = []
        for i in array:
            if len(i) == 5: out.append(Converter.convertArray2Hexa(i))# structure
            else: out.append(i)
        out = T.join(out)
    else:
        if len(array) == 5: out = Converter.convertArray2Hexa(array) # structure
        else: out = array
    out = G.close(out)
    if out[3] != "BAR": out = T.reorder(out,(1,))
    res = post.sharpEdges(out, alphaRef)
    res = G.close(res)
    try: res = T.splitConnexity(res); return res
    except: return [res]

#==============================================================================
# return shape of unstructured mesh
# IN: array
# IN: vector: vector of visualization
# OUT: res: array of shapes
#==============================================================================
def silhouette(array, vector):
    """Detect shape of an unstructured surface."""
    try: import Generator as G
    except: raise ImportError("silhouette: requires Generator module.")

    if isinstance(array[0], list): # liste d'arrays
        res = []
        for a in array:
            if len(a) == 5: #structure
                a = Converter.convertArray2Hexa(a); a = G.close(a)
            else:
                a = G.close(a)
            r = post.silhouette(a,vector)
            if r is not None: res.append(r)
    else:
        if len(array) == 5: #structure
            array = Converter.convertArray2Hexa(array); array = G.close(array)
        else:
            array = G.close(array)
        res = post.silhouette(array, vector)
    return res

#==============================================================================
# modify variable names
#==============================================================================
def renameVars(array, varsPrev, varsNew):
    """Rename variables names in varsPrev with names defined by varsNew."""
    if isinstance(array[0], list): # liste d'arrays
        res = []
        for a in array:
            b = a[:]
            varsb = b[0]; varsb = varsb.split(',')
            for nov in range(len(varsPrev)):
                try:
                    pos = varsb.index(varsPrev[nov])
                    varsb[pos] = varsNew[nov]
                except: pass
            b[0] = ','.join(varsb)
            res.append(b)
        return res
    else:
        res = array[:]
        varsb = res[0]; varsb = varsb.split(',')
        for nov in range(len(varsPrev)):
            try:
                pos = varsb.index(varsPrev[nov])
                varsb[pos] = varsNew[nov]
            except: pass
        res[0] = ','.join(varsb)
    return res

## [AJ - KEEP FOR NOW - FROM MASTER]
## This function needs to be further tested & validated and should be used at your own risk
##      Further devs might occur upon further discussion with other developers
def computeIndicatorField_AMR(octreeHexa, indicVal, nbTargetPts=-1, bodies=[],
                              refineFinestLevel=1, coarsenCoarsestLevel=1,valMin=0,valMax=1,isOnlySmallest=False):
    """Compute the indicator -1, 0 or 1 for each element of the HEXA octree
    with respect to the indicatorValue field located at element centers. The
    bodies fix the indicator to 0 in the vicinity of bodies. nbTargetPts
    controls the number of points after adaptation.
    If refineFinestLevel=1, the finest levels are refined.
    If coarsenCoarsestLevel=1, the coarsest levels are coarsened wherever possible.
    Returns the indicator field.
    Usage: computeIndicatorField(octreeHexa, indicVal, nbTargetPts, bodies, refinestLevel, coarsenCoarsestLevel)"""
    try: import Generator as G
    except: raise ImportError("computeIndicatorField: requires Generator module.")
    npts = octreeHexa[1][0].shape[0]
    valName = indicVal[0]; nelts = indicVal[2]
    indicVal[1] = numpy.absolute(indicVal[1])
    indicator = Converter.initVars(indicVal, 'indicator', 0.)
    indicator = Converter.extractVars(indicator, ['indicator'])

    if bodies != []:
        bodies = Converter.convertArray2Tetra(bodies)
        indicator = post.enforceIndicatorNearBodies(indicator, octreeHexa, bodies)
    if refineFinestLevel == 0:
        indicator = post.enforceIndicatorForFinestLevel(indicator, octreeHexa)
        if isOnlySmallest: return indicator
    if coarsenCoarsestLevel == 1:
        indicator = post.enforceIndicatorForCoarsestLevel(indicator, octreeHexa)


    indicator1 = computeIndicatorFieldForBounds(indicator, indicVal,valMin,valMax,isAMR=True)
    res        = G.adaptOctree(octreeHexa, indicator1)
    nptsfin    = len(res[1][0])
    print('Number of points: Pre %d | Post %d | Increase: %f'%(npts, nptsfin, nptsfin/npts))
    if nptsfin < nbTargetPts:
        return indicator1

    count = 0
    while count < 100 and nptsfin/nbTargetPts > 1.02:
        valMean    = 0.5*(valMin+valMax)
        valMax    += 0.25*abs(valMean-valMin)
        valMin    += 0.25*abs(valMean-valMin)
        indicator1 = computeIndicatorFieldForBounds(indicator, indicVal,valMin,valMax,isAMR=True)
        res        = G.adaptOctree(octreeHexa, indicator1)
        nptsfin    = len(res[1][0])
        count     += 1
        print('Number of points: Pre %d | Post %d | Increase: %f | Lower Threshold: %f | Upper Threshold: %f | Limits Modif. Loop Cnt: %d'%(npts, nptsfin, count))

    return indicator1
