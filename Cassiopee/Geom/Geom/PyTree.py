"""Geometry definition module.
"""
import Geom
__version__ = Geom.__version__

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Converter
except ImportError:
    raise ImportError("Geom.PyTree: requires Converter.PyTree module.")
from .Offset import offsetSurface

def point(P):
    """Create a point.
    Usage: point((x,y,z))"""
    a = Geom.point(P)
    return C.convertArrays2ZoneNode('point', [a])

def cloud(arr):
    """Create a point cloud.
    Usage: cloud([(x1,x2,...,xn),(y1,y2,...,yn),(z1,z2,...,zn)])"""
    a = Geom.cloud(arr)
    return C.convertArrays2ZoneNode('cloud', [a])

def naca(e, N=101, sharpte=True):
    """Create a naca profile of N points.
    Usage: naca(e, N)"""
    a = Geom.naca(e, N, sharpte)
    zname = 'naca'
    if isinstance(e, str): zname ='NACA'+e
    return C.convertArrays2ZoneNode(zname, [a])

def line(P1, P2, N=100):
    """Create a line of N points. Usage: line( (x1,y1,z1), (x2,y2,z2), N )"""
    a = Geom.line(P1, P2, N)
    return C.convertArrays2ZoneNode('line', [a])

def polyline(Pts):
    """Create a polyline of N points.
    Usage: polyline( [(x1,y1,z1),....,(xn,yn,zn)])"""
    a = Geom.polyline( Pts )
    return C.convertArrays2ZoneNode('polyline', [a])

def circle(Center, R, tetas=0., tetae=360., N=100):
    """Create a portion of circle of N points and of center C,
    radius R, between angle tetas and tetae.
    Usage: circle((xc,yc,zc), R, tetas, tetae, N)"""
    a = Geom.circle(Center, R, tetas, tetae, N)
    return C.convertArrays2ZoneNode('circle', [a])

def bezier(t, N=100, M=100, density=-1):
    """Create a a Bezier curve controlled by an array of control points.
    Usage: bezier(tc, N)"""
    return C.TZGC1(t, 'nodes', True, Geom.bezier, N, M, density)

def spline(t, order=3, N=100, M=100, density=-1):
    """Create a spline of N points.
    Usage: spline(ctrlsPts, order, N)"""
    return C.TZGC1(t, 'nodes', True, Geom.spline, order, N, M, density)

def nurbs(t, weight='weight', order=3, N=100, M=100, density=-1):
    """Create a nurbs of N points.
    Usage: nurbs(ctrlPts, order, N)"""
    w = C.getField(weight, t)[0]
    Pts = C.getFields(Internal.__GridCoordinates__, t)[0]
    Pts = Converter.addVars([Pts,w])
    surf = Geom.nurbs(Pts, weight, order, N, M, density)
    return C.convertArrays2ZoneNode('nurbs', [surf])

def curve(f, N=100):
    """Create a curve from a user defined parametric function.
    Usage: curve(f, N)"""
    a = Geom.curve(f, N)
    return C.convertArrays2ZoneNode('curve', [a])

def cone(center, Rb, Rv, H, N=100):
    """Create cone of NxNh points and of center C, basis radius Rb,
    vertex radius Rv and height H.
    Usage: cone((xc,yc,zc), Rb, Rv, H, N)"""
    a = Geom.cone(center, Rb, Rv, H, N)
    return C.convertArrays2ZoneNode('cone', [a])

def sphere(center, R, N=100):
    """Create a sphere of Nx2N points and of center C and radius R.
    Usage: sphere((xc,yc,zc), R, N)"""
    a = Geom.sphere(center, R, N)
    return C.convertArrays2ZoneNode('sphere', [a])

def sphere6(center, R, N=100, ntype='STRUCT'):
    """Create a shpere of NxN points and of center C and radius R.
    Usage: sphere((xc,yc,zc), R, N)"""
    A = Geom.sphere6(center, R, N, ntype)
    if ntype == 'STRUCT':
        return [C.convertArrays2ZoneNode('sphere-part1', [A[0]]),
                C.convertArrays2ZoneNode('sphere-part2', [A[1]]),
                C.convertArrays2ZoneNode('sphere-part3', [A[2]]),
                C.convertArrays2ZoneNode('sphere-part4', [A[3]]),
                C.convertArrays2ZoneNode('sphere-part5', [A[4]]),
                C.convertArrays2ZoneNode('sphere-part6', [A[5]])]
    else:
        return C.convertArrays2ZoneNode('sphere', [A])

def sphereYinYang(center, R, N=100, ntype='STRUCT'):
    """Create a sphere of center C and radius R made of two overlapping zones.
    Usage: sphereYinYang((xc,yc,zc), R, N)"""
    A = Geom.sphereYinYang(center, R, N)
    if ntype == 'STRUCT':
        return [C.convertArrays2ZoneNode('sphere-part1', [A[0]]),
                C.convertArrays2ZoneNode('sphere-part2', [A[1]])]
    else:
        return C.convertArrays2ZoneNode('sphere', [A])

def disc(center, R, N=100, ntype='STRUCT'):
    """Create a disc of center C and radius R.
    Usage: disc((xc,yc,zc), R, N)"""
    A = Geom.disc(center, R, N, ntype)
    if ntype == 'STRUCT':
        return [C.convertArrays2ZoneNode('disc-part1', [A[0]]),
                C.convertArrays2ZoneNode('disc-part2', [A[1]]),
                C.convertArrays2ZoneNode('disc-part3', [A[2]]),
                C.convertArrays2ZoneNode('disc-part4', [A[3]]),
                C.convertArrays2ZoneNode('disc-part5', [A[4]])]
    else:
        return C.convertArrays2ZoneNode('disc', [A])

def box(Pmin, Pmax, N=100, ntype='STRUCT'):
    """Create a box passing by Pmin and Pmax.
    Usage: box(Pmin, Pmax, N)"""
    A = Geom.box(Pmin, Pmax, N, ntype)
    if ntype == 'STRUCT':
        return [C.convertArrays2ZoneNode('box-part1', [A[0]]),
                C.convertArrays2ZoneNode('box-part2', [A[1]]),
                C.convertArrays2ZoneNode('box-part3', [A[2]]),
                C.convertArrays2ZoneNode('box-part4', [A[3]]),
                C.convertArrays2ZoneNode('box-part5', [A[4]]),
                C.convertArrays2ZoneNode('box-part6', [A[5]])]
    else:
        return C.convertArrays2ZoneNode('box', [A])

def cylinder(center, R, H, N=100, ntype='STRUCT'):
    """Create a cylinder of center C, radius R and hieght H.
    Usage: cylinder((xc,yc,zc), R, H, N)"""
    A = Geom.cylinder(center, R, H, N, ntype)
    if ntype == 'STRUCT':
        return [C.convertArrays2ZoneNode('cyl-part1', [A[0]]),
                C.convertArrays2ZoneNode('cyl-part2', [A[1]]),
                C.convertArrays2ZoneNode('cyl-part3', [A[2]]),
                C.convertArrays2ZoneNode('cyl-part4', [A[3]]),
                C.convertArrays2ZoneNode('cyl-part5', [A[4]]),
                C.convertArrays2ZoneNode('cyl-part6', [A[5]]),
                C.convertArrays2ZoneNode('cyl-part7', [A[6]]),
                C.convertArrays2ZoneNode('cyl-part8', [A[7]]),
                C.convertArrays2ZoneNode('cyl-part9', [A[8]]),
                C.convertArrays2ZoneNode('cyl-part10', [A[9]]),
                C.convertArrays2ZoneNode('cyl-part11', [A[10]])]
    else:
        return C.convertArrays2ZoneNode('cyl', [A])

def torus(center, R, r, alphas=0., alphae=360.,
          betas=0., betae=360., NR=100, Nr=100):
    """Create NRxNr points lying on a torus of center C and radii R (main)
    and r (tube) between the angles alphas and alphae (XY-plane) and
    between betas and betae (RZ-plane).
    Usage: torus((xc,yc,zc), R, r, NR, Nr, alphas, alphae, betas, betae)"""
    a = Geom.torus(center, R, r, alphas, alphae, betas, betae, NR, Nr)
    return C.convertArrays2ZoneNode('torus', [a])

def triangle(P1, P2, P3, N=0, ntype='TRI'):
    """Create a single triangle with points P1, P2, P3.
    Usage: triangle((x1,y,1,z1), (x2,y2,z2), (x3,y3,z3))"""
    a = Geom.triangle(P1, P2, P3, N, ntype)
    if ntype == 'STRUCT':
        return [C.convertArrays2ZoneNode('tri-part1', [a[0]]),
                C.convertArrays2ZoneNode('tri-part2', [a[1]]),
                C.convertArrays2ZoneNode('tri-part3', [a[2]])]
    else:
        return C.convertArrays2ZoneNode('triangle', [a])

def quadrangle(P1, P2, P3, P4, N=0, ntype='QUAD'):
    """Create a single quadrangle with points P1, P2, P3, P4.
    Usage: quadrangle((x1,y,1,z1), (x2,y2,z2), (x3,y3,z3), (x4,y4,z4))"""
    a = Geom.quadrangle(P1, P2, P3, P4, N, ntype)
    if ntype == 'STRUCT':
        return [C.convertArrays2ZoneNode('quad-part1', [a[0]]),
                C.convertArrays2ZoneNode('quad-part2', [a[1]]),
                C.convertArrays2ZoneNode('quad-part3', [a[2]]),
                C.convertArrays2ZoneNode('quad-part2', [a[3]])]
    else:
        return C.convertArrays2ZoneNode('quadrangle', [a])

def surface(f, N=100):
    """Create a surface from a user defined parametric function.
    Usage: surface(f, N)"""
    a = Geom.surface(f, N)
    return C.convertArrays2ZoneNode('surface', [a])

def getLength(t):
    """Return the length of 1D array(s) defining a mesh.
    Usage: getLength(t)"""
    coords = C.getFields(Internal.__GridCoordinates__, t)
    return Geom.getLength(coords)

def getDistantIndex(t, ind, l):
    """Return the index of 1D array defining a mesh located at a
    distance l of ind.
    Usage: getDistantIndex(t, ind, l)"""
    a = C.getFields(Internal.__GridCoordinates__, t)[0]
    return Geom.getDistantIndex(a, ind, l)

def getNearestPointIndex(t, pointList):
    """Return the nearest index of points in array.
    Usage: getNearestPointIndex(t, pointList)"""
    a = C.getFields(Internal.__GridCoordinates__, t)[0]
    return Geom.getNearestPointIndex( a, pointList )

def getCurvatureHeight(t):
    """Return the curvature height for each point.
    Usage: getCurvatureHeight(t)"""
    return C.TZGC1(t, 'nodes', True, Geom.getCurvatureHeight)

def _getCurvatureHeight(t):
    """Return the curvature height for each point."""
    return C._TZGC1(t, 'nodes', False, Geom.getCurvatureHeight)

def getCurvatureRadius(t):
    """Return the curvature radius for each point.
    Usage: getCurvatureRadius(t)"""
    return C.TZGC1(t, 'nodes', True, Geom.getCurvatureRadius)

def _getCurvatureRadius(t):
    """Return the curvature radius for each point."""
    return C._TZGC1(t, 'nodes', False, Geom.getCurvatureRadius)

def getCurvatureAngle(t):
    """Return the curvature angle for each point.
    Usage: getCurvatureAngle(t)"""
    return C.TZGC1(t, 'nodes', True, Geom.getCurvatureAngle)

def _getCurvatureAngle(t):
    """Return the curvature angle for each point."""
    return C._TZGC1(t, 'nodes', False, Geom.getCurvatureAngle)

def getSharpestAngle(t):
    """Return the sharpest angle for each point of a surface based on the sharpest angle
    between adjacent element to which the point belongs to.
    Usage: getSharpestAngle(a)"""
    return C.TZGC1(t, 'nodes', True, Geom.getSharpestAngle)

def _getSharpestAngle(t):
    """Return the sharpest angle for each point of a surface based on the sharpest angle."""
    return C._TZGC1(t, 'nodes', False, Geom.getSharpestAngle)

def getCurvilinearAbscissa(t):
    """Return the curvilinear abscissa for each point.
    Usage: getCurvilinearAbscissa(t)"""
    return C.TZGC1(t, 'nodes', True, Geom.getCurvilinearAbscissa)

def _getCurvilinearAbscissa(t):
    """Return the curvilinear abscissa for each point.
    Usage: getCurvilinearAbscissa(t)"""
    return C._TZGC1(t, 'nodes', False, Geom.getCurvilinearAbscissa)

def getDistribution(t):
    """Return the curvilinear abscissa for each point as coordinates.
    Usage: getDistribution(t)"""
    return C.TZGC1(t, 'nodes', True, Geom.getDistribution)

def getTangent(t):
    """Makes the tangent of a 1D curve. The input argument shall be a structured
    1D curve. Each node of the output represents the unitary tangent vector, 
    pointing towards the tangent direction of the input 1D curve.
    Usage: b = getTangent(t)"""
    tp = Internal.copyRef(t)
    C._deleteFlowSolutions__(tp)
    C._TZGC1(tp, 'nodes', True, Geom.getTangent)
    return tp

def addSeparationLine(t, line0):
    """Add a separation line defined in line0 to a mesh defined in t.
    Usage: addSeparationLine(t, line0)"""
    al = C.getFields(Internal.__GridCoordinates__, line0, api=3)[0]
    at = C.getFields(Internal.__GridCoordinates__, t, api=3)[0]
    arrays = Geom.addSeparationLine(at, al)
    zones = []
    for i in arrays:
        zone = C.convertArrays2ZoneNode(t[0], [i])
        zones.append(zone)
    return zones

# Obsolete
def lineGenerate(t, line):
    """Obsolete function."""
    return lineDrive(t, line)

def lineDrive(t, line):
    """Generate a surface mesh by using 1D array (defining a mesh)
    and following the curve defined in line.
    Usage: lineDrive(t, line)"""
    al = C.getFields(Internal.__GridCoordinates__, line)
    if len(al) == 1: al = al[0]
    al2 = Converter.node2Center(al)
    # Attention les coord. des centres ne sont pas justes! mais
    # elles ne sont pas utilisees dans la fonction
    return C.TZAGC(t, 'both', 'both', True, Geom.lineDrive,
                   Geom.lineDrive, al, al2)

def orthoDrive(t, line, mode=0):
    """Generate a surface mesh by using 1D array (defining a mesh)
    and following orthogonally the curve defined in line.
    Usage: orthoDrive(t, line)"""
    al = C.getFields(Internal.__GridCoordinates__, line)
    if len(al) == 1: al = al[0]
    al2 = Converter.node2Center(al)
    # Attention les coord. des centres ne sont pas justes! mais
    # elles ne sont pas utilisees dans la fonction
    return C.TZAGC(t, 'both', 'both', True, Geom.orthoDrive,
                   Geom.orthoDrive, al, mode, al2, mode)

def axisym(t, center, axis, angle=360., Ntheta=360, rmod=None):
    """Create an axisymmetric mesh given an azimuthal surface mesh.
    Usage: axisym(t, (xo,yo,zo), (nx,ny,nz), teta, Nteta, rmod)"""
    # Attention en centres, les coord. des centres ne sont pas justes! mais
    # elles ne sont pas utilisees dans la fonction
    if rmod is not None: rmod = C.getFields(Internal.__GridCoordinates__, rmod)[0]
    return C.TZAGC(t, 'both', 'both', True, Geom.axisym, Geom.axisym,
                   center, axis, angle, Ntheta, rmod,
                   center, axis, angle, Ntheta-1, rmod)

def _axisym(t, center, axis, angle=360., Ntheta=360, rmod=None):
    if rmod is not None: rmod = C.getFields(Internal.__GridCoordinates__, rmod)[0]
    return C._TZAGC(t, 'both', 'both', True, Geom.axisym, Geom.axisym,
                    center, axis, angle, Ntheta, rmod,
                    center, axis, angle, Ntheta-1, rmod)

def volumeFromCrossSections(t):
    """Generate a 3D volume from cross sections contours in (x,y) planes.
    Usage: volumeFromCrossSections(t)"""
    tp = Internal.copyRef(t)
    nodes = Internal.getZones(tp)
    coords = []
    for z in nodes:
        coords.append(C.getFields(Internal.__GridCoordinates__, z)[0])
    coordp = Geom.volumeFromCrossSections(coords)
    zone = C.convertArrays2ZoneNode('Volume', [coordp])
    return zone

def text1D(string, font='vera', smooth=0, offset=0.5):
    """Create a 1D text. offset is the space between letters,
    font indicates used font, smooth indicates the smoothing intensity.
    Usage: text1D(string, font, smooth, offset)"""
    a = Geom.text1D(string, font, smooth, offset)
    zones = []
    for i in a:
        zone = C.convertArrays2ZoneNode(string, [i])
        zones.append(zone)
    return zones

def text2D(string, font='vera', smooth=0, offset=0.5):
    """Create a 2D text. offset is the space between letters.
    font indicates used font, smooth indicates the smoothing intensity.
    Usage: text2D(string, font, smooth, offset)"""
    a = Geom.text2D(string, font, smooth, offset)
    return C.convertArrays2ZoneNode(string, [a])

def text3D(string, font='vera', smooth=0, offset=0.5):
    """Create a 3D text. offset is the space between letters.
    font indicates used font, smooth indicates the smoothing intensity.
    Usage: text3D(string, font, smooth, offset)"""
    a = Geom.text3D(string, font, smooth, offset)
    return C.convertArrays2ZoneNode(string, [a])

def connect1D(curves, sharpness=0, N=10, lengthFactor=1.):
    """Connect curves with sharp or smooth junctions.
    Usage: a = connect1D(A, sharpness=0)"""
    a = C.getFields(Internal.__GridCoordinates__, curves)
    z = Geom.connect1D(a, sharpness, N, lengthFactor)
    return C.convertArrays2ZoneNode('connected', [z])

def uniformize(a, N=100, h=-1, factor=-1., density=-1., sharpAngle=30.):
    """Uniformize a 1D curve."""
    ap = Internal.copyRef(a)
    _uniformize(ap, N, h, factor, density, sharpAngle)
    return ap

def _uniformize(a, N=100, h=-1., factor=-1, density=-1, sharpAngle=30.):
    C._deleteFlowSolutions__(a)
    ar = C.getFields(Internal.__GridCoordinates__, a)
    ar = Geom.uniformize(ar, N, h, factor, density, sharpAngle)
    C.setFields(ar, a, 'nodes')
    return None

def refine(a, N=10, factor=-1, sharpAngle=30.):
    """Refine a 1D curve."""
    ap = Internal.copyRef(a)
    _refine(ap, N, factor, sharpAngle)
    return ap

def _refine(a, N=10, factor=-1, sharpAngle=30.):
    """Refine a 1D curve."""
    C._deleteFlowSolutions__(a)
    C._TZGC1(a, 'nodes', True, Geom.refine, N, factor, sharpAngle)
    return None

def distrib1(a, h, normalized=True):
    """Enforce h everywhere in line. Return distribution."""
    ap = Internal.copyRef(a)
    _distrib1(ap, h, normalized)
    return ap

def _distrib1(a, h, normalized=True):
    """Enforce h everywhere in line. Return distribution."""
    C._TZA1(a, 'nodes', 'nodes', True, Geom.distrib1, h, normalized)
    return None

def distrib2(a, h1, h2, add=20, forceAdd=False, normalized=True, algo=0):
    """Enforce h1,h2 in line. Return distribution."""
    ap = Internal.copyRef(a)
    _distrib2(ap, h1, h2, add, forceAdd, normalized, algo)
    return ap

def _distrib2(a, h1, h2, add=20, forceAdd=False, normalized=True, algo=0):
    """Enforce h1,h2 in line. Return distribution."""
    C._TZA1(a, 'nodes', 'nodes', True, Geom.distrib2, h1, h2, add, forceAdd, normalized, algo)
    return None

def enforceh(a, N=100, h=-1.):
    """Remesh a 1D curve with imposed steps."""
    ap = Internal.copyRef(a)
    _enforceh(ap, N, h)
    return ap

def _enforceh(a, N=100, h=-1.):
    """Remesh a 1D curve with imposed steps."""
    C._TZA1(a, 'nodes', 'nodes', True, Geom.enforceh, N, h)
    C._deleteFlowSolutions__(a, loc='centers')
    C._rmVars(a, ['h','f']) # h,f is no longer correct
    return None

def setH(a, ind, h):
    """Set h step indicator for enforceh."""
    _setH(a, ind, h)

def _setH(a, ind, h):
    """Set h step indicator for enforceh."""
    zones = Internal.getZones(a)
    for z in zones:
        if C.isNamePresent(z, 'h') == -1: C._initVars(a, 'h', 0.)
        C.setValue(a, 'h', ind, h)

def setF(a, ind, f):
    """Set f factor indicator for enforceh."""
    _setF(a, ind, f)

def _setF(a, ind, f):
    """Set f factor indicator for enforceh."""
    zones = Internal.getZones(a)
    for z in zones:
        if C.isNamePresent(z, 'f') == -1: C._initVars(a, 'f', 0.)
        C.setValue(a, 'f', ind, f)

def smooth(a, eps, niter):
    """Smooth distribution on a curve."""
    ap = Internal.copyRef(a)
    _smooth(ap, eps, niter)
    return ap

def _smooth(a, eps, niter):
    """Smooth distribution on a curve."""
    C._TZGC1(a, 'nodes', True, Geom.smooth, eps, niter)
    return None

def getUV(a, normalDeviationWeight=2., texResolution=1920, fields=None):
    """Return uv of surface and atlas."""
    b = Internal.getZones(a)[0] # only first zone for now
    array = C.getFields('nodes', b, api=3)[0]
    ret = Geom.getUV(array, normalDeviationWeight, texResolution, fields)
    z0 = C.convertArrays2ZoneNode(b[0], [ret[0]])
    z1 = C.convertArrays2ZoneNode('colorAtlas', [ret[1]])
    z2 = C.convertArrays2ZoneNode('bumpAtlas', [ret[2]])
    if len(ret) > 3: z3 = C.convertArrays2ZoneNode('speedAtlas', [ret[3]])
    return z0, z1, z2

def getUVFromIJ(t):
    tp = Internal.copyRef(t)
    _getUVFromIJ(tp)
    return tp

def _getUVFromIJ(a):
    zones = Internal.getZones(a)
    for z in zones:
        C._initVars(z, '_u_=0.')
        C._initVars(z, '_v_=0.')
        pu = Internal.getNodeFromName2(z, '_u_')
        pv = Internal.getNodeFromName2(z, '_v_')
        pu = pu[1].ravel('k')
        pv = pv[1].ravel('k')
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured':
            ni = dim[1]; nj = dim[2]
            for j in range(nj):
                for i in range(ni):
                    pu[i+j*ni] = j*1./(nj-1)
                    pv[i+j*ni] = i*1./(ni-1)
    return None