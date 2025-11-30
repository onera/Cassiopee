"""Geometry definition module.
"""
__version__ = '4.1'
__author__ = "Stephanie Peron, Christophe Benoit, Pascal Raud, Sam Landier"

from . import geom
import numpy
import KCore.Vector as Vector

from .MapEdge import enforceh, uniformize, refine, setH, setF, enforce, distrib1, distrib2, smooth, mapCurvature, enforceh3D

# - Basic entities -
def point(P):
    """Create a point. 
    Usage: a = point((x,y,z))"""
    a = numpy.zeros((3, 1), numpy.float64)
    a[0,0] = P[0]; a[1,0] = P[1]; a[2,0] = P[2]
    c = numpy.ones((1, 0), numpy.int32)
    return ['x,y,z', a, c, 'NODE']

def cloud(arr):
    """Create a point cloud. 
    Usage: a = cloud([(x1,x2,...,xn),(y1,y2,...,yn),(z1,z2,...,zn)])"""
    if not isinstance(arr, numpy.ndarray):
        arr = numpy.asarray(arr, numpy.float64)
    c = numpy.ones((arr.shape[1], 0), numpy.int32)
    return ['x,y,z', arr, c, 'NODE']

def naca(e, N=101, sharpte=True):
    """Create a NACA00xx profile of N points and thickness e. 
    Usage: a = naca(e, N)"""
    im = -1; ip = -1; it = -1; ith = -1; iq = -1
    if isinstance(e, str): # digits
        if len(e) == 4:
            im = int(e[0:1])
            ip = int(e[1:2])
            it = int(e[2:4])
        elif len(e) == 5:
            im = int(e[0:1]) # il
            ip = int(e[1:2])
            iq = int(e[2:3])
            it = int(e[3:5])
        elif len(e) == 7:
            im = int(e[0:1])
            ip = int(e[1:2])
            ith = int(e[2:4])
            iq = int(e[5:6]) # ii
            it = int(e[6:7])
        e = 0.
    if not sharpte: sharpte = 0
    else: sharpte = 1
    return geom.naca(e, N, im, ip, it, ith, iq, sharpte)

def profile(name=None):
    """Create a wing profile mesh."""
    import os
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import KCore.installPath
    prodname = os.getenv("ELSAPROD")
    libpath = KCore.installPath.libPath
    pos = libpath.rfind(prodname)
    filepath = os.path.join(libpath[:pos], prodname, "UIUCAirfoils.cgns")
    t = C.convertFile2PyTree(filepath)
    bases = Internal.getBases(t)
    if name is None:
        # print a dictionary
        for b in bases:
            for z in Internal.getZones(b):
                print("%s/%s\n"%(b[0],z[0]))
        return None
    name = name.split('/')
    if len(name) == 1:
        # print a dictionary
        b = Internal.getNodeFromName1(t, name[0])
        if b is not None:
            for z in Internal.getZones(b):
                print("%s/%s\n"%(b[0],z[0]))
        return None
    if len(name) != 2:
        raise ValueError('profile: name must be base/name.')
    b = Internal.getNodeFromName1(t, name[0])
    if b is None:
        raise ValueError('profile: base name not found.')
    z = Internal.getNodeFromName1(b, name[1])
    if z is None:
        raise ValueError('profile: zone name not found.')
    a = C.getAllFields(z, 'nodes', api=1)[0]
    return a

def line(P1, P2, N=100):
    """Create a line of N points. 
    Usage: a = line((x1,y1,z1), (x2,y2,z2), N)"""
    if len(P1) == 2: P1 = (P1[0],P1[1],0.)
    if len(P2) == 2: P2 = (P2[0],P2[1],0.)
    return geom.line(P1, P2, N)

def spline(Pts, order=3, N=100, M=100, density=-1):
    """Create a spline of N points. 
    Usage: a = spline(Pts, order, N, M, density)"""
    return geom.spline(Pts, order, N, order, M, density)

def nurbs(Pts, weight='weight', order=3, N=100, M=100, density=-1):
    """Create a nurbs of N points.
    Usage: a = nurbs(ctrlsPts, order, N)"""
    try:
        import Converter
        Weights = Converter.extractVars(Pts,[weight])
    except:
        ni = Pts[2]; nj = Pts[3]; nk = Pts[4]
        Weights = [weight, numpy.ones((1,ni*nj*nk), numpy.float64), ni, nj, nk]
    return geom.nurbs(Pts, Weights, order, N, order, M, density)

def polyline(Pts):
    """Create a polyline of N points. 
    Usage: a = polyline([(x1,y1,z1),....,(xn,yn,zn)])"""
    return geom.polyline(Pts)

def circle(C, R, tetas=0., tetae=360., N=100):
    """Create a portion of circle of N points and of center C, radius R, between angle tetas and tetae.
    Usage: a = circle((xc,yc,zc), R, tetas, tetae, N)"""
    if len(C) == 2: C = (C[0],C[1],0.)
    return geom.circle(C, R, tetas, tetae, N)

def sphere(C, R, N=100):
    """Create a sphere of Nx2N points and of center C and radius R.
    Usage: a = sphere((xc,yc,zc), R, N)"""
    return geom.sphere(C, R, N)

def sphere6(C, R, N=100, ntype='STRUCT'):
    """Create a sphere of 6NxN points and of center C and radius R, made of 6 parts.
    Usage: a = sphere6((xc,yc,zc), R, N)"""
    try: import Transform as T; import Generator as G
    except ImportError:
        raise ImportError("sphere6: requires Transform and Generator modules.")
    s = sphere(C, R, 2*N)

    (x0,y0,z0) = (-R/2.+C[0],-R/2.+C[1],-R/2.+C[2])
    (hx,hy,hz) = (R/(N-1.),R/(N-1.),R/(N-1.))
    b1 = G.cart((x0,y0,z0), (hx,hy,hz), (N,N,1)) # k=1
    b2 = G.cart((x0,y0,z0), (hx,hy,hz), (1,N,N)) # i=1
    b3 = G.cart((x0,y0,z0), (hx,hy,hz), (N,1,N)) # j=1
    b4 = G.cart((x0+(N-1)*hx,y0,z0), (hx,hy,hz), (1,N,N)) # i=imax
    b5 = G.cart((x0,y0,z0+(N-1)*hz), (hx,hy,hz), (N,N,1)) # k=kmax
    b6 = G.cart((x0,y0+(N-1)*hy,z0), (hx,hy,hz), (N,1,N)) # j=jmax

    c1 = T.projectRay(b1, [s], C); c1 = T.reorder(c1, (-1,2,3))
    c2 = T.projectRay(b2, [s], C); c2 = T.reorder(c2, (3,2,1))
    c3 = T.projectRay(b3, [s], C); c3 = T.reorder(c3, (1,3,2))
    c4 = T.projectRay(b4, [s], C); c4 = T.reorder(c4, (3,1,2))
    c5 = T.projectRay(b5, [s], C)
    c6 = T.projectRay(b6, [s], C); c6 = T.reorder(c6, (-1,3,2))
    m = [c1, c2, c3, c4, c5, c6]
    return export__(m, ntype)

def sphereYinYang(C, R, N=100, ntype='STRUCT'):
    """Create a sphere of center C and radius R made of two overset parts.
    Usage: a = sphereYinYang((xc,yc,zc), R, N)"""
    try: import Transform as T
    except ImportError:
        raise ImportError("sphereYinYang: requires Transform module.")
    Ni = 4*(N//2)
    a = sphere(C, R, N=Ni)
    a = T.subzone(a, (Ni//4-2,Ni//4-2,1), (3*Ni//4+2,7*Ni//4+2,1))
    b = T.rotate(a, (0,0,0), (0,1,0), 90.)
    b = T.rotate(b, (0,0,0), (0,0,1), 180.)
    m = [a, b]
    return export__(m, ntype)

# IN: liste de maillages struct
def export__(a, ntype='STRUCT'):
    try: import Converter as C; import Transform as T; import Generator as G
    except ImportError:
        raise ImportError("export: requires Converter, Generator and Transform modules.")
    if ntype == 'STRUCT': return a
    elif ntype == 'QUAD':
        a = C.convertArray2Hexa(a)
        a = T.join(a)
        a = G.close(a)
        return a
    elif ntype == 'TRI':
        a = C.convertArray2Tetra(a)
        a = T.join(a)
        a = G.close(a)
        return a

def disc(C, R, N=100, ntype='STRUCT'):
    """Create a disc of center C and radius R made of 5 parts.
    Usage: a = disc((xc,yc,zc), R, N)"""
    try: import Generator as G; import Transform as T; import math
    except ImportError:
        raise ImportError("disc: requires Generator and Transform module.")
    coeff = R*math.sqrt(2.)*0.25
    x = C[0]; y = C[1]; z = C[2]
    c = circle(C, R, tetas=-45., tetae=45., N=N)
    l1 = line((x+coeff,y-coeff,z), (x+coeff,y+coeff,z), N=N)
    l2 = line((x+coeff,y-coeff,z), (x+2*coeff,y-2*coeff,z), N=N)
    l3 = line((x+coeff,y+coeff,z), (x+2*coeff,y+2*coeff,z), N=N)
    m1 = G.TFI([c, l1, l2, l3])

    c = circle(C, R, tetas=45., tetae=45.+90., N=N)
    l1 = line((x+coeff,y+coeff,z), (x-coeff,y+coeff,z), N=N)
    l2 = line((x+coeff,y+coeff,z), (x+2*coeff,y+2*coeff,z), N=N)
    l3 = line((x-coeff,y+coeff,z), (x-2*coeff,y+2*coeff,z), N=N)
    m2 = G.TFI([c, l1, l2, l3])

    c = circle(C, R, tetas=45.+90, tetae=45.+180., N=N)
    l1 = line((x-coeff,y+coeff,z), (x-coeff,y-coeff,z), N=N)
    l2 = line((x-coeff,y+coeff,z), (x-2*coeff,y+2*coeff,z), N=N)
    l3 = line((x-coeff,y-coeff,z), (x-2*coeff,y-2*coeff,z), N=N)
    m3 = G.TFI([c, l1, l2, l3])

    c = circle(C, R, tetas=45.+180, tetae=45.+270., N=N)
    l1 = line((x-coeff,y-coeff,z), (x+coeff,y-coeff,z), N=N)
    l2 = line((x-coeff,y-coeff,z), (x-2*coeff,y-2*coeff,z), N=N)
    l3 = line((x+coeff,y-coeff,z), (x+2*coeff,y-2*coeff,z), N=N)
    m4 = G.TFI([c, l1, l2, l3])

    h = 2*coeff/(N-1)
    m5 = G.cart((x-coeff,y-coeff,z), (h,h,h), (N, N, 1))
    m5 = T.reorder(m5, (-1,2,3))
    m = [m1,m2,m3,m4,m5]
    return export__(m, ntype)

def quadrangle(P1, P2, P3, P4, N=0, ntype='QUAD'):
    """Create a single quadrangle with points P1, P2, P3, P4.
    Usage: a = quadrangle((x1,y,1,z1), (x2,y2,z2), (x3,y3,z3), (x4,y4,z4))"""
    try: import Generator as G
    except ImportError:
        raise ImportError("quadrangle: requires Generator module.")
    if N == 0 and ntype == 'QUAD': return geom.quadrangle(P1, P2, P3, P4)
    l1 = line(P1, P2, N)
    l2 = line(P2, P3, N)
    l3 = line(P3, P4, N)
    l4 = line(P4, P1, N)
    m = [G.TFI([l1, l2, l3, l4])]
    return export__(m, ntype)

def triangle(P0, P1, P2, N=0, ntype='TRI'):
    """Create a triangle made of 3 parts.
    Usage: a = triangle(P1, P2, P3, N)"""
    if len(P0) == 2: P0 = (P0[0],P0[1],0.)
    if len(P1) == 2: P1 = (P1[0],P1[1],0.)
    if len(P2) == 2: P2 = (P2[0],P2[1],0.)
    if N == 0 and ntype == 'TRI': return geom.triangle(P0, P1, P2)
    try: import Generator as G; import Transform as T
    except ImportError:
        raise ImportError("triangle: requires Generator and Transform module.")
    C01 = (0.5*(P0[0]+P1[0]), 0.5*(P0[1]+P1[1]), 0.5*(P0[2]+P1[2]))
    C12 = (0.5*(P1[0]+P2[0]), 0.5*(P1[1]+P2[1]), 0.5*(P1[2]+P2[2]))
    C02 = (0.5*(P0[0]+P2[0]), 0.5*(P0[1]+P2[1]), 0.5*(P0[2]+P2[2]))
    C = (1./3.*(P0[0]+P1[0]+P2[0]),
         1./3.*(P0[1]+P1[1]+P2[1]),
         1./3.*(P0[2]+P1[2]+P2[2]))

    l1 = line(P0, C01, N)
    l2 = line(C01, C, N)
    l3 = line(C, C02, N)
    l4 = line(C02, P0, N)
    m1 = G.TFI([l1, l2, l3, l4])
    m1 = T.reorder(m1, (-1,2,3))

    l1 = line(C01, P1, N)
    l2 = line(P1, C12, N)
    l3 = line(C12, C, N)
    l4 = line(C, C01, N)
    m2 = G.TFI([l1, l2, l3, l4])
    m2 = T.reorder(m2, (-1,2,3))

    l1 = line(C, C12, N)
    l2 = line(C12, P2, N)
    l3 = line(P2, C02, N)
    l4 = line(C02, C, N)
    m3 = G.TFI([l1, l2, l3, l4])
    m3 = T.reorder(m3, (-1,2,3))
    m = [m1, m2, m3]
    return export__(m, ntype)

def box(Pmin, Pmax, N=100, ntype='STRUCT'):
    """Create a box passing by Pmin and Pmax (axis aligned)."""
    try: import Generator as G; import Transform as T
    except ImportError:
        raise ImportError("box: requires Generator and Transform module.")
    N = max(N, 2)
    (xmin,ymin,zmin) = Pmin
    (xmax,ymax,zmax) = Pmax
    hx = abs(xmax-xmin)/(N-1.)
    hy = abs(ymax-ymin)/(N-1.)
    hz = abs(zmax-zmin)/(N-1.)

    s1 = G.cart(Pmin, (hx,hy,hz), (N, N, 1))
    s1 = T.reorder(s1, (-1,2,3))
    s2 = G.cart((xmin,ymin,zmax), (hx,hy,hz), (N, N, 1))
    s3 = G.cart( Pmin, (hx,hy,hz), (N, 1, N) )
    s3 = T.reorder(s3, (2,1,3))
    s4 = G.cart((xmin,ymax,zmin), (hx,hy,hz), (N, 1, N))
    s4 = T.reorder(s4, (-2,1,3))
    s5 = G.cart(Pmin, (hx,hy,hz), (1, N, N))
    s5 = T.reorder(s5, (1,-2,3))
    s6 = G.cart((xmax,ymin,zmin), (hx,hy,hz), (1, N, N))
    s = [s1, s2, s3, s4, s5, s6]
    return export__(s, ntype)

def cylinder(C, R, H, N=100, ntype='STRUCT'):
    """Create a cylinder of center C, radius R and height H."""
    try: import Transform as T
    except ImportError:
        raise ImportError("cylinder: requires Generator and Transform module.")
    (x0,y0,z0) = C
    m0 = disc(C, R, N)
    m1 = disc((x0,y0,z0+H), R, N)
    m1 = T.reorder(m1, (-1,2,3))
    m2 = circle(C, R, tetas=-45, tetae=-45+360, N=4*N-3)
    l = line(C, (x0,y0,z0+H), N=N)
    m2 = lineDrive(m2, l)
    s = m0 + m1 + [m2]
    return export__(s, ntype)

def cone(C, Rb, Rv, H, N=100):
    """Create a cone of NxN points and of center C, basis radius Rb, vertex radius Rv and height H.
    Usage: a = cone((xc,yc,zc), Rb, Rv, H, N)"""
    return geom.cone(C, Rb, Rv, H, N)

def torus(C, R, r, alphas=0., alphae=360.,
          betas=0., betae=360., NR=100, Nr=100):
    """Create a surface mesh of a torus made by NRxNr points.
    Usage: a = torus((xc,yc,zc), R, r, alphas, alphae, betas, betae, NR, Nr)"""
    return geom.torus(C, R, r, alphas, alphae, betas, betae, NR, Nr)

def bezier(controlPts, N=100, M=100, density=-1):
    """Create a a Bezier curve defined by an array of control points controlPts.
    Usage: a = bezier(controlPts, N, M)"""
    return geom.bezier(controlPts, N, M, density)

def curve(f, N=100):
    """Create a curve from a user defined parametric function or a formula.
    Usage: a = curve(f, N)"""
    if isinstance(f, str): return curve_(f, N)
    else: return curve__(f, N)

# Courbe parametree a partir d'une formule
def curve_(f, N):
    import Converter; import Generator
    a = Generator.cart( (0,0,0), (1./(N-1),1,1), (N,1,1))
    a[0] = 't,y,z'
    a = Converter.initVars(a, f)
    a = Converter.extractVars(a, ['x','y','z'])
    return a

# Courbe parametree a partir d'une fonction
def curve__(f, N):
    a = numpy.zeros((3, N), dtype=numpy.float64)
    r = f(0)
    if len(r) != 3:
        print("Warning: curve: parametric function must return a (x,y,z) tuple.")
        return ['x,y,z', a, N, 1, 1]
    for i in range(N):
        t = 1.*i/(N-1)
        r = f(t)
        a[0,i] = r[0]; a[1,i] = r[1]; a[2,i] = r[2]
    return ['x,y,z', a, N, 1, 1]

def surface(f, N=100):
    """Create a surface from a user defined parametric function or a formula.
    Usage: a = surface(f, N)"""
    if isinstance(f, str): return surface_(f, N)
    else: return surface__(f, N)

# Surface parametree a partir d'une formule
def surface_(f, N):
    import Converter; import Generator
    a = Generator.cart( (0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
    a[0] = 't,u,z'
    a = Converter.initVars(a, f)
    a = Converter.extractVars(a, ['x','y','z'])
    return a

# Surface parametree a partir d'une fonction
def surface__(f, N):
    a = numpy.zeros((3, N*N), dtype=numpy.float64)
    r = f(0,0)
    if len(r) != 3:
        print("Warning: surface: parametric function must return a (x,y,z) tuple.")
        return ['x,y,z', a, N, N, 1]
    for j in range(N):
        u = 1.*j/(N-1)
        for i in range(N):
            ind = i + j*N
            t = 1.*i/(N-1)
            r = f(t,u)
            a[0,ind] = r[0]
            a[1,ind] = r[1]
            a[2,ind] = r[2]
    return ['x,y,z', a, N, N, 1]

# - informations -
def getLength(a):
    """Return the length of 1D-mesh.
    Usage: l = getLength(a"""
    if isinstance(a[0], list):
        l = 0.
        for i in a: l += geom.getLength(i)
        return l
    else: return geom.getLength(a)

def getDistantIndex(a, ind, l):
    """Return the index of 1D-mesh located at a distance l of ind.
    Usage: ind = getDistantIndex(a, ind, l)"""
    return geom.getDistantIndex(a, ind, l)

def getNearestPointIndex(a, pointList):
    """Return the nearest index of points in array.
    Usage: getNearestPointIndex(a, pointList)"""
    if isinstance(pointList, tuple): pL = [pointList]
    else: pL = pointList

    if isinstance(a[0], list):
        # keep nearest
        npts = len(pL)
        res0 = [(0,1.e6) for i in range(npts)]
        noi = 0
        for i in a:
            res = geom.getNearestPointIndex(i, pL)
            for j in range(npts):
                if res0[j][1] > res[j][1]:
                    res0[j] = (res[j][0], res[j][1])
            noi += 1
        if isinstance(pointList, tuple): return res0[0]
        else: return res0
    else:
        res = geom.getNearestPointIndex(a, pL)
        if isinstance(pointList, tuple): return res[0]
        else: return res

def getCurvatureRadius(a):
    """Return the curvature radius for each point.
    Usage: getCurvatureRadius(a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(geom.getCurvatureRadius(i))
        return b
    else:
        return geom.getCurvatureRadius(a)

def getCurvatureAngle(a):
    """Return the curvature angle for each point.
    Usage: getCurvatureAngle(a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(geom.getCurvatureAngle(i))
        return b
    else:
        return geom.getCurvatureAngle(a)

def getCurvatureHeight(a):
    """Return the curvature height for each node in a 2D or 1D mesh.
    Usage: getCurvatureHeight(a)"""
    if isinstance(a[0], list):
        b = []
        for i in a: b.append(geom.getCurvatureHeight(i))
        return b
    else:
        return geom.getCurvatureHeight(a)

def getSharpestAngle(a):
    """Return the sharpest angle for each point of a surface based on the sharpest angle
    between adjacent element to which the point belongs to.
    Usage: getSharpestAngle(a)"""
    if isinstance(a[0], list):
        out = []
        for i in a: out.append(geom.getSharpestAngle(i))
        return out
    else: return geom.getSharpestAngle(a)

def getCurvilinearAbscissa(a):
    """Return the curvilinear abscissa for each point.
    Usage: getCurvilinearAbscissa(a)"""
    if isinstance(a[0], list):
        b = []
        for i in a: b.append(geom.getCurvilinearAbscissa(i))
        return b
    else:
        return geom.getCurvilinearAbscissa(a)

def getDistribution(a):
    """Return the curvilinear abscissa for each point as X coordinate.
    Usage: getDistribution(a)"""
    if isinstance(a[0], list):
        b = []
        for i in a:
            if i[-1] == 'BAR': raise TypeError("getDistribution: only for structured array.")
            c = line((0,0,0),(1,0,0),i[2])
            c[1][0] = geom.getCurvilinearAbscissa(i)[1]
            b.append(c)
        return b
    else:
        if a[-1] == 'BAR': raise TypeError("getDistribution: only for structured arrays.")
        c = line((0,0,0),(1,0,0),a[2])
        c[1][0] = geom.getCurvilinearAbscissa(a)[1]
        return c

def getTangent(a):
    """Return the unit tangent vector of a 1D curve as coordinates. 
    The input argument shall be an array. Each node of the output represents the unitary tangent 
    vector, pointing towards the tangent direction of the input 1D curve.
    Usage: getTangent(a)"""
    if not isinstance(a[0], list): Arrays = [a]
    else: Arrays = a
    b = []
    for a in Arrays:
        t = ['x,y,z',0,a[2],1,1]
        # Central difference
        n = a[1]
        OrientationAbsolute = 0.5*(numpy.diff(n[:,:-1],axis=1)+numpy.diff(n[:,1:],axis=1))
        # Not centered for the bounds
        OrientationAbsolute = numpy.hstack(((n[:,1]-n[:,0])[numpy.newaxis].T,
                                            OrientationAbsolute,
                                            (n[:,-1]-n[:,-2])[numpy.newaxis].T))
        #Norm = numpy.linalg.norm(OrientationAbsolute, axis=0)
        Norm = numpy.sqrt(numpy.sum(OrientationAbsolute*OrientationAbsolute, axis=0))
        OrientationRelative = OrientationAbsolute/Norm
        t[1] = OrientationRelative
        b.append(t)
    if len(b)==1: return b[0]
    else: return b

# Obsolete (use lineDrive)
def lineGenerate(a, d):
    return lineDrive(a, d)

def lineDrive(a, d):
    """Generate a surface mesh starting from a curve and a driving curve defined by d.
    Usage: lineDrive(a, d)"""
    if isinstance(d[0], list): # set of driving curves
        if isinstance(a[0], list):
            b = []
            for i in a:
                b.append(lineGenerate2__(i, d))
            return b
        else:
            return lineGenerate2__(a, d)
    else: # one driving curve
        if isinstance(a[0], list):
            b = []
            for i in a:
                b.append(geom.lineGenerate(i, d))
            return b
        else:
            return geom.lineGenerate(a, d)

def lineGenerate2__(array, drivingCurves):
    import Converter; import Generator
    # Copie la distribution de 0 sur les autres courbes
    d = []
    ref = drivingCurves[0]; d += [ref]
    #l = getLength(ref)
    distrib = getCurvilinearAbscissa(ref)
    distrib[0] = 'x'; distrib = Converter.addVars(distrib, ['y','z'])
    for i in drivingCurves[1:]:
        d += [Generator.map(i, distrib)]
    return geom.lineGenerate2(array, d)

# Ortho drive avec copy ou avec stack
# IN: a et d doivent etre orthogonals
# IN: mode=0 (stack), mode=1 (copy)
def orthoDrive(a, d, mode=0):
    """Generate a surface mesh starting from a curve and a driving orthogonally to curve defined by d.
    Usage: orthoDrive(a, d)"""
    try: import Generator as G; import Transform as T
    except ImportError:
        raise ImportError("orthoDrive: requires Generator module.")
    coord = d[1]
    center = (coord[0,0],coord[1,0],coord[2,0])
    coordA = a[1]
    P0 = [coordA[0,0],coordA[1,0],coordA[2,0]]
    P1 = [coordA[0,1],coordA[1,1],coordA[2,1]]
    xg = G.barycenter(a)
    v0 = Vector.sub(P0, xg)
    v1 = Vector.sub(P1, xg)
    S = Vector.cross(v0, v1)
    S = Vector.normalize(S)
    if abs(S[1]) > 1.e-12 and abs(S[2]) > 1.e-12: # x,S plane
        alpha = -S[0]
        U = Vector.mul(alpha, S)
        U = Vector.add(U, [1,0,0])
    else: # y,S plane
        alpha = -S[0]
        U = Vector.mul(alpha, S)
        U = Vector.add(U, [0,1,0])
    V = Vector.cross(U, S)

    n = d[2]
    all = []

    e2p = None
    P0 = [coord[0,0],coord[1,0],coord[2,0]]
    for i in range(n):
        if i == n-1:
            Pi = [coord[0,i-1],coord[1,i-1],coord[2,i-1]]
            Pip = [coord[0,i],coord[1,i],coord[2,i]]
            v = Vector.sub(Pip, P0)
        else:
            Pi = [coord[0,i],coord[1,i],coord[2,i]]
            Pip = [coord[0,i+1],coord[1,i+1],coord[2,i+1]]
            v = Vector.sub(Pi, P0)
        # vecteur e1 (transformation de S)
        e1 = Vector.sub(Pip, Pi)
        e1 = Vector.normalize(e1)
        # vecteur e2 (intersection plan)
        # intersection du plan normal a e1 avec le plan x,y
        if abs(S[1]) > 1.e-12 and abs(S[2]) > 1.e-12: # x,S plane
            alpha = -e1[0]
            e2 = Vector.mul(alpha, e1)
            e2 = Vector.add(e2, [1,0,0])
        else:
            alpha = -e1[0]
            e2 = Vector.mul(alpha, e1)
            e2 = Vector.add(e2, [0,2,0])
        e2 = Vector.normalize(e2)
        if e2p is not None:
            if Vector.dot(e2,e2p) < -0.9:
                e2 = Vector.mul(-1,e2)
        e2p = e2
        e3 = Vector.cross(e2,e1)

        b2 = T.rotate(a, center, (S,U,V), (e1,e2,e3))
        b2 = T.translate(b2, v)
        all.append(b2)
    if mode == 0: all = G.stack(all)
    return all

def addSeparationLine(array, array2):
    """Add a separation line defined in array2 to a mesh defined in array.
    Usage: addSeparationLine(array, array2)"""
    return geom.addSeparationLine(array, array2)

def axisym(a, C, axis, angle=360., Ntheta=180, rmod=None):
    """Create an axisymetrical mesh given an azimuthal 1D or 2D mesh.
    Usage: axisym(array, (xo,yo,zo), (nx,ny,nz), teta, Nteta, rmod)"""
    try:
        import Converter
        if rmod is not None: rmod = Converter.convertBAR2Struct(rmod)
    except: pass
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(geom.axisym(i, C, axis, angle, Ntheta, rmod))
        return b
    else:
        return geom.axisym(a, C, axis, angle, Ntheta, rmod)

def volumeFromCrossSections(contours):
    """Generate a 3D volume from cross sections contours in (x,y) planes.
    Usage: volumeFromCrossSections(contours)"""
    try:
        import KCore.kcore as KCore
        import Converter as C
        import Transform as T
        import Generator as G
    except ImportError:
        raise ImportError("volumeFromCrossSections: require Converter, Transform and Generator.")

    c = {}
    # Dictionnaire des contours suivant z
    for i in contours:
        posz = KCore.isNamePresent(i, "z")
        if posz == -1:
            posz = KCore.isNamePresent(i, "Z")
            if posz == -1:
                posz = KCore.isNamePresent(i, "CoordinateZ")
                if posz == 1:
                    raise TypeError("volumeFromCrossSections: Z coordinates not found in an array.")
        z = C.getValue(i, 0)[posz]
        if z in c:
            d = C.convertArray2Tetra(i)
            b = c[z]; f = T.join(d, b); c[z] = f
        else:
            d = C.convertArray2Tetra(i)
            c[z] = C.convertArray2Tetra(d)

    sort = sorted(c.items())

    # Delaunay constrained
    DT = []; CT = []
    for i in sort:
        d = G.close(i[1], 1.e-6); CT.append(d)
        #m = G.constrainedDelaunay(d, 1.e-10, 1)
        m = G.T3mesher2D(d); m = T.translate(m, (0,0,i[0]))
        DT.append(m)
    if len(sort) < 2:
        raise ValueError("volumeFromCrossSections: require at least two cross sections.")

    l = 0; vol = None
    for i in sort:
        if l == 0:
            vol = geom.volumeFromCrossSections(DT[l], DT[l+1], CT[l], CT[l+1])
        elif l < len(sort)-1:
            vol2 = geom.volumeFromCrossSections(DT[l], DT[l+1], CT[l], CT[l+1])
            vol = T.join(vol, vol2)
        l += 1
    return vol

# - text functions -
def text1D(string, font='vera', smooth=0, offset=0.5):
    """Create a 1D text.
    Usage: text1D(string, font, smooth, offset)"""
    if font == 'text1': from . import text1 as Text
    elif font == 'vera': from . import vera as Text
    elif font == 'chancery': from . import chancery as Text
    elif font == 'courier': from . import courier as Text
    elif font == 'nimbus': from . import nimbus as Text
    else: from . import text1 as Text
    try: import Transform
    except ImportError:
        raise ImportError("text1D: requires Transform.")
    retour = []
    offx = 0.; offy = 0.; s = 6
    for i in string:
        if i == 'A': a, s = Text.A()
        elif i == 'a': a, s = Text.a()
        elif i == 'B': a, s = Text.B()
        elif i == 'b': a, s = Text.b()
        elif i == 'C': a, s = Text.C()
        elif i == 'c': a, s = Text.c()
        elif i == 'D': a, s = Text.D()
        elif i == 'd': a, s = Text.d()
        elif i == 'E': a, s = Text.E()
        elif i == 'e': a, s = Text.e()
        elif i == 'F': a, s = Text.F()
        elif i == 'f': a, s = Text.f()
        elif i == 'G': a, s = Text.G()
        elif i == 'g': a, s = Text.g()
        elif i == 'H': a, s = Text.H()
        elif i == 'h': a, s = Text.h()
        elif i == 'I': a, s = Text.I()
        elif i == 'i': a, s = Text.i()
        elif i == 'J': a, s = Text.J()
        elif i == 'j': a, s = Text.j()
        elif i == 'K': a, s = Text.K()
        elif i == 'k': a, s = Text.k()
        elif i == 'L': a, s = Text.L()
        elif i == 'l': a, s = Text.l()
        elif i == 'M': a, s = Text.M()
        elif i == 'm': a, s = Text.m()
        elif i == 'N': a, s = Text.N()
        elif i == 'n': a, s = Text.n()
        elif i == 'O': a, s = Text.O()
        elif i == 'o': a, s = Text.o()
        elif i == 'P': a, s = Text.P()
        elif i == 'p': a, s = Text.p()
        elif i == 'Q': a, s = Text.Q()
        elif i == 'q': a, s = Text.q()
        elif i == 'R': a, s = Text.R()
        elif i == 'r': a, s = Text.r()
        elif i == 'S': a, s = Text.S()
        elif i == 's': a, s = Text.s()
        elif i == 'T': a, s = Text.T()
        elif i == 't': a, s = Text.t()
        elif i == 'U': a, s = Text.U()
        elif i == 'u': a, s = Text.u()
        elif i == 'V': a, s = Text.V()
        elif i == 'v': a, s = Text.v()
        elif i == 'W': a, s = Text.W()
        elif i == 'w': a, s = Text.w()
        elif i == 'X': a, s = Text.X()
        elif i == 'x': a, s = Text.x()
        elif i == 'Y': a, s = Text.Y()
        elif i == 'y': a, s = Text.y()
        elif i == 'Z': a, s = Text.Z()
        elif i == 'z': a, s = Text.z()
        elif i == '0': a, s = Text.C0()
        elif i == '1': a, s = Text.C1()
        elif i == '2': a, s = Text.C2()
        elif i == '3': a, s = Text.C3()
        elif i == '4': a, s = Text.C4()
        elif i == '5': a, s = Text.C5()
        elif i == '6': a, s = Text.C6()
        elif i == '7': a, s = Text.C7()
        elif i == '8': a, s = Text.C8()
        elif i == '9': a, s = Text.C9()
        elif i == '.': a, s = Text.POINT()
        elif i == ',': a, s = Text.COMMA()
        elif i == ';': a, s = Text.POINTCOMMA()
        elif i == ':': a, s = Text.TWOPOINTS()
        elif i == '!': a, s = Text.EXCLAMATION()
        elif i == '?': a, s = Text.INTERROGATION()
        elif i == '+': a, s = Text.PLUS()
        elif i == '-': a, s = Text.MINUS()
        elif i == '=': a, s = Text.EQUAL()
        elif i == '(': a, s = Text.LEFTBRACE()
        elif i == ')': a, s = Text.RIGHTBRACE()
        elif i == 'é': a, s = Text.EACUTE()
        elif i == 'è': a, s = Text.ELOW()
        elif i == 'à': a, s = Text.ALOW()
        elif i == 'ç': a, s = Text.CCEDILLE()
        elif i == '\n':
            offy = offy - 8 - offset
            offx = -6 - offset
            a = []
        else:
            a = []
        if a != []:
            a = Transform.translate(a, (offx,offy,0))
            retour += a
        offx += s + offset

    if smooth != 0:
        try: import Generator
        except ImportError:
            raise ImportError("text1D: requires Generator for smooth option.")
        if smooth == 1:
            nmap = 40; hdensify = 8./100
        elif smooth == 2:
            nmap = 40; hdensify = 8./10.
        elif smooth == 3:
            nmap = 40; hdensify = 8./5
        else:
            nmap = 40; hdensify = 8./2
        d = Generator.cart((0,0,0), (1./nmap,1,1), (nmap+1,1,1))
        c = 0
        for i in retour:
            b = Generator.densify(i, hdensify)
            retour[c] = Generator.map(b, d); c += 1

    return retour

def text2D(string, font='vera', smooth=0, offset=0.5):
    """Create a 2D text. 
    Usage: text2D(string, font, smooth, offset)"""
    try:
        import Generator
        import Transform
        import Converter
    except ImportError:
        raise ImportError("text2D: requires Generator, Transform, Converter.")
    a = text1D(string, font, smooth, offset)
    a = Converter.convertArray2Tetra(a)
    b = Transform.join(a)
    b = Generator.constrainedDelaunay(b)
    b = Generator.selectInsideElts(b, a)
    return b

def text3D(string, font='vera', smooth=0, offset=0.5, thickness=8.):
    """Create a 3D text.
    Usage: text3D(string, font, smooth, offset)"""
    try:
        import Transform
        import Converter
    except ImportError:
        raise ImportError("text3D: requires Generator, Transform, Converter.")
    a = text1D(string, font, smooth, offset)
    l = line((0,0,0),(0,0,thickness),2)
    a = lineDrive(a, l)
    a = Converter.convertArray2Tetra(a)
    b = Transform.join(a)

    a = text2D(string, font, smooth, offset)
    a = Transform.translate(a, (0,0,-0.0001))
    c = Transform.translate(a, (0,0,thickness+.0002))
    a = Transform.join([a, b, c])
    return a

#======================================================================
# connect 1D curves
# IN: sharpness=0: par des lignes, =1 par des splines
# IN: N: nbre de pts dans les raccords
# IN: lengthFactor: enleve les raccords trop longs
#======================================================================
def connect1D(curves, sharpness=0, N=10, lengthFactor=1.):
    """Connect 1D curves in a single curve.
    Usage: a = connect1D(A, sharpness, N, lengthFactor)"""
    try:
        import Transform as T
        import Converter as C
        import Generator as G
    except ImportError:
        raise ImportError("connect1D requires Transform, Converter, Generator modules.")
    #curves = T.splitTBranch(curves)
    curves = C.convertBAR2Struct(curves)
    ncurves = len(curves)

    Pts = []; PtsM= []
    for i in curves:
        ni = i[2]
        e1 = C.getValue(i, 0)
        e2 = C.getValue(i, ni-1)
        e1M = C.getValue(i, 1)
        e2M = C.getValue(i, ni-2)
        Pts.append([e1, e2])
        PtsM.append([e1M, e2M])

    added = []
    for c in range(ncurves):
        lcurve = getLength(curves[c]) * lengthFactor
        P1 = Pts[c][0]
        minDist, P2, d, ext = findNearest__(P1, Pts, c)
        n1 = Vector.sub(P1, PtsM[c][0])
        n1 = Vector.normalize(n1)
        n2 = Vector.sub(P2, PtsM[d][ext])
        n2 = Vector.normalize(n2)
        PI = intersectionPoint__(P1,n1,P2,n2)
        if sharpness == 0: # sharp
            la = line(P1,PI, N=N)
            lb = line(PI,P2, N=N)
            if getLength(la) < lcurve and getLength(lb) < lcurve:
                added += [la,lb]
        elif sharpness == 1: # spline
            controlPts = polyline([P1,PI,P2])
            sp = spline(controlPts, N=N)
            if getLength(sp) < lcurve: added += [sp]

        P1 = Pts[c][1]
        minDist, P2, d, ext = findNearest__(P1, Pts, c)
        n1 = Vector.sub(P1, PtsM[c][0])
        n1 = Vector.normalize(n1)
        n2 = Vector.sub(P2, PtsM[d][ext])
        n2 = Vector.normalize(n2)
        PI = intersectionPoint__(P1,n1,P2,n2)
        if sharpness == 0: # sharp
            la = line(P1,PI, N=N)
            lb = line(PI,P2, N=N)
            if getLength(la) < lcurve and getLength(lb) < lcurve:
                added += [la,lb]
        elif sharpness == 1: # spline
            controlPts = polyline([P1,PI,P2])
            sp = spline(controlPts, N=N)
            if getLength(sp) < lcurve: added += [sp]
    out = C.convertArray2Hexa(curves)
    if added != []:
        added = C.convertArray2Hexa(added)
        out = T.join(out+added)
    out = G.close(out)
    return out

# Pt d'intersection par minimal distance
def intersectionPoint__(P1,n1,P2,n2):
    s = Vector.dot(n1, n2)
    s2 = max(1.-s*s, 1.e-12)
    dP1P2 = Vector.sub(P1, P2)
    sn1 = Vector.mul(s, n1)
    p2 = Vector.sub(n2, sn1)
    q = Vector.dot(dP1P2, p2)
    tp = q / s2
    t = Vector.dot(Vector.sub(P2,P1),n1)+s*tp
    PI1 = Vector.add(P2, Vector.mul(tp,n2))
    PI2 = Vector.add(P1, Vector.mul(t, n1))
    PI = Vector.add(PI1,PI2)
    PI = Vector.mul(0.5,PI)
    return PI

# trouve le pt le plus proche de Pt dans Pts mais different de c
def findNearest__(Pt, Pts, c):
    minDist = 1.e6; nearest = None; dmin = -1; ext=0
    for d in range(len(Pts)):
        if d <= c: # possible sur lui meme !!
            e2a = Pts[d][0]; e2b = Pts[d][1]
            d1 = Vector.squareDist(Pt, e2a)
            d2 = Vector.squareDist(Pt, e2b)
            if d1 < minDist and d1 > 1.e-12:
                minDist = d1; nearest = e2a; ext=0; dmin = d
            if d2 < minDist and d2 > 1.e-12:
                minDist = d2; nearest = e2b; ext=1; dmin = d
    return minDist, nearest, dmin, ext

def getUV(a, normalDeviationWeight=2., texResolution=1920, fields=None):
    """Return uv of surface and atlas."""
    import Converter
    a = Converter.initVars(a, '_u_', 0.)
    a = Converter.initVars(a, '_v_', 0.)
    return geom.getUV(a, normalDeviationWeight, texResolution, fields)

# Init _u_ et _v_ from i,j (struct surface)
def getUVFromIJ(a):
    """Return uv of structured surface."""
    import Converter
    a = Converter.initVars(a, '_u_', 0.)
    a = Converter.initVars(a, '_v_', 0.)
    if isinstance(a[1], list): # array2/3
        pu = a[1][3].ravel('k'); pv = a[1][4].ravel('k')
    else: # array1
        pu = a[1][3].ravel('k'); pv = a[1][4].ravel('k')
    ni = a[2]; nj = a[3]
    for j in range(nj):
        for i in range(ni):
            pu[i+j*ni] = j*1./(nj-1)
            pv[i+j*ni] = i*1./(ni-1)
    return a
